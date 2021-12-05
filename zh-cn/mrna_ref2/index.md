# 转录组有参分析 （二)


## 介绍

本篇记录转录组有参分析主要基于STAR(比对)+featureCount(定量)。

## 1st, 前期准备
```shell
# 部署工作环境
mkdir -p ${HOME}/project/mRNA/{bin,data/{genes,genome,index,samples},output/{ballgown,hisat2,tmp}}

##bin，存放分析脚本
##data/genes，存放基因组gff
##data/genome,存放基因组fasta
##data/index，存放基因组索引文件index
##data/samples，存放样品数据reads fastq

# 准备数据
ln -s /path-to-reference-genome/ref.fa ${HOME}/project/mRNA/data/genome/
ln -s /path-to-reference-genome/ref.gff ${HOME}/project/mRNA/data/genes/

##/path-to-reference-genome/为基因组存放路径
##ref.fa为基因组fasta
##ref.gff3为基因组结构注释文件gff3

for i in `find /path-to-samples-fastq/ -type f -name "*.fq.gz"`;
do
	ln -s ${i} ${HOME}/project/mRNA/data/samples/;
done

##/path-to-samples-fastq/为样品fastq存放路径
```

## 2nd, 比对
```shell
cd ${HOME}/project/mRNA/bin
```
Copy the codes below to **02alignment.pbs**

```shell
#!/bin/bash -x

# 环境参数配置
NUMCPUS=28
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"
GENOMEIDX="${BASEDIR}/data/index/ref_tran"
GTFFILE="${BASEDIR}/data/genes/ref.gtf"
GENOMEFA="${BASEDIR}/data/genome/ref.fa"
WRKDIR="${BASEDIR}/output"
ALIGNLOC="${WRKDIR}/hisat2"
BALLGOWNLOC="${WRKDIR}/ballgown"
TEMPLOC="${WRKDIR}/tmp"

##NUMCPUS 使用CPU数，根据计算资源调整；
##THREADS，线程数为CPU数x2
##BASEDIR，proj基本目录
##FASTQLOC，reads fastq路径
##GENOMEIDX，基因组ref的index文件
##GTFFILE，基因组ref的结构注释文件gtf，可从gff文件转换
##GENOMEFA，基因组ref的fasta
##WRKDIR，proj输出目录
##ALIGNLOC，比对结果目录
##BALLGOWNLOC，定量等结果目录
##TEMPLOC，临时文件


# 将基因组注释文件gff3转为gtf

cd ${BASEDIR}/data/genes/
gffread ref.gff -T -o ${GTFFILE}

# 提取剪切位点和外显子信息
cd ${BASEDIR}/data/index
python2 /path-to-hisat_install_dir/extract_splice_sites.py ${GTFFILE} > ref.ss
python2 /path-to-hisat_install_dir/extract_exons.py ${GTFFILE} > ref.exon

##/path-to-hisat_install_dir/为HISAT2安装路径

# 建立基因组index
hisat2-build -p ${THREADS} --ss ref.ss --exon ref.exon ${GENOMEFA} ref_tran

# 提取reads信息
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")

###样品reads pair-end fastq必须是以"_1.*"或"_2.*"结尾，如"gill01_1.fastq.gz","gill01_2.fastq.gz"

# 开始比对

cd ${WRKDIR}

stime=`date +"%Y-%m-%d %H:%M:%S"`

echo "[$stime] #> START: HISAT2 Alignment"

for ((i=0; i<=${#reads1[@]}-1; i++ )); do
	sample="${reads1[$i]%%.*}"
	sample="${sample%_*}"
##根据reads1获得sample

	echo "[$stime] Processing sample: $sample"

	echo [$stime] "   * Alignment of reads to genome (HISAT2)"
	hisat2 -p $NUMCPUS \
		--dta \
		-x ${GENOMEIDX} \
		-1 ${FASTQLOC}/${reads1[$i]} \
		-2 ${FASTQLOC}/${reads2[$i]} \
		-S ${TEMPLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats

##-p，设置线程数
##--dta,
##--x，指定基因组index
##--1，输入reads1.fastq
##--2，输入reads2.fastq
##--S，指定sam

	echo [$stime] "   * Alignments conversion (SAMTools)"
	samtools view -S -b ${TEMPLOC}/${sample}.sam | \
		samtools sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam -
	rm ${TEMPLOC}/${sample}.sam

	echo [$stime] "   * Assemble transcripts (StringTie)"
	stringtie -p $NUMCPUS \
		-G ${GTFFILE} \
		-o ${ALIGNLOC}/${sample}.gtf \
		-l ${sample} ${ALIGNLOC}/${sample}.bam

done

echo [$stime] "#> Step Two : Alignment Finished!"
```

## 3rd, 转录本定量

```shell
#!/bin/bash -x

# 环境参数配置
NUMCPUS=28
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"
GENOMEIDX="${BASEDIR}/data/index/ref_tran"
GTFFILE="${BASEDIR}/data/genes/ref.gtf"
GENOMEFA="${BASEDIR}/data/genome/ref.fa"
WRKDIR="${BASEDIR}/output"
ALIGNLOC="${WRKDIR}/hisat2"
BALLGOWNLOC="${WRKDIR}/ballgown"
TEMPLOC="${WRKDIR}/tmp"

stime=`date +"%Y-%m-%d %H:%M:%S"`

echo [$stime] "#> START: STRINGTIE Estimation"

echo [$stime] "#> Merge all transcripts (StringTie)"

ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt

stringtie --merge \
	-p $NUMCPUS \
	-G  ${GTFFILE} \
	-o ${BALLGOWNLOC}/stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt
##
##
##
##

echo [$stime] "#> Estimate abundance for each sample (StringTie)"

for ((i=0; i<=${#reads1[@]}-1; i++ )); do

	sample="${reads1[$i]%%.*}"
	sample="${sample%_*}"

	if [ ! -d ${BALLGOWNLOC}/${sample} ]; then
		mkdir -p ${BALLGOWNLOC}/${sample}
	fi

# Estimate transcript abundance depend on the referenge genome with ignoring novel transcripts
	stringtie -e -B \
		-p $NUMCPUS \
		-G ${GTFFILE} \
		-o ${BALLGOWNLOC}/${sample}/${sample}.gtf ${ALIGNLOC}/${sample}.bam

##
##
##

done

# Generate transcript_count_matrix.csv and gene_count_matrix.csv using the script under stringtie_install_dir

python2 /public/tools/rna_seq/stringtie-1.3.4/prepDE.py -l 150

echo [$stime] "#> Step Three: Estimation Finished!"
```

后续差异基因分析和功能分析，参考RNA-seq analysis (有参分析 一)

## Reference
-	https://www.nature.com/articles/nprot.2016.095
-	http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
-	https://guangchuangyu.github.io/software/clusterProfiler/


