# 转录组有参分析（二)


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

## 2nd, 第一次比对
```shell
#!/bin/bash -x

# 环境参数配置
NUMCPUS=28
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"

GTFFILE="${BASEDIR}/data/genes/ref.gtf"
GENOMEFA="${BASEDIR}/data/genome/ref.fa"
GENOMEIDXDIR="$BASEDIR/data/index"
GENOMEIDX_TPDIR="$BASEDIR/data/index_twopass"
WRKDIR="${BASEDIR}/output"
INDEX_TP=`ls ${WORKDIR}/*SJ.out.tab`
LOGFILE=${WRKDIR}/star_align.log


##NUMCPUS 使用CPU数，根据计算资源调整；
##THREADS，线程数为CPU数x2
##BASEDIR，proj基本目录
##FASTQLOC，reads fastq路径
##GTFFILE，基因组ref的结构注释文件gtf，可从gff文件转换
##GENOMEFA，基因组ref的fasta
##GENOMEIDXDIR，基因组ref的index路径
##GENOMEIDX_TPDIR，基因组ref的index路径
##WRKDIR，proj输出目录
##INDEX_TP
##LOGFILE

stime=`date +"%Y-%m-%d %H:%M:%S"`


set -e
#set -x
SCRIPTARGS="$@"


# 将基因组注释文件gff3转为gtf

cd ${BASEDIR}/data/genes/
gffread ref.gff -T -o ${GTFFILE}


# 建立基因组index
cd $WRKDIR

first_index(){
    echo [$stime] "#> START: " $0 $SCRIPTARGS
    echo [$stime] "   * Generating genome indexes (STAR)"

    STAR --runThreadN ${NUMCPUS} \
    --runMode genomeGenerate \
    --genomeDir ${GENOMEIDX} \
    --genomeFastaFiles ${GENOMEFA} \
    --sjdbGTFfile ${GTFFILE} \
    --sjdbOverhang 149

##
##
##
##
##

    echo [$stime] "#> Genome Index DONE."
}

first_index 2>&1 | tee ${LOGFILE}

# 提取reads信息
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")

###样品reads pair-end fastq必须是以"_1.*"或"_2.*"结尾，如"gill01_1.fastq.gz","gill01_2.fastq.gz"

# 开始比对

cd ${WRKDIR}

first_align() {
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"

    echo "[$stime] Processing sample: $sample"
    echo [$stime] "   * Running mapping jobs (STAR)"

    STAR --runThreadN ${NUMCPUS} \
    	--genomeDir ${GENOMEIDX} \
    	--readFilesIn ${FASTQLOC}/${reads1[$i]} ${FASTQLOC}/${reads2[$i]} \
	 	--readFilesCommand zcat \
	 	--outSAMtype BAM Unsorted \
	 	--outFileNamePrefix ${sample}

    echo [$stime] "   * $sample Finished "
done

echo [$stime]" 1st mapping jobs (STAR) Finished."
}

first_align 2>&1 | tee $LOGFILE

```

## 3rd, 第二次比对

```shell
#!/bin/bash -x

# 环境参数配置
NUMCPUS=28
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"

GTFFILE="${BASEDIR}/data/genes/ref.gtf"
GENOMEFA="${BASEDIR}/data/genome/ref.fa"
GENOMEIDXDIR="$BASEDIR/data/index"
GENOMEIDX_TPDIR="$BASEDIR/data/index_twopass"
WRKDIR="${BASEDIR}/output"
INDEX_TP=`ls ${WORKDIR}/*SJ.out.tab`
LOGFILE=${WRKDIR}/star_align.log

stime=`date +"%Y-%m-%d %H:%M:%S"`


two-pass_index(){
    echo $stime] "#> START: " ${0} ${SCRIPTARGS}
    echo [$stime] "   * Generating 2-pass genome indexes (STAR)"
    
    STAR --runThreadN ${NUMCPUS} \
        --runMode genomeGenerate \
        --genomeDir ${GENOMEIDX_TP} \
        --genomeFastaFiles ${GENOMEFA} \
        --sjdbGTFfile ${GTFFILE} \
        --sjdbFileChrStartEnd ${INDEX2PASS} \
        --sjdbOverhang 149

	echo [$stime] "#> 2-pass index DONE."

}

two-pass_index 2>&1 | tee $LOGFILE


two-pass_align() {

for ((i=0; i<=${#reads1[@]}-1; i++ )); do
    sample="${reads1[$i]%%.*}"
    sample="${sample%_*}"

    echo "[$stime] Processing sample: $sample"
    echo [$stime] "   * Running mapping jobs (STAR)"
    
    STAR --runThreadN ${NUMCPUS} \
        --genomeDir ${GENOMEIDX_TP} \
        --readFilesIn ${FASTQLOC}/${reads1[$i]} ${FASTQLOC}/${reads2[$i]} \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ./align_2pass/${sample}_2-pass \
    echo [$stime] "   * $sample Finished "
done

echo [$stime] " 2-pass mode mapping jobs (STAR) Finished."
}

two-pass_align 2>&1 | tee $LOGFILE
```

## 4th, 转录本定量

```shell

```

后续差异基因分析和功能分析，参考RNA-seq analysis (有参分析 一)

## Reference



