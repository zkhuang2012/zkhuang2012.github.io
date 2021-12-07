# 转录组有参分析（一)


## 介绍

常规转录组分析：**数据指控（QC）-> 比对（Alignment）-> 定量（Quantity）-> 差异基因分析（DEG）-> 功能富集（Enrichment）**

本篇记录的转录组有参分析流程主要包括：
- HISAT2（比对）
- Stringtie（定量）
- DESeq2（差异分析）
- clusterProfiler（富集分析）

## 1st, 前期准备
```shell
# 部署工作环境
mkdir -p ${HOME}/project/mRNA/{bin,data/{genes,genome,index,samples},output/{hisat2_dir,stringtie_dir,tmp}}

##bin，存放分析脚本
##data/genes，存放基因组gff
##data/genome,存放基因组fasta
##data/index，存放基因组索引文件index
##data/samples，存放样品数据reads fastq
##output/hisat2_dir,存放比对结果，如bam
##output/stringtie_dir，存放定量结果等

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
ALIGNLOC="${WRKDIR}/hisat2_dir"
STRINGTIELOC="${WRKDIR}/stringtie_dir"
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
##STRINGTIELOC，定量等结果目录
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

##-p，
##-G，
##-o，
##-l，


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
ALIGNLOC="${WRKDIR}/hisat2_dir"
STRINGTIELOC="${WRKDIR}/stringtie_dir"
TEMPLOC="${WRKDIR}/tmp"

stime=`date +"%Y-%m-%d %H:%M:%S"`

echo [$stime] "#> START: STRINGTIE Estimation"

# 转录本重构
echo [$stime] "#> Merge all transcripts (StringTie)"

ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt

stringtie --merge \
	-p $NUMCPUS \
	-G  ${GTFFILE} \
	-o ${STRINGTIELOC}/stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt

##--merge，Transcript merge mode,StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. 
##-G，基因组ref的gtf文件,StringTie will assemble the transfrags from the input GTF files with the reference transcripts.
##-o，指定输出样品gtf文件

# 转录本定量
echo [$stime] "#> Estimate abundance for each sample (StringTie)"

for ((i=0; i<=${#reads1[@]}-1; i++ )); do

	sample="${reads1[$i]%%.*}"
	sample="${sample%_*}"

	if [ ! -d ${STRINGTIELOC}/${sample} ]; then
		mkdir -p ${STRINGTIELOC}/${sample}
	fi

	stringtie -e -B \
		-p $NUMCPUS \
		-G ${GTFFILE} \
		-o ${STRINGTIELOC}/${sample}/${sample}.gtf ${ALIGNLOC}/${sample}.bam

##-e，expression estimation mode; this limits the processing of read alignments to estimating the coverage of the transcripts given with the -G option
##-B，enables the output of Ballgown input table files (*.ctab) containing coverage data for the reference transcripts given with the -G option. 


done

# 提取read count矩阵

python2 /public/tools/rna_seq/stringtie-1.3.4/prepDE.py -l 150

echo [$stime] "#> Step Three: Estimation Finished!"
```


## 4th, 差异基因分析
```R
library("DESeq2")

path=getwd()
countData <- as.matrix(read.csv(paste(path,"/transcript_count_matrix.csv",sep=""), row.names="transcript_id"))
colData <- read.csv(paste(path,"/phenodata.csv",sep=""), sep=",", row.names=1)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)
dds <- DESeq(dds)

# filter count<10
dds_filter <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds_filter)

# plot PCA
rld <- rlog(dds)
pdf("pca_t.pdf")
plotPCA(rld, intgroup="group")
dev.off()

# take group H_Red and H_YX for examples
res<- results(dds,contrast=c("group","H_Red","C_Red"), alpha=0.05)
resOrdered <- res[order(res$padj), ]
diff_gene_deseq2 <-subset(resOrdered,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- merge(as.data.frame(diff_gene_deseq2), as.data.frame(counts(dds,normalize=TRUE)), by="row.names",sort=FALSE)
write.csv(diff_gene_deseq2, file="./Red_HvC.sig.csv", row.names = F)

EOF

```

## 5th, 功能富集
```R

cat <<EOF >>enrichment.R

sig_gene <- read.csv("./Red_HvC.sig.csv",header = T)
universe <- read.delim("/path-to-reference-genome/ref_universe.txt",sep="\t",header=FALSE,stringsAsFactors = F)
go2gene <- read.delim("/path-to-reference-genome/ref_go2gene.txt", sep="\t",header = FALSE,stringsAsFactors = F)
go2name <- read.delim("/path-to-reference-genome/ref_go2term.txt", sep="\t",header = FALSE,stringsAsFactors = F)
map2gene <- read.delim("/path-to-reference-genome/ref_map2gene.txt", sep="\t",header = FALSE,stringsAsFactors = F)
map2name <- read.delim("/path-to-reference-genome/ref_map2term.txt", sep="\t",header = FALSE,stringsAsFactors = F)

library("clusterProfiler")

# ge enrichment
go=enricher(gene = sig_gene[,1],
	universe = universe[,1],
	TERM2GENE= go2gene,
	TERM2NAME = go2name[,1:2])

pdf("go_enrich.pdf")
barplot(go,showcatergory=20)
#View(data.frame(go))
dev.off()

# kegg enrichment
kegg=enricher(gene = sig_gene[,1],
	universe = universe[,1],
	TERM2GENE= map2gene,
	TERM2NAME = map2name)

pdf("kegg_enrich.pdf")
barplot(kegg,showcatergory=20)
dev.off()

EOF
```

## Reference
-	https://www.nature.com/articles/nprot.2016.095
-	http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
-	https://guangchuangyu.github.io/software/clusterProfiler/


