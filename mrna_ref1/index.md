# RNA-seq Ref-align Analysis (1)


## Introduction

This RNA-seq pipeline was built with HISAT2 + Stringtie + DESeq2 + clusterProfiler.

## 1st, Prepariation
```shell
# Directory Deployment
mkdir -p ${HOME}/project/mRNA/{bin,data/{genes,genome,index,samples},output/{ballgown,hisat2,tmp}}


# Data link
cd ${HOME}/../Data_benthic/practice/mRNA
ln -s ${PWD}/Hdh.fa ${HOME}/project/mRNA/data/genome/
ln -s ${PWD}/Hdh.gff3 ${HOME}/project/mRNA/data/genes/

for i in `find ${PWD} -type f -name "*.fq.gz"`;
do
	ln -s ${i} ${HOME}/project/mRNA/data/samples/;
done

```

## 2nd, Alignment
```shell
cd ${HOME}/project/mRNA/bin
```
Copy the codes below to **02alignment.pbs**

```shell
#!/bin/bash -x
#PBS -N ref_mRNA_02alignment
#PBS -o ref_mRNA_02alignment.out
#PBS -e ref_mRNA_02alignment.err
#PBS -q high
#PBS -j oe
#PBS -l nodes=1:ppn=28

NUMCPUS=`wc -l <$PBS_NODEFILE`
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"
GENOMEIDX="${BASEDIR}/data/index/Hdh_tran"
GTFFILE="${BASEDIR}/data/genes/Hdh.gtf"
GENOMEFA="${BASEDIR}/data/genome/Hdh.fa"
WRKDIR="${BASEDIR}/output"
ALIGNLOC="${WRKDIR}/hisat2"
BALLGOWNLOC="${WRKDIR}/ballgown"
TEMPLOC="${WRKDIR}/tmp"

# transfer gff to gtf
cd ${BASEDIR}/data/genes/
gffread Hdh.gff -T -o ${GTFFILE}

# Extract splice sites and exon info from gtf using the script under hisat_install_dir
cd ${BASEDIR}/data/index
python2 /public/tools/rna_seq/hisat2-2.1.0/extract_splice_sites.py ${GTFFILE} > Hdh.ss
python2 /public/tools/rna_seq/hisat2-2.1.0/extract_exons.py ${GTFFILE} > Hdh.exon

# Build genome index for hisat2
hisat2-build -p 50 --ss Hdh.ss --exon Hdh.exon ${GENOMEFA} Hdh_tran

## sample pair-end fastq must be named in this format “_1.*”"_2.*". eg, "gill01_1.fastq.gz","gill01_2.fastq.gz"
reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")


cd ${WRKDIR}

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: HISAT2 Alignment"
for ((i=0; i<=${#reads1[@]}-1; i++ )); do
	sample="${reads1[$i]%%.*}"
	sample="${sample%_*}"

	stime=`date +"%Y-%m-%d %H:%M:%S"`
	echo "[$stime] Processing sample: $sample"
	echo [$stime] "   * Alignment of reads to genome (HISAT2)"

	hisat2 -p $NUMCPUS \
		--dta \
		-x ${GENOMEIDX} \
		-1 ${FASTQLOC}/${reads1[$i]} \
		-2 ${FASTQLOC}/${reads2[$i]} \
		-S ${TEMPLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
	samtools view -S -b ${TEMPLOC}/${sample}.sam | \
		samtools sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam -

	rm ${TEMPLOC}/${sample}.sam

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
	stringtie -p $NUMCPUS \
		-G ${GTFFILE} \
		-o ${ALIGNLOC}/${sample}.gtf \
		-l ${sample} ${ALIGNLOC}/${sample}.bam

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Step Two : Alignment Finished!"
```

Then run via PBS:

```shell
qsub 02alignment.pbs

```

## 3rd, Estimate transcript abundance

Copy the codes below to **03estimation.pbs**

```shell
#!/bin/bash -x
#PBS -N ref_mRNA_03estimation
#PBS -o ref_mRNA_03estimation.out
#PBS -e ref_mRNA_03estimation.err
#PBS -q high
#PBS -j oe
#PBS -l nodes=1:ppn=28

NUMCPUS=`wc -l <$PBS_NODEFILE`
THREADS=`expr ${NUMCPUS} \* 2`
BASEDIR="${HOME}/project/mRNA"
FASTQLOC="${BASEDIR}/data/samples"
GENOMEIDX="${BASEDIR}/data/index/Hdh_tran"
GTFFILE="${BASEDIR}/data/genes/Hdh.gtf"
GENOMEFA="${BASEDIR}/data/genome/Hdh.fa"
WRKDIR="${BASEDIR}/output"
ALIGNLOC="${WRKDIR}/hisat2"
BALLGOWNLOC="${WRKDIR}/ballgown"
TEMPLOC="${WRKDIR}/tmp"

reads1=(${FASTQLOC}/*_1.*)
reads1=("${reads1[@]##*/}")
reads2=("${reads1[@]/_1./_2.}")

cd ${WRKDIR}

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: STRINGTIE Estimation"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge all transcripts (StringTie)"
ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt

stringtie --merge \
	-p $NUMCPUS \
	-G  ${GTFFILE} \
	-o ${BALLGOWNLOC}/stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt

## Estimate transcript abundance
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Estimate abundance for each sample (StringTie)"
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

done

# Generate transcript_count_matrix.csv and gene_count_matrix.csv using the script under stringtie_install_dir

python2 /public/tools/rna_seq/stringtie-1.3.4/prepDE.py -l 150

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Step Three: Estimation Finished!"
```

Then run via PBS:

```shell
qsub 03estimation.pbs

```


## 4th, Differential expression analysis
```R

ln -s /public/home/benthic/Data_benthic/practice/mRNA/phenodata.csv ${HOME}/project/mRNA/output/phenodata.csv

cd ${HOME}/project/mRNA/output

cat <<EOF >>deseq.r

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

qsub -I -N DESeq2 -l nodes=1:ppn=2 -l walltime=1000:00:00 -q high `Rscript deseq.R`

```

## 5th, Funtional enrichment analysis
```R

cat <<EOF >>enrichment.R

sig_gene <- read.csv("./Red_HvC.sig.csv",header = T)
universe <- read.delim("/public/home/benthic/Data_benthic/practice/mRNA/Hdh_universe.txt",sep="\t",header=FALSE,stringsAsFactors = F)
go2gene <- read.delim("/public/home/benthic/Data_benthic/practice/mRNA/Hdh_go2gene.txt", sep="\t",header = FALSE,stringsAsFactors = F)
go2name <- read.delim("/public/home/benthic/Data_benthic/practice/mRNA/Hdh_go2term.txt", sep="\t",header = FALSE,stringsAsFactors = F)
map2gene <- read.delim("/public/home/benthic/Data_benthic/practice/mRNA/Hdh_map2gene.txt", sep="\t",header = FALSE,stringsAsFactors = F)
map2name <- read.delim("/public/home/benthic/Data_benthic/practice/mRNA/Hdh_map2term.txt", sep="\t",header = FALSE,stringsAsFactors = F)

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

qsub -I -N enrichment -l nodes=1:ppn=2 -l walltime=1000:00:00 -q high `Rscript enrichment.R`

```

## Reference
-	https://www.nature.com/articles/nprot.2016.095
-	http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
-	https://guangchuangyu.github.io/software/clusterProfiler/


