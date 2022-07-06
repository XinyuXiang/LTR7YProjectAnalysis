### Software summary
<< EOF
  1. cellranger -  v3.1.0
  2. Fastqc -  v0.11.8
  3. Cutadapt - v2.9
  4. STAR - v2.7.0e
  5. Samtools - v1.2
  6. bamCoverage - v3.3.1
  7. featureCounts - v2.0.0
  8. HTseq - v1.99.2

EOF

### Reference genome
<< EOF
Reference genome: Homo_sapiens.GRCh38.97.dna.primary_assembly.fa
Gene annotation: Homo_sapiens.GRCh38.97.gtf
TE annotation: GRCh38_rmsk_TE.gtf
EOF


### Cellranger for GSE207431 
# Make reference
cellranger mkref \
--genome=refdata-cellranger-GRCh38-gene-3.1.0 \
--fasta=/home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa \
--genes=/home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.gtf

cellranger mkref \
--genome=refdata-cellranger-GRCh38-TE-3.1.0 \
--fasta=/home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa \
--genes=/home/zje610/workspace/refData/repeatmasker/GRCh38_rmsk_TE.gtf

# Align & Count
/home/zje610/biosoft/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger count \
--id=scFINEmerge2_gene \
--transcriptome=/home/zje610/workspace/refData/cellranger/refdata-cellranger-GRCh38-gene-3.1.0 \
--fastqs=/home/zje610/workspace/labW/XXY/FINE_scRNA \
--localcores 48 \
--chemistry=SC3Pv3 \
--sample=scFINE_gene

/home/zje610/biosoft/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger count \
--id=scFINEmerge2_TE \
--transcriptome=/home/zje610/workspace/refData/cellranger/refdata-cellranger-GRCh38-TE-3.1.0 \
--fastqs=/home/zje610/workspace/labW/XXY/FINE_scRNA \
--localcores 48 \
--chemistry=SC3Pv3 \
--sample=scFINE_TE





### STAR + featureCounts for GSE36552
# Trim data 
cutadapt -u 10 -u -40 -o ${seqfile}_trimmed.fq.gz ${seqfile}.fq.gz

# Align
# STAR   --runMode genomeGenerate   --runThreadN 8   --genomeDir /home/zje610/workspace/refData/ensemble/STAR_index_GRCh38   --genomeFastaFiles /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa      --sjdbGTFfile /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.gtf
STAR --runThreadN 20 \
--genomeDir /home/zje610/workspace/refData/ensemble/STAR_index_GRCh38/ \
--readFilesCommand zcat \
--readFilesIn ./fq/${seqfile}_trimmed.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1000 \
--outSAMmultNmax 1 \
--outFileNamePrefix ${seqfile}

samtools index ${seqfile}Aligned.sortedByCoord.out.bam

# Count
featureCounts -a ~/scratch/refData/Homo_sapiens.GRCh38.97.PCG.gtf \
-B -P -C -T 64 \
-t gene -g gene_id \
-o ${seqfile}.featureCounts_PCG.txt ${seqfile}Aligned.sortedByCoord.out.bam

featureCounts -a ~/scratch/refData/GRCh38_rmsk_TE.gtf \
-B -P -C -M -T 64 \
-t exon -g gene_id \
-o ${seqfile}.featureCounts_exon_TEcy.txt ${seqfile}Aligned.sortedByCoord.out.bam

# Bam to bigwig
bamCoverage -b ${seqfile}Aligned.sortedByCoord.out.bam -bs 10 --normalizeUsing RPKM -o ${seqfile}Aligned.sortedByCoord.bw





### HTseq for GSE109555
# Dump SRA and split fq
fastq-dump --gzip --split-3 ${seqfile}.sra
mkdir ${seqfile}
perl s01.Barcode_UMI_QC_per1w_V2.pl ${seqfile}_1.fastq.gz ${seqfile}_2.fastq.gz ../GSE109555_RAW/ ./${seqfile}

# Quality control
fastqc ${file%%.*}.R1.clean.fq.gz -o ./

# Align
STAR --runThreadN 64 --genomeDir /u/home/y/yutaomcd/scratch/refdata/STAR_index_GRCh38/ \
--readFilesCommand zcat \
--readFilesIn ${file%%.*}.R1.clean.fq.gz \
--outFilterMultimapNmax 1000 --outFilterMismatchNmax 3 --outSAMmultNmax 1 \
--outFileNamePrefix ${file%%.*}

samtools view -bS ${file%%.*}\Aligned.out.sam > ${file%%.*}\.bam
samtools sort ${file%%.*}\.bam -o ${file%%.*}\_srt.bam
samtools index ${file%%.*}\_srt.bam

# Bam to bigwig
bamCoverage -b ${file%%.*}\_srt.bam -bs 10 --normalizeUsing RPKM -o ${file%%.*}\.bw

rm ${file%%.*}\Aligned.out.sam;
rm ${file%%.*}\.bam

# Count by HTseq and collapse UMI
count.py \
    -s no -f bam -a 10 \
    --samout-format=sam \
    -t exon -i gene_name \
    ${file_head}_srt.bam  \
    /u/home/y/yutaomcd/scratch/refData/Homo_sapiens.GRCh38.97.gtf \
    -o ${file_head}_gene.sam \
    > ${file_head}_gene.txt  

python3.7 UMI_HTseq.py ${file_head}_gene.sam  ${file_head}_gene_umi.txt # UMI_HTseq.py was from https://github.com/WRui/Post_Implantation/



count.py \
    -s no -f bam -a 10 \
    --samout-format=sam \
    -t exon -i gene_id \
    ${file_head}_srt.bam  \
    /u/home/y/yutaomcd/scratch/refData/GRCh38_rmsk_TE.gtf \
    -o ${file_head}_TE.sam \
    > ${file_head}_TE.txt  
   
python3.7 UMI_HTseq.py ${file_head}_TE.sam ${file_head}_TE_umi.txt # UMI_HTseq.py was from https://github.com/WRui/Post_Implantation/










