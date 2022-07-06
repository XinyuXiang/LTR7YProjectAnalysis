### Software summary
<< EOF
  1. Fastqc -  v0.11.8
  2. Cutadapt - v2.9
  3. STAR - v2.7.0e
  4. Samtools - v1.2
  5. bamCoverage - v3.3.1
  6. featureCounts - v2.0.0
EOF

### Reference genome
<< EOF
Reference genome: Homo_sapiens.GRCh38.97.dna.primary_assembly.fa
Gene annotation: Homo_sapiens.GRCh38.97.gtf
TE annotation: GRCh38_rmsk_TE.gtf
EOF


### Quality control
fastqc ${seqfile}/${seqfile}_1.fq.gz ${seqfile}/${seqfile}_2.fq.gz -o ../fastqc_result/

### Trim data 
cutadapt -u 10 -u -90 -U 10 -U -90 -o ${seqfile}_1_trimmed.fq.gz -p ${seqfile}_2_trimmed.fq.gz ${seqfile}_1.fq.gz ${seqfile}_2.fq.gz 
# or cutadapt -u 10 -U 10 -o ${seqfile}_1_trimmed.fq.gz -p ${seqfile}_2_trimmed.fq.gz ${seqfile}_1.fq.gz ${seqfile}_2.fq.gz

### STAR alignment 
# STAR   --runMode genomeGenerate   --runThreadN 8   --genomeDir /home/zje610/workspace/refData/ensemble/STAR_index_GRCh38   --genomeFastaFiles /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa      --sjdbGTFfile /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.gtf

STAR --runThreadN 20 \
--genomeDir /home/zje610/workspace/refData/ensemble/STAR_index_GRCh38/ \
--readFilesCommand zcat \
--readFilesIn ./fq/${seqfile}_trimmed_read1.fastq.gz ./fq/${seqfile}_trimmed_read2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1000 \
--outFilterMismatchNmax 3 \
--outSAMmultNmax 1 \
--outFileNamePrefix $seqfile

samtools index ${seqfile}Aligned.sortedByCoord.out.bam


### Bam to bigwig
bamCoverage -b ${seqfile}Aligned.sortedByCoord.out.bam -bs 10 --normalizeUsing RPKM -o ${seqfile}Aligned.sortedByCoord.bw

### Count
# paired end library (stranded) 
# count gene
# awk '$3=="gene"' Homo_sapiens.GRCh38.97.gtf | grep protein_coding > Homo_sapiens.GRCh38.97.PCG.gtf
featureCounts -a /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.PCG.gtf \
-s 2 -p -B -P -C -T 20 \
-t gene -g gene_id \
-o ${seqfile}_gene.featureCounts_PCG.txt ${seqfile}Aligned.sortedByCoord.out.bam
# count TE
featureCounts -a /home/zje610/workspace/refData/repeatmasker/GRCh38_rmsk_TE.gtf \
-s 2 -p -B -P -C -M -T 20 \
-t exon -g transcript_id \
-o ${seqfile}_TE.featureCounts_exon_TEcy.txt ${seqfile}Aligned.sortedByCoord.out.bam


# paired end library (unstranded)
# count gene
featureCounts -a /home/zje610/workspace/refData/ensemble/Homo_sapiens.GRCh38.97.PCG.gtf \
-p -B -P -C -T 20 \
-t gene -g gene_id \
-o ${seqfile}_gene.featureCounts_PCG.txt ${seqfile}Aligned.sortedByCoord.out.bam
# count TE
featureCounts -a /home/zje610/workspace/refData/repeatmasker/GRCh38_rmsk_TE.gtf \
-p -B -P -C -M -T 20 \
-t exon -g transcript_id \
-o ${seqfile}_TE.featureCounts_exon_TEcy.txt ${seqfile}Aligned.sortedByCoord.out.bam





