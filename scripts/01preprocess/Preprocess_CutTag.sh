### Software summary
<< EOF
  1. Fastqc -  v0.11.8
  2. Cutadapt - v2.9
  3. STAR - v2.7.0e
  4. Samtools - v1.2
  5. bamCoverage - v3.3.1
  6. macs2 - v2.2.7.1
  7. homer - v4.7
EOF

### Reference genome
<< EOF
Reference genome: Homo_sapiens.GRCh38.97.dna.primary_assembly.fa
EOF


### Quality control
fastqc ${seqfile}/${seqfile}_1.fq.gz ${seqfile}/${seqfile}_2.fq.gz -o ../fastqc_result/

### Trim data 
cutadapt -m 30 -j 8 -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG -o ${seqfile}_1_trimmed.fq.gz -p ${seqfile}_2_trimmed.fq.gz ${seqfile}_1.fq.gz ${seqfile}_2.fq.gz

### STAR alignment 
# STAR   --runMode genomeGenerate   --genomeDir STAR_index_GRCh38_withoutGTF/   --genomeFastaFiles Homo_sapiens.GRCh38.97.dna.primary_assembly.fa   

STAR --runThreadN 20 \
--genomeDir /home/zje610/workspace/refData/ensemble/STAR_index_GRCh38_withoutGTF/ \
--readFilesCommand zcat \
--readFilesIn ./fq/${seqfile}_1_trimmed.fq.gz ./fq/${seqfile}_2_trimmed.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1000 \
--outFilterMismatchNmax 3 \
--alignIntronMax 1 \
--outSAMmultNmax 1 \
--outFileNamePrefix $seqfile

samtools rmdup ${seqfile}Aligned.sortedByCoord.out.bam ${seqfile}Aligned.sortedByCoord_rmdup.bam
samtools index ${seqfile}Aligned.sortedByCoord_rmdup.bam


### Bam to bigwig
bamCoverage -b ${seqfile}Aligned.sortedByCoord_rmdup.bam -bs 10 --normalizeUsing RPKM -o ${seqfile}Aligned.sortedByCoord_rmdup.bw

### Call peak
macs2 callpeak -t ${seqfile}Aligned.sortedByCoord_rmdup.bam -f BAM -g hs -n ${seqfile}Aligned.sortedByCoord --call-summits

### Motif
sed 's/^/chr&/g' ${seqfile}_summits.bed > tmp
mv tmp ${seqfile}_summits.bed
    
annotatePeaks.pl ${seqfile}_summits.bed hg38 -size given -annStats ${seqfile}_annStats.txt > ${seqfile}_Annotation.txt
findMotifsGenome.pl ${seqfile}_summits.bed hg38 ${seqfile}_motif -size 200 -mask




