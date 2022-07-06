### Software summary
<< EOF
  1. deeptools - v3.2.1 
EOF


### plot heatmap over TE regions
computeMatrix scale-regions \
-S bigwigfile1.bw bigwigfile2.bw \
-R TEregion1.bed TEregion2.bed \
-p 20 \
--beforeRegionStartLength 5000 --regionBodyLength 1000 --afterRegionStartLength 5000 --skipZeros \
--outFileSortedRegions TEregion_ordered.bed \
--sortRegions descend \
-o heatmap_TE.mat.gz

plotHeatmap -m heatmap_TE.mat.gz \
-out heatmap_TE.pdf \
--sortRegions keep \
--heatmapHeight 10 \
--boxAroundHeatmaps no \
--colorMap Blues \
--missingDataColor "#FFFAF4"

### motif enrichment
# obtain motif region
scanMotifGenomeWide.pl nanog.motif hg38 -bed > nanog.sites.hg38.bed # motif file is from homer: Homer/data/knownTFs/motifs
scanMotifGenomeWide.pl klf5.motif hg38 -bed > klf5.sites.hg38.bed # motif file is from homer: Homer/data/knownTFs/motifs

# obtain genome size
# samtools faidx Homo_sapiens.GRCh38.97.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.97.dna.primary_assembly.fa.fai > Homo_sapiens.GRCh38.97.dna.primary_assembly.size

# convert motif.sites.bed to bw
bedtools merge -i nanog.sites.hg38.bed -c 5 -o mean > nanog.sites.hg38.merge.bed
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' nanog.sites.hg38.merge.bed > nanog.sites.hg38.bedgraph
sort -k1,1 -k2,2n nanog.sites.hg38.bedgraph > nanog.sites.hg38_srt.bedgraph
bedGraphToBigWig nanog.sites.hg38_srt.bedgraph Homo_sapiens.GRCh38.97.dna.primary_assembly.size nanog.sites.hg38.bw 
# the bw file is used to plot heatmap




### plot heatmap over peak regions
computeMatrix reference-point \
-S bigwigfile1.bw bigwigfile2.bw \
-R Peakregion1.bed Peakregion2.bed \
-p 20 \
-a 5000 -b 5000 --skipZeros \
-o heatmap_peak.mat.gz

plotHeatmap -m heatmap_peak.mat.gz \
-out heatmap_peak.pdf \
--boxAroundHeatmaps no \
--colorMap Blues \
--missingDataColor "#b5ddff"











