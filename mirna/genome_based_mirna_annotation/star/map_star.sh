STAR --genomeDir /groups/bcbio/bcbio/genomes/Hsapiens/hg19/star --readFilesIn $1 --outFilterMultimapNmax 50 --outSAMattributes NH HI NM --alignIntronMax 1
samtools view -Sbh Aligned.out.sam >| sim.20.hsa.bam
bedtools intersect -bed -wo -s -f 0.80 -abam sim.20.hsa.bam -b hsa.gff3 >| sim.20.hsa.bam.anno
