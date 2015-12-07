seqcluster cluster -a align/star_sim.bam -m sim.ma -o star_res -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
seqcluster cluster -a align/bowtie_sim.bam -m sim.ma -o bowtie_res -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
seqcluster cluster -a align/bowtie2_sim.bam -m sim.ma -o bowtie2_res -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
seqcluster cluster -a align/hisat2_sim.bam -m sim.ma -o hisat2_res -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
