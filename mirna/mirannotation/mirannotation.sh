#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

. update_path.bash

RNAB_DB="../tools/sRNAbenchDB"
RNAB="../tools/"

echo "download miRBase files"
if [ ! e hairpin.hsa.fa ] ; then
    wget -qN ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O hairpin.fa.gz
    wget -qN ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz -O mature.fa.gz
    wget -qN ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -O miRNA.str.gz

    gunzip -f hairpin.fa.gz miRNA.str.gz

    awk '{if ($0~/>/){print name"\t"seq;name=$1;seq=""}else{seq=seq$0}}' hairpin.fa | grep "^>hsa" | sed "s/U/T/g" | sed 's/\t/\n/' > hairpin.hsa.fa
fi

if [ ! -e sim.21.hsa.fa ] ; then

    wget https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/miRNA.simulator.py -O miRNA.simulator.py
    python miRNA.simulator.py -f hairpin.hsa.fa -m miRNA.str -n 10 -s hsa -p sim.21.hsa.fa

fi

if [ ! -e  "miraligner/sim.21.hsa.mirna" ]; then
    SECONDS=0
    mkdir -p miraligner
    cd miraligner
     miraligner  -sub 1 -trim 3 -add  3 -s hsa -i ../sim.21.hsa.fa -db ../ -o sim.21.hsa -pre
    cd ..
    duration=$SECONDS
    echo $(($duration / 60)) minutes elapsed
fi

if [ ! -e  "bowtie2/sim.21.hsa.sam" ]; then
    SECONDS=0
    mkdir -p bowtie2
    cd bowtie2
     bowtie2-build ../hairpin.hsa.fa mirbowtie2
     bowtie2 -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 -f -x mirbowtie2 -U ../sim.21.hsa.fa -S sim.21.hsa.sam
    cd ..
    duration=$SECONDS
    echo $(($duration / 60)) minutes elapsed
fi

if [ ! -e "bowtie/sim.21.hsa.sam" ]; then
    mkdir -p bowtie
    cd bowtie
     bowtie-build ../hairpin.hsa.fa mirbowtie
     bowtie  -a --best --strata  -n 1 -l 16 -f mirbowtie  ../sim.21.hsa.fa -S sim.21.hsa.sam
    cd ..
fi

if [ ! -e  "novo/sim.21.hsa.sam" ]; then
    mkdir -p novo
    cd novo
     novoindex -k 14 -s 2 -m mirnovo ../hairpin.hsa.fa  
     novoalign -d mirnovo -F FA -f ../sim.21.hsa.fa -a -l 15 -t 30 -o SAM > sim.21.hsa.sam
    cd ..
fi

if [ ! -e  "gem/sim.21.hsa.sam" ]; then
    mkdir -p gem
    cd gem
     gem-indexer -i ../hairpin.hsa.fa  -o mirgem --complement 'yes' 
     gem-mapper -e 0.1 -m 0.1 -I mirgem.gem -i ../sim.21.hsa.fa -o sim.21.hsa 
     gem-2-sam   -I mirgem.gem -i sim.21.hsa.map -o sim.21.hsa.sam
    cd ..
fi

if [ ! -e  "srnabench/hairpin.parsed" ]; then
    mkdir -p srnabench
    cd srnabench
    java -jar $RNAB/sRNAbench.jar  dbPath=$RNAB_DB microRNA=hsa input=../sim.21.hsa.fa output=srnabench isoMiR=true
    cp srnabench/hairpin.parsed .
    cd ..
fi

if [ ! -e  "microrazer/sim.21.hsa.res" ]; then
    mkdir -p microrazer
    cd microrazer
    micro_razers -o sim.21.hsa.res -sL 12 ../hairpin.hsa.fa ../sim.21.hsa.fa
    cd ..
fi

if [ ! -e  "razer/sim.21.hsa.sam" ]; then
    mkdir -p razer
    cd razer
    razers3 -dr 0 -i 80 -rr 90 -f -o sim.21.hsa.sam ../hairpin.hsa.fa ../sim.21.hsa.fa
    cd ..
fi

if [ ! -e  "mirexpress/out/hsa.sim.21" ]; then
    mkdir -p mirexpress
    cd mirexpress
    # if [ ! -e data2 ]; then
    # wget -N http://mirexpress.mbc.nctu.edu.tw/Download/mirbase21_precursor.tar.gz %% tar -xvcf mirbase21_precursor.tar.gz
    awk '{if ($0!~/^>/){print 1"\t"$0}}' ../sim.21.hsa.fa > sim.21.hsa.trim
    mkdir -p out
    alignmentSIMD -r ../tools/data_v2/hsa_precursor.txt -i sim.21.hsa.trim -o out
    analysis -r ../tools/data_v2/hsa_precursor.txt -m ../tools/data_v2/hsa_miRNA.txt -d out -o /hsa.sim.21 -t hsa.sim.21.profile
    cd ..
fi

if [ ! -e  "star/Aligned.out.sam" ]; then
    mkdir -p star
    cd star
    mkdir -p genome
    STAR --runMode genomeGenerate --genomeDir genome --genomeFastaFiles ../hairpin.hsa.fa
    STAR --genomeDir genome --readFilesIn ../sim.21.hsa.fa --alignIntronMax 1  --outFilterMultimapNmax 1000 --outSAMattributes NH HI NM --outSAMtype SAM --sysShell /bin/bash --readFilesCommand cat
    cd ..
fi

if [ ! -e  "bwa/sim.21.hsa.sam" ]; then
    mkdir -p bwa
    cd bwa
    ln -sf ../hairpin.hsa.fa .
    bwa index hairpin.hsa.fa -a is
    bwa aln  -l 15 -M 1 -R 100 -n 3  hairpin.hsa.fa ../sim.21.hsa.fq  > sim.21.hsa.sai
    bwa  samse hairpin.hsa.fa sim.21.hsa.sai  ../sim.21.hsa.fq > sim.21.hsa.sam
    cd ..
fi

if [ ! -e  "miraligner-python/res/sim.21.hsa.mirna" ]; then
    mkdir -p miraligner-python
    cd miraligner-python
    seqcluster seqbuster --mirna ../miRNA.str --hairpin ../hairpin.hsa.fa --sps hsa -o res ../sim.21.hsa.fa
    cd ..
fi

if [ ! -e  "chimira_blast/sim.hsa.blast_out.txt" ]; then
    mkdir -p chimira_blast
    cd chimira_blast
    ln -fs ../hairpin.hsa.fa .
    formatdb -p F -i hairpin.hsa.fa
    blastall -p blastn -e 0.001 -d hairpin.hsa.fa -i ../sim.21.hsa.fa -m 8 -o sim.hsa.blast_out.txt -v 1 -b 3 -F F -W 7 
    cd ..
fi

if [ ! -e "isomir-sea/sim.out.txt" ]; then
    mkdir -p isomir-sea/data/out
    cd isomir-sea
    awk '{if ($0!~/^>/){print $0"\t"10}}' ../sim.21.hsa.fa > data/input.txt 
    cp ../mature.fa data/mature_21.txt
    ../tools/isomiR-SEA_OS/Debian_linux5_x86_64/isomiR-SEA_1_6 -s hsa -l 16 -b 4 -i data -p data/out -ss 6 -h 11 -m mature_21 -t input
    cat data/out/out_result_mature_21_tag_ambigue_selected.txt data/out/out_result_mature_21_tag_unique.txt > sim.out.txt
    cd ..
fi

python check_mirna.py > stats

echo "now time to run stats.rmd"

Rscript -e 'library(rmarkdown);render("stats.rmd")'
Rscript -e 'library(knitr);knit("stats.rmd")'
Rscript -e 'library(rmarkdown);render("stats_isomirs.rmd")'
Rscript -e 'library(knitr);knit("stats_isomirs.rmd")'
