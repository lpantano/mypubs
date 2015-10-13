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
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O hairpin.fa.gz
wget -q ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -O miRNA.str.gz

gunzip -f hairpin.fa.gz miRNA.str.gz

awk '{if ($0~/>/){print name"\t"seq;name=$1;seq=""}else{seq=seq$0}}' hairpin.fa | grep "^>hsa" | sed "s/U/T/g" | sed 's/\t/\n/' > hairpin.hsa.fa

if [ ! -e sim.20.hsa.fa ] ; then

    wget https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/miRNA.simulator.py -O miRNA.simulator.py
    python miRNA.simulator.py -f hairpin.hsa.fa -m miRNA.str -n 10 -s hsa > sim.20.hsa.fa

fi

mkdir -p miraligner
cd miraligner
 miraligner  -sub 1 -trim 3 -add  3 -s hsa -i ../sim.20.hsa.fa -db ../ -o sim.20.hsa -pre
cd ..

mkdir -p bowtie2
cd bowtie2
 bowtie2-build ../hairpin.hsa.fa mirbowtie2
 bowtie2 -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 -f -x mirbowtie2 -U ../sim.20.hsa.fa -S sim.20.hsa.sam
cd ..

mkdir -p bowtie
cd bowtie
 bowtie-build ../hairpin.hsa.fa mirbowtie
 bowtie  -a --best --strata  -n 1 -l 16 -f mirbowtie  ../sim.20.hsa.fa -S sim.20.hsa.sam
cd ..

mkdir -p novo
cd novo
 novoindex -k 14 -s 2 -m mirnovo ../hairpin.hsa.fa  
 novoalign -d mirnovo -F FA -f ../sim.20.hsa.fa -a -l 15 -t 30 -o SAM >sim.20.hsa.sam
cd ..

mkdir -p gem
cd gem
 gem-indexer -i ../hairpin.hsa.fa  -o mirgem --complement 'yes' 
 gem-mapper -e 0.1 -m 0.1 -I mirgem.gem -i ../sim.20.hsa.fa -o sim.20.hsa 
 gem-2-sam   -I mirgem.gem -i sim.20.hsa.map -o sim.20.hsa.sam
cd ..

mkdir -p srnabench
cd srnabench
java -jar $RNAB/sRNAbench.jar  dbPath=$RNAB_DB microRNA=hsa input=../sim.20.hsa.fa output=srnabench isoMiR=true
cd ..

mkdir -p microrazer
cd microrazer
micro_razers -o sim.20.hsa.res -sL 12 ../hairpin.hsa.fa ../sim.20.hsa.fa
cd ..


mkdir -p razer
cd razer
razers3 -dr 0 -i 80 -rr 90 -f -o sim.20.hsa.sam ../hairpin.hsa.fa ../sim.20.hsa.fa
cd ..

mkdir -p mirexpress
cd mirexpress
awk '{if ($0!~/^>/){print 1"\t"$0}}' ../sim.20.hsa.fa > sim.20.hsa.trim
mkdir -p out
alignmentSIMD -r ../data_v2/hsa_precursor.txt -i sim.20.hsa.trim -o out
analysis -r ../data_v2/hsa_precursor.txt -m ../data_v2/hsa_miRNA.txt -d out -o /hsa.sim.20 -t hsa.sim.20.profile
cd ..

mkdir -p star
cd star
mkdir -p genome
STAR --runMode genomeGenerate --genomeDir genome --genomeFastaFiles ../hairpin.hsa.fa
STAR --genomeDir genome --readFilesIn ../sim.20.hsa.fa --alignIntronMax 1  --outFilterMultimapNmax 1000 --outSAMattributes NH HI NM --outSAMtype SAM SortedByCoordinate
cd ..

mkdir -p miraligner-python
cd miraligner-python
seqcluster seqbuster --mirna ../miRNA.str --hairpin ../hairpin.hsa.fa --sps hsa -o res ../sim.20.hsa.fa
cd ..

mkdir -p chimira_blast
cd chimira_blast
ln -fs ../hairpin.hsa.fa .
formatdb -p F -i hairpin.hsa.fa
blastall -p blastn -e 0.001 -d hairpin.hsa.fa -i ../sim.20.hsa.fa -m 8 -o chimira_input.gz.blast_out.txt -v 1 -b 3 -F F -W 7 
cd ..

python check.align.py > stats

echo "now time to run stats.rmd"


