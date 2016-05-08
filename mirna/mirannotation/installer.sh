# All tools not listed heres are installed by brew at homebrew-cbl and homebrew-science

mkdir -p  tools ; cd tools
wget http://barnaserver.com/gemtools/releases/GEMTools-static-i3-1.7.1.tar.gz
tar xzvf GEMTools-static-i3-1.7.1.tar.gz

wget http://bioinfo5.ugr.es/sRNAbench/sRNAbench.jar
wget http://bioinfo5.ugr.es/sRNAbench/sRNAbenchDB.tgz
tar xzvf sRNAbenchDB.tgz

wget https://github.com/seqan/seqan/releases/download/seqan-v1.4.2/seqan-apps-1.4.2-Linux-x86_64.tar.bz2
tar xjfv seqan-apps-1.4.2-Linux-x86_64.tar.bz2

wget http://mirexpress.mbc.nctu.edu.tw/Download/miRExpress.2.1.4.tar.gz
tar xzvf miRExpress.2.1.4.tar.gz
cd miRExpress
./configure
make
cd ..

wget -nc  http://mirexpress.mbc.nctu.edu.tw/Download/mirbase21_precursor.tar.gz
tar xzvf mirbase21_precursor.tar.gz

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/netblast-2.2.26-x64-linux.tar.gz
tar xzvf blast-2.2.26-x64-linux.tar.gz
