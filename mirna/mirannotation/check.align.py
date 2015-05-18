#!/usr/bin/python

#load simulated reads
sim=open('sim.20.hsa.fa','r')
data={}
dataseq={}
lendata={}
name=""
for line in sim:
    line=line.strip()
    if line.startswith(">"):
        name=line.replace(">","")
        slot=name.split("_")[0].split("-")
        data[name]="-".join(slot[0:3])
    else:
        lendata[name]=len(line)
        #print name
        dataseq[line]=name
sim.close()

#load miraligner results
check={}
mir=open('miraligner/sim.20.hsa.mirna')
for line in mir:
    cols=line.split("\t")
    slot=cols[2].split("-")
    add=line.find("add:null")
    mut=line.find("mut:null")
    if (line.find("hsa")>=0):
        if (line.find("miRNA")>=0) and (not check.has_key(cols[1])):
                check[cols[1]]=1
                print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(cols[1],cols[2],cols[1].split("_")[1],add,mut,lendata[cols[1]])
        elif (line.find("precursor")>=0) and (not check.has_key(cols[1])):
                check[cols[1]]=1
                slot=cols[2].split("-")
                print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(cols[1],cols[2],cols[1].split("_")[0],add,mut,lendata[cols[1]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmiraligner" % (k,data[k],add,mut,lendata[k])

#load bwotie2 results
check={}
mir=open('bowtie2/sim.20.hsa.sam')
for line in mir:
    cols=line.split("\t")
    if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
            check[cols[0]]=1
            slot=cols[2].split("-")
            add=cols[0].find("add:null")
            mut=cols[0].find("mut:null")
            print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tbowtie2" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tbowtie2" % (k,data[k],add,mut,lendata[k])

#load bwotie results
check={}
mir=open('bowtie/sim.20.hsa.sam')
for line in mir:
    cols=line.split("\t")
    if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
            check[cols[0]]=1
            slot=cols[2].split("-")
            add=cols[0].find("add:null")
            mut=cols[0].find("mut:null")
            print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tbowtie" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tbowtie" % (k,data[k],add,mut,lendata[k])


#load novoaligner results
check={}
mir=open('novo/sim.20.novo.sam')
for line in mir:
    cols=line.split("\t")
    if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
            slot=cols[2].split("-")
            add=cols[0].find("add:null")
            mut=cols[0].find("mut:null")
            check[cols[0]]=1
            print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tnovocraft" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tnovocraft" % (k,data[k],add,mut,lendata[k])


#load srnabench results
check={}
mir=open('srnabench/reads_orig.fa')
srna={}
name=""
for line in mir:
    line=line.strip()
    if (line.find(">")>=0):
            name=line.replace(">","")
    else:
            srna[name]=line
mir.close()

mir=open('srnabench/hairpin.parsed')
for line in mir:
    cols=line.split("\t")
    seq=srna[cols[0]]
    if (cols[2].find("hsa")>=0) and (not check.has_key(dataseq[seq])):
            slot=cols[2].split("-")
            add=dataseq[seq].find("add:null")
            mut=dataseq[seq].find("mut:null")
            check[dataseq[seq]]=1
            print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tsrnabench" %(dataseq[seq],cols[2],dataseq[seq].split("_")[0],add,mut,lendata[dataseq[seq]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tsrnabench" % (k,data[k],add,mut,lendata[k])


#load GEM results
check={}
mir=open('gem/sim.20.hsa.sam')
for line in mir:
    if line.find("@")<0:
            cols=line.split("\t")
            if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
                    slot=cols[2].split("-")
                    add=cols[0].find("add:null")
                    mut=cols[0].find("mut:null")
                    check[cols[0]]=1
                    print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tGEM" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tGEM" % (k,data[k],add,mut,lendata[k])

#load microrazer results
check={}
mir=open('microrazer/sim.20.hsa.res')
for line in mir:
    cols=line.split("\t")
    if (cols[4].find("hsa")>=0) and (not check.has_key(cols[0])):
        slot=cols[2].split("-")
        add=cols[0].find("add:null")
        mut=cols[0].find("mut:null")
        check[cols[0]]=1
        print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmicrorazer" %(cols[0],cols[4],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
            add=k.find("add:null")
            mut=k.find("mut:null")
            print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmicrorazer" % (k,data[k],add,mut,lendata[k])


#load STAR results
check={}
mir=open('star/Aligned.out.sam')
for line in mir:
    cols=line.split("\t")
    if not line.startswith("@"):
        if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
            check[cols[0]]=1
            slot=cols[2].split("-")
            add=cols[0].find("add:null")
            mut=cols[0].find("mut:null")
            print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tstar" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tstar" % (k,data[k],add,mut,lendata[k])

# load miRExpress results
check={}
mir=open('mirexpress/out/read_align_premiRNA.txt')
for line in mir:
    if line.startswith("hsa"):
        ref = line.strip()
        mir.next()
        mir.next()
        mir.next()
        continue
    cols = line.split()
    if (line.find("hsa")<0) and (line.find("10") >= 0) and (not check.has_key(cols[0])):
        name = dataseq[cols[0]]
        slot=name.split("-")
        add=name.find("add:null")
        mut=name.find("mut:null")
        check[cols[0]] = 1
        print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmirexpress" %(name,ref,name.split("_")[0],add,mut,lendata[name])
mir.close()

# if not aligned, print them here
for k in dataseq.keys():
    if not check.has_key(k):
        name = dataseq[k]
        add=name.find("add:null")
        mut=name.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmirexpress" % (name,data[name],add,mut,lendata[name])


