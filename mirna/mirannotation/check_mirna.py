#!/usr/bin/python

#load simulated reads
sim=open('sim.21.hsa.fa','r')
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
check, save = {}, {}
mir=open('miraligner/sim.21.hsa.mirna')
for line in mir:
    cols=line.split("\t")
    slot=cols[2].split("-")
    add=line.find("add:null")
    mut=line.find("mut:null")
    score=len(cols[5]) - (len(cols[5]) - 1) + len(cols[6]) - (len(cols[6]) - 1)
    name = cols[1]
    if name not in check:
        check[name] = 100
    if (line.find("hsa")>=0):
        if line.find("miRNA") >=0 and score < check[name]:
                check[name]=score
                save[name] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(name,cols[2],name.split("_")[1],add,mut,lendata[name])
        elif (line.find("precursor")>=0) and score < check[name]:
                check[name]=score
                slot=cols[2].split("-")
                save[name] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(name,cols[2],name.split("_")[0],add,mut,lendata[name])
mir.close()

for k in save:
    print save[k]

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmiraligner" % (k,data[k],add,mut,lendata[k])

#load isomir-sea results
check, save = {}, {}
mir=open('isomir-sea/sim.out.txt')
for line in mir:
    cols=line.split("\t")
    ref = cols[5][1:].split(" ")[0]
    add=line.find("add:null")
    mut=line.find("mut:null")
    score=0
    if (line.find("hsa")>=0):
        seq=cols[1].replace("U", "T")
        name=dataseq[seq]
        add=name.find("add:null")
        mut=name.find("mut:null")
        check[seq]=1
        print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tisomir-sea" %(name,ref,name.split("_")[1],add,mut,lendata[name])
mir.close()

# if not aligned, print them here
for k in dataseq.keys():
    if not check.has_key(k):
        name = dataseq[k]
        add=name.find("add:null")
        mut=name.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tisomir-sea" % (name,data[name],add,mut,lendata[name])


#load miraligner-python results
check, save = {}, {}
mir=open('miraligner-python/res/sim.21.hsa.mirna')
h = mir.next()
for line in mir:
    cols=line.split("\t")
    slot=cols[3].split("-")
    add=line.find("add:null")
    mut=line.find("mut:null")
    score=len(cols[6]) - (len(cols[6]) - 1) + len(cols[7]) - (len(cols[7]) - 1)
    name = cols[1]
    # if name not in check:
    #    check[name] = 100
    if (line.find("hsa")>=0) and name not in check:
        check[name] = score
        save[name] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner-python" %(name,cols[3],name.split("_")[1],add,mut,lendata[name])
mir.close()

for k in save:
    print save[k]

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmiraligner-python" % (k,data[k],add,mut,lendata[k])

#load chimira_blast results
check, save = {}, {}
mir=open('chimira_blast/sim.hsa.blast_out.txt')
for line in mir:
    cols=line.split("\t")
    slot=cols[1].split("-")
    add=line.find("add:null")
    mut=line.find("mut:null")
    score=float(cols[-2])
    name = cols[0]
    if name not in check:
        check[name] = 100
    if (line.find("hsa")>=0) and score < check[name]:
        check[name]=score
        save[name] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\tblast" %(name,cols[1],name.split("_")[1],add,mut,lendata[name])
mir.close()

for k in save:
    print save[k]

#if not aligned, print them here
for k in data.keys():
    if not check.has_key(k):
        add=k.find("add:null")
        mut=k.find("mut:null")
        print "%s\t%s\tno\tNA\t%s\t%s\t%s\tblast" % (k,data[k],add,mut,lendata[k])


def _sam(fn, tool):
    check, save = {}, {}
    with open(fn) as in_handle:
        for line in in_handle:
            cols=line.split("\t")
            if line.startswith("@"):
                continue
            if cols[2].find("hsa") == -1:
                continue
            if cols[0] not in check:
                check[cols[0]] = 100
            score = [len(flag) for flag in cols if flag.startswith('MD')]
            if score:
                score = score[0]
            else:
                score = 0
            if (cols[2].find("hsa")>-1) and score < check[cols[0]]:
                    check[cols[0]]=score
                    slot=cols[2].split("-")
                    add=cols[0].find("add:null")
                    mut=cols[0].find("mut:null")
                    save[cols[0]] = "%s\t%s\tyes\t%s\t%s\t%s\t%s\t%s" %(cols[0],cols[2],cols[0].split("_")[0],add,mut,lendata[cols[0]], tool)

    for k in save:
        print save[k]
    #if not aligned, print them here
    for k in data.keys():
        if k not in check:
            add=k.find("add:null")
            mut=k.find("mut:null")
            print "%s\t%s\tno\tNA\t%s\t%s\t%s\t%s" % (k,data[k],add,mut,lendata[k], tool)


#load razer3 results
_sam('razer/sim.21.hsa.sam', 'razer3')

#load bwotie2 results
_sam('bowtie2/sim.21.hsa.sam', 'bowtie2')

#load bwotie results
_sam('bowtie/sim.21.hsa.sam', 'bowtie')

#load novoaligner results
_sam('novo/sim.21.hsa.sam', 'novo')

#load GEM results
_sam('gem/sim.21.hsa.sam', 'GEM')

#load STAR results
_sam('star/Aligned.out.sam', 'STAR')

#load STAR results
_sam('bwa/sim.21.hsa.sam', 'bwa')

#load srnabench results
check={}
mir=open('srnabench/srnabench/reads_orig.fa')
srna={}
name=""
for line in mir:
    line=line.strip()
    if (line.find(">")>=0):
            name=line.replace(">","")
    else:
            srna[name]=line
mir.close()

mir=open('srnabench/srnabench/hairpin.parsed')
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



#load microrazer results
check={}
mir=open('microrazer/sim.21.hsa.res')
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



# load miRExpress results
check={}
mir=open('mirexpress/out/hsa.sim.21')
for line in mir:
    if line.startswith("hsa"):
        ref = line.strip()
        mir.next()
        mir.next()
        # mir.next()
        continue
    cols = line.split()
    if (line.find("hsa")<0) and (line.find("1") >= 0) and (not check.has_key(cols[0])):
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


