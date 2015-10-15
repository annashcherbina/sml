data=open('/users/annashch/bed/ceu.bed','r').read().split('\n') 
while '' in data: 
    data.remove('') 
print "read in our BED file"
bed_dict=dict() 
for line in data: 
    line=line.split('\t') 
    chrom=line[0] 
    pos=line[1] 
    if chrom not in bed_dict: 
        bed_dict[chrom]=dict() 
    bed_dict[chrom][pos]=1
print "built bed dictionary"

ref=open('/users/annashch/bed/phyloP100way.bedgraph','r')
outf=open('intersection.ceu','w') 
c=0 
for line in ref: 
    c+=1
    if c%1000000==0: 
        print str(c) 
    tokens=line.split('\t') 
    chrom=tokens[0] 
    pos=tokens[1] 
    if chrom not in bed_dict: 
        continue 
    if pos in bed_dict[chrom]: 
        outf.write(line+'\n') 

