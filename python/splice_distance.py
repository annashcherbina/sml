data=open("/users/annashch/vcfchunks/full",'r').read().split('\n') 
while '' in data: 
    data.remove('') 
gtf=open('/users/annashch/genes.gtf','r').read().split('\n') 
while '' in gtf: 
    gtf.remove('') 
exon_boundaries=dict() 
binsize=500  
for line in gtf: 
    line=line.split("\t") 
    chrom=line[0] 
    exon_start=int(line[3]) 
    exon_end=int(line[4]) 
    if chrom not in exon_boundaries: 
        exon_boundaries[chrom]=dict() 
    #hash exon start and end sites to bins of 100 bases 
    start_bin=exon_start/binsize 
    end_bin=exon_end/binsize 
    if start_bin not in exon_boundaries[chrom]: 
        exon_boundaries[chrom][start_bin]=dict() 
    if end_bin not in exon_boundaries[chrom]: 
        exon_boundaries[chrom][end_bin]=dict() 
    exon_boundaries[chrom][start_bin][exon_start]=1
    exon_boundaries[chrom][end_bin][exon_end]=1 
print "build dictionary of exon ends" 
outf=open("/users/annashch/vcfchunks/splice.full","w")
header=data[0] 
outf.write('SpliceDistance\t'+header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t')
    #print str(tokens) 
    chrom=tokens[5] 
    try:
        pos=int(tokens[6])
    except: 
        continue 
    #print "chrom:"+str(chrom) 
    #print "pos:"+str(pos) 
    pos_bin=pos/binsize 
    bins=[pos_bin-1,pos_bin,pos_bin+1] 
    min_dist=float("inf") 
    for b in bins: 
        if b not in exon_boundaries[chrom]: 
            continue 
        candidates=exon_boundaries[chrom][b] 
        #print str(candidates) 
        for c in candidates: 
            dist=abs(c-pos)
            #print str(dist) 
            if dist < min_dist: 
                min_dist=dist                 
    outf.write(str(min_dist)+'\t'+str(line)+'\n') 

        

