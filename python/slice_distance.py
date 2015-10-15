data=open("/users/annashch/vcfchunks/full",'r').read().split('\n') 
while '' in data: 
    data.remove('') 
gtf=open('/users/annashch/genes.gtf','r').read().split('\n') 
while '' in gtf: 
    gtf.remove('') 
exon_boundaries=dict() 
binsize=100 
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
    
