snps=open('/users/annashch/out.INFO','r').read().split('\n') 
if '' in snps: 
    snps.remove('') 
print "Read in snps"
sites=open('/users/annashch/exons.gtf','r').read().split('\n') 
if '' in sites: 
    sites.remove('') 
splice_pos=dict() 
binsize=50
thresh=50 
for line in sites: 
    line=line.split('\t')
    print str(line) 
    chrom=line[0]  
    if chrom not in splice_pos: 
        splice_pos[chrom]=dict() 
    startpos=int(line[3])
    endpos=int(line[4]) 
    start_bin=startpos/binsize 
    end_bin=endpos/binsize 
    if start_bin not in splice_pos[chrom]: 
        splice_pos[chrom][start_bin]=dict() 
    if end_bin not in splice_pos[chrom]: 
        splice_pos[chrom][end_bin]=dict() 
    splice_pos[chrom][start_bin][startpos]=1 
    splice_pos[chrom][end_bin][endpos]=1 
print "built dictionary of splice positions!" 
outf=open('/users/annashch/out.INFO.splice','w') 
for line in sites: 
    tokens=line.split('\t') 
    chrom=tokens[0] 
    pos=int(tokens[1]) 
    posbin=pos/binsize 
    hit=False 
    for b in [posbin-1, posbin,posbin+1]: 
        if hit==True: 
            break 
        if b not in splice_pos[chrom]: 
            continue 
        candidates=splice_pos[chrom][b] 
        for c in candidates: 
            if abs(c-pos)<thresh: 
                hit=True 
                break 
    if hit==True: 
        outf.write(line+'\n') 
