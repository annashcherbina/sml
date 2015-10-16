#ADDS Phylop and CADD conservation scores 
cadd=open('/users/annashch/CEUFAMsnp012.012.pos.CADD.txt','r').read().split('\n') 
while '' in cadd: 
    cadd.remove('') 
cadd_dict=dict() 
for line in cadd: 
    tokens=line.split('\t') 
    chrom=tokens[0] 
    pos=tokens[1] 
    ref=tokens[2] 
    alt=tokens[3] 
    cons=tokens[5] 
    if chrom not in cadd_dict: 
        cadd_dict[chrom]=dict() 
    if pos not in cadd_dict[chrom]: 
        cadd_dict[chrom][pos]=dict() 
        cadd_dict[chrom][pos]['mean']=0 
    if ref not in cadd_dict[chrom][pos]: 
        cadd_dict[chrom][pos][ref]=dict() 
    cadd_dict[chrom][pos][ref][alt]=cons 
    cadd_dict[chrom][pos]['mean']+=float(cons) 
    #print str(cadd_dict[chrom][pos])

        
print "build cadd dictionary" 

phylop=open('/users/annashch/phylop.scores','r').read().split('\n') 
while '' in phylop: 
    phylop.remove('') 
phylop_dict=dict() 
for line in phylop: 
    tokens=line.split('\t') 
    chrom=tokens[0].replace('chr','') 
    pos=tokens[1] 
    score=tokens[3] 
    if chrom not in phylop_dict: 
        phylop_dict[chrom]=dict() 
    phylop_dict[chrom][pos]=score 
    
print "build phylop dictionary" 

data=open('/users/annashch/vcfchunks/MAF.DNAse.TF.splice.full','r').read().split('\n') 
while '' in data: 
    data.remove('') 
outf=open('/users/annashch/vcfchunks/CONSERVED.MAF.DNAse.TF.splice.full','w')
header=data[0] 
outf.write('CADD\tPhylop\t'+header+'\n')
for line in data[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[9]
    pos=tokens[10]
    ref=tokens[12] 
    alt=tokens[13] 

    #print "chrom:"+str(chrom) 
    #print "pos:"+str(pos) 
    cadd="NA" 
    phylop="NA"
    if (chrom in cadd_dict) and (pos in cadd_dict[chrom]): 
        if (ref in cadd_dict[chrom][pos]) and (alt in cadd_dict[chrom][pos][ref]):
            cadd=cadd_dict[chrom][pos][ref][alt]
        else: 
            cadd=str(cadd_dict[chrom][pos]['mean']/3)
    if (chrom in phylop_dict) and (pos in phylop_dict[chrom]): 
        phylop=phylop_dict[chrom][pos] 
    outf.write(cadd+'\t'+phylop+'\t'+line+'\n') 



