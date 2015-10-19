#APPENDS MAF INFORMATION TO ANNOTATION FILE 
maf_data=open("/users/annashch/out.INFO",'r').read().split('\n') 
while '' in maf_data: 
    maf_data.remove('') 
maf_dict=dict() 
for line in maf_data[1::]: 
    line=line.split('\t') 
    chrom=line[0] 
    if chrom not in maf_dict: 
        maf_dict[chrom]=dict() 
    pos=line[1] 
    maf=line[4] 
    if maf=="?": 
        maf=line[5] # use global maf as a fallback if MAF for European is unknown 
    maf_dict[chrom][pos]=maf 
    #print "chrom:"+str(chrom) 
    #print "pos:"+str(pos) 
    #print "maf:"+str(maf) 
print "build MAF dictionary" 
outf=open('/users/annashch/vcfchunks/MAF.DNAse.TF.splice.full','w') 
data=open('/users/annashch/vcfchunks/DNAse.TF.splice.full','r').read().split('\n') 
while '' in data: 
    data.remove('') 
header=data[0]
header_parts=header.split('\t') 
chrom_index=header_parts.index("CHROM") 
pos_index=header_parts.index("POS") 
outf.write("MAF\t"+header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[chrom_index] 
    pos=tokens[pos_index] 
    if chrom in maf_dict: 
        if pos in maf_dict[chrom]: 
            outf.write(maf_dict[chrom][pos]+'\t'+line+'\n') 
        else: 
            outf.write('NA\t'+line+'\n') 
    else: 
        outf.write('NA\t'+line+'\n') 

    
