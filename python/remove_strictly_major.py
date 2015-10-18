#REMOVES SNPS WHERE EVERYONE IN THE FAMILY HAS ONLY MAJOR ALLELES! 
data=open('/users/annashch/vcfchunks/motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
outf=open('/users/annashch/vcfchunks/NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','w')
while '' in data: 
    data.remove('') 

header=data[0] 
outf.write(header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t') 
    if len(tokens)<38: 
        print str(tokens) 
        continue 
    maf=tokens[3] 
    subjects_with_alt=int(tokens[37]) 
    if (maf=="NA") and (subjects_with_alt > 0): 
        outf.write(line+'\n')
    elif (maf < 0.5) and (subjects_with_alt >0): 
        outf.write(line+'\n') 
    else: 
        # THE ALTERNATIVE ALLELE IS THE MAJOR ALLELE! 
        if (subjects_with_alt < 17): #NOT EVERYONE HAS THE MAJOR ALLELE! 
            outf.write(line+'\n') 
