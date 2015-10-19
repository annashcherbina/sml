#REMOVES SNPS WHERE EVERYONE IN THE FAMILY HAS ONLY MAJOR ALLELES! 
data=open('/users/annashch/vcfchunks/motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
outf=open('/users/annashch/vcfchunks/NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','w')
while '' in data: 
    data.remove('') 

header=data[0] 
header_parts=header.split('\t') 
maf_index=header_parts.index("MAF") 
sum_index=header_parts.index("SubjectWithVar") 
outf.write(header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t') 
    if len(tokens)<(sum_index+1): 
        print str(tokens) 
        continue 
    maf=tokens[maf_index] 
    subjects_with_alt=int(tokens[sum_index]) 
    print "maf:"+str(maf) 
    print "subjects_with_alt:"+str(subjects_with_alt) 
    if maf=="NA": 
        if subjects_with_alt > 0: 
            outf.write(line+'\n') 
        else: 
            print "skipping!" 
    else: 
        maf=float(maf) 
        if (float(maf) < 0.5) and (subjects_with_alt >0): 
            outf.write(line+'\n') 
        elif (float(maf) > 0.5) and (subjects_with_alt < 17): 
            outf.write(line+'\n') 
        else: 
            print "skipping!" 
