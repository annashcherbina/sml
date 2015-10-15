import sys 
#SOURCE FILE OF VARIANT INFORMATION 
data=open("/users/annashch/vcfchunks/FULL","r").read().split('\n') 
header=data[0] 

#APPENDS MAF INFORMATION TO ANNOTATION FILE 
maf_data=open("/users/annashch/out.INFO",'r').read().split('\n') 
print "read in maf file" 
while '' in maf_data: 
    maf_data.remove('') 
maf_dict=dict() 
c=0 
total=str(len(maf_data)) 
for line in maf_data[1::]: 
    c+=1 
    if c%1000==0: 
        print str(c)+":"+total 
    line=line.split('\t') 
    chrom=line[0] 
    pos=line[1] 
    maf=line[4] 
    if maf=="?": 
        maf=line[5] # use global maf as a fallback if MAF for European is unknown 
    if chrom not in maf_dict: 
        maf_dict[chrom]=dict() 
    maf_dict[chrom][pos]=maf 
print "built dictionary of minor allele frequencies" 

outf=open('/users/annashch/vcfchunks/annotated.txt','w') 
outf.write('MAF\tCADD\tPhylop\tMOTIF\tSPLICE_REGION\tTF\tDNASEI\t'+header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t') 
    varchrom=tokens[5] 
    varpos=tokens[6] 

    #GET THE MINOR ALLELE FREQUENCY 
    if varchrom not in maf_dict: 
        print "CHROM NOT FOUND!:"+str(line) 
        continue 
    if varpos not in maf_dict: 
        print "POS NOT FOUND!:"+str(line) 
        continue 
    mafval=maf_dict[varchrom][varpos]


    outf.write(mafval+'\t'+line+'\n') 

