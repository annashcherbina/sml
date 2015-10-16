#APPENDS MAF INFORMATION TO ANNOTATION FILE 
maf_data=open("/users/annashch/out.INFO",'r').read().split('\n') 
while '' in maf_data: 
    maf_data.remove('') 
maf_dict=dict() 
for line in maf_data[1::]: 
    line=line.split('\t') 
    chrom=line[0] 
    pos=line[1] 
    maf=line[4] 
    if maf=="?": 
        maf=line[5] # use global maf as a fallback if MAF for European is unknown 
        
