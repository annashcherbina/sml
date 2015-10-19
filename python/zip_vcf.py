#KEEPS ONLY THE VARIANT CALLS FROM INPUT VCF FILE  
import sys 
data=open(sys.argv[1],'r').read().split('\n') 
while '' in data: 
    data.remove('') 
outf=open(sys.argv[2],'w')
#WRITE THE HEADER! 
outf.write('uniqID\tGene\tcis-eQTL_PVAL\tRsquare\tBETA\tRsquare>Population\tBETA>Population\tTSS\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNA12877-200-37-ASM\tNA12878-200-37-ASM\tNA12879-200-37-ASM\tNA12880-200-37-ASM\tNA12881-200-37-ASM\tNA12882-200-37-ASM\tNA12883-200-37-ASM\tNA12884-200-37-ASM\tNA12885-L2-200-37-ASM\tNA12886-L2-200-37-ASM\tNA12887-L2-200-37-ASM\tNA12888-200-37-ASM\tNA12889-L2-200-37-ASM\tNA12890-200-37-ASM\tNA12891-200-37-ASM\tNA12892-L2-200-37-ASM\tNA12893-200-37-ASM\tSubjectWithVar\n')  
first_person_index=14 
for line in data: 
    if line.startswith('#'): 
        continue
    line=line.split('\t') 
    if len(line)<9: 
        continue
    outstring='\t'.join(line[0:first_person_index-1])
    #JUST STORE THE GENOTYPE! 
    has1=0 
    for alleles in line[first_person_index::]: 
        alleles=alleles.split(':')[0] 
        if alleles.__contains__('1'): 
            has1+=1 
        outstring=outstring+'\t'+alleles 
    outstring=outstring+'\t'+str(has1)+'\n'
    outf.write(outstring) 
    
