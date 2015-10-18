import sys 
from helpers import * 
#start_pos=int(sys.argv[1]) 
#end_pos=int(sys.argv[2]) 
#GROUP BY MAF & PROPORTION OF FAMILY EFFECT  > POPULATION  
data=open('/users/annashch/vcfchunks/NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
#data=open('NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n')

while '' in data: 
    data.remove('')
print "Read in source data" 
#BUILD SUMMARY DICTIONARIES OF WHICH SNPs SATISFY WHICH REQUIREMENTS

terms=['dnase','tf','splice50','splice30','tss','phylop1','phylop2','cadd10','motif']
special=build_categories(terms)
for key in special:
    for binnum in range(7):
        special[key][binnum]=[0,0]
print "set up summary dict"

beta_pop_cutoff=95 
cl=0 
for line in data[1::]: 
    cl+=1 
    if cl%10000==0: 
        print str(cl) 
    tokens=line.split('\t') 
    motif=int(tokens[0])
    cadd=tokens[1]
    if cadd=="NA":
        cadd=float("-inf")
    else:
        cadd=float(tokens[1])
    phylop=tokens[2]
    if phylop=="NA":
        phylop=float("-inf")
    else: 
        phylop=float(tokens[2]) 
    maf=tokens[3] 
    maf_bin=get_maf_bin(maf) 
    dnase=int(tokens[4]) 
    tf=int(tokens[5]) 
    splice_dist=tokens[6] 
    if splice_dist=="inf": 
        splice_dist=float("inf") 
    else: 
        splice_dist=float(splice_dist) 
    tss_dist=int(tokens[11])
    
    beta_pop=float(tokens[10]) 
    if beta_pop > beta_pop_cutoff: 
        beta_add=1
    else: 
        beta_add=0 
    matches=[] 
    #MOTIF! 
    if motif==1:
        matches.append('motif')
    if cadd>10:
        matches.append('cadd10')
    if phylop > 1:
        matches.append('phylop1')
    if phylop > 2:
        matches.append('phylop2')
    if dnase==1:
        matches.append('dnase')
    if tf==1:
        matches.append('tf')
    if splice_dist <50:
        matches.append('splice50')
    if splice_dist <30:
        matches.append('splice30')
    if tss_dist <5000:
        matches.append('tss')
    #print str(matches) 
    for key in special:
        hasall=True
        for group in key:
            if group not in matches:
                hasall=False
                break
        if hasall:
            special[key][maf_bin][0]+=beta_add
            special[key][maf_bin][1]+=1
#print str(special) 
print "writing output file:"
outf=open('/users/annashch/summary/summary.tsv','w')
#outf=open('summary.tsv','w') 
outf.write('Factor1\tFactor2\tFactor3\tFactor4\tFactor5\tFactor6\tFactor7\tFactor8\tFactor9\tMAF_0_.05\tMAF_.05_.10\tMAF_.10_.15\tMAF_.15_.20\tMAF.20_\tMAF_NA\n')

#SUMMARIZE THE RESULTS!
for key in special:
    #key=list(key)
    for i in range(9):
        if i<len(key):
            outf.write(key[i]+'\t')
        else:
            outf.write('\t') 
    for binval in range(7):
        if special[key][binval][1]==0:
            fract="NA"
            outf.write('\tNA')
        else: 
            fract=special[key][binval][0]/float(special[key][binval][1])
            outf.write('\t'+str(round(fract,3)))
    outf.write('\n')
    
