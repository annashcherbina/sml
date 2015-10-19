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
    for binnum in range(5):
        special[key][binnum]=[0,0]
print "set up summary dict"

beta_pop_cutoff=95 
cl=0 
for line in data[int(sys.argv[1]):int(sys.argv[2])]: 
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
    if beta_pop >= beta_pop_cutoff: 
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
categories=['dnase','tf','splice50','splice30','tss','phylop1','phylop2','cadd10','motif']
categories.sort() 
header='NumAnnotations\t'+'\t'.join(categories) +'\tMAF<0.01_Total\tMAF<0.01_Sig\tMAF<0.01_Prop\tMAF_0.01_0.05_Total\tMAF_0.01_0.05_Sig\tMAF_0.01_0.05_Prop\t'
header=header+'MAF_0.05_0.20_Total\tMAF_0.05_0.20_Sig\tMAF_0.05_0.20_Prop\t'
header=header+'MAF_>0.20_Total\tMAF_>0.20_Sig\tMAF_>0.20_Prop\t'
header=header+'MAF_NA_Total\tMAF_NA_Sig\tMAF_NA_Prop\t'
outf.write(header+'\n') 

#SUMMARIZE THE RESULTS!
for key in special:
    print str(key) 
    numentries=len(key) 
    outf.write(str(numentries)) 
    for c in categories: 
        if c in key: 
            outf.write('\t1') 
        else: 
            outf.write('\t0') 
    for binval in range(5):
        totalval=special[key][binval][1] 
        sigval=special[key][binval][0] 
        if totalval==0: 
            prop=0
        else: 
            prop=float(sigval)/totalval 
        outf.write('\t'+str(totalval)+'\t'+str(sigval)+'\t'+str(prop))
    outf.write('\n')
    
