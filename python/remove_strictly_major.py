#REMOVES SNPS WHERE EVERYONE IN THE FAMILY HAS ONLY MAJOR ALLELES! 
data=open('/users/annashch/vcfchunks/motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
outf=open('/users/annashch/vcfchunks/NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','w')
header=data[0] 
outf.write(header+'\n') 

