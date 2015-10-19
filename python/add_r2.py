data=open('/users/annashch/vcfchunks/NOMAJOR.motifs.CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
if '' in data: 
    data.remove('') 
ref=open('/users/annashch/effects/CEU.family1463.gene.cis-eQTL.txt','r').read().split('\n') 
ref_dict=dict() 
for line in ref[1::]: 
    tokens=line.split('\t') 
    ref_dict[tokens[
