#ADDS MOTIF INFORMATION FOR SNP ANNOTATION 
motifs=open('/users/annashch/CEUfamily.MOTIFS.flattened.txt','r').read().split('\n') 
if '' in motifs: 
    motifs.remove('')
motif_dict=dict() 
for line in motifs: 
    tokens=line.split('\t') 
    chrom=tokens[0]
    if chrom not in motif_dict: 
        motif_dict[chrom]=dict() 
    pos=tokens[1]
    motif_dict[chrom][pos]=1 


data=open('/users/annashch/vcfchunks/CONSERVED.MAF.DNAse.TF.splice.full','r').read().split('\n') 
while '' in data: 
    data.remove('') 
header=data[0] 
outf=open('/users/annashch/vcfchunks/motifs.CONSERVED.MAF.DNAse.TF.splice.full','w')
outf.write('MOTIF\t'+header+'\n') 
for line in data[1::]: 
    tokens=line.split('\t') 
    chrom=tokens[11]
    pos=tokens[12] 
    motif=0 
    if chrom in motif_dict: 
        if pos in motif_dict[chrom]: 
            motif=1
    outf.write(str(motif)+'\t'+line+'\n') 
