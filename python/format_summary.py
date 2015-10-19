data=open('/users/annashch/summary/summary.tsv','r').read().split('\n')
data.remove('') 
outf=open('/users/annashch/summary/formatted.summary.tsv','w')
categories=['dnase','tf','splice50','splice30','tss','phylop1','phylop2','cadd10','motif']
categories.sort() 
header='NumAnnotations\t'+'\t'.join(categories) +'\tMAF_0.01\tMAF_0_.05\tMAF_.05_.10\tMAF_.10_.15\tMAF_.15_.20\tMAF.20_\tMAF_NA\n'
outf.write(header) 
for line in data[1::]: 
    tokens=line.split('\t') 
    present_cat=tokens[0:9] 
    numentries=9-present_cat.count('') 
    outf.write(str(numentries)) 
    for c in categories: 
        if c in present_cat: 
            outf.write('\t1') 
        else: 
            outf.write('\t') 
    outf.write('\t'+'\t'.join(tokens[10::])+'\n')

