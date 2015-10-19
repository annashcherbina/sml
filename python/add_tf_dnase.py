tfbinding=open("/users/annashch/tf_binding_sites","r").read().split('\n') 
dnase_peaks=open("/users/annashch/ensembl_open_chromatin","r").read().split('\n') 

while ''in tfbinding: 
    tfbinding.remove('') 
while '' in dnase_peaks: 
    dnase_peaks.remove('') 

tf_dict=dict()
dnase_dict=dict() 
 
binsize=100
for line in tfbinding: 
    line=line.split('\t') 
    chrom=line[0].replace('chr','') 
    startpos=int(line[3]) 
    endpos=int(line[4]) 
    if chrom not in tf_dict: 
        tf_dict[chrom]=dict() 
    binval=startpos/binsize 
    if binval not in tf_dict[chrom]: 
        tf_dict[chrom][binval]=dict() 
    tf_dict[chrom][binval][startpos]=endpos 
print "built transcription factor dictionary" 

for line in dnase_peaks: 
    line=line.split('\t') 
    chrom=line[0].replace('chr','') 
    startpos=int(line[3])  
    endpos=int(line[4]) 
    if chrom not in dnase_dict: 
        dnase_dict[chrom]=dict() 
    binval=startpos/binsize 
    if binval not in dnase_dict[chrom]: 
        dnase_dict[chrom][binval]=dict() 
    dnase_dict[chrom][binval][startpos]=endpos 
print "built dnase_peak dictionary" 

data=open("/users/annashch/vcfchunks/splice.full",'r').read().split('\n') 
while '' in data: 
    data.remove('') 
print 'read in splice data' 
outf=open('/users/annashch/vcfchunks/DNAse.TF.splice.full','w') 
outf.write('DNAse\tTF\t'+data[0]+'\n') 
header=data[0].split('\t') 
chrom_index=header.index("CHROM") 
pos_index=header.index("POS") 
countval=0 
for line in data[1::]: 
    countval+=1 
    if countval%1000==0: 
        print str(countval) 
    tokens=line.split('\t') 
    chrom=tokens[chrom_index]
    base=int(tokens[pos_index]) 
    base_bin=base/binsize 
    bins=[base_bin-1,base_bin,base_bin+1] 
    tf_bind=0
    dnase_peak=0
    if chrom in tf_dict: 
        for b in bins: 
            if (tf_bind==0) and (b in tf_dict[chrom]): 
                candidates=tf_dict[chrom][b] 
                for c in candidates: 
                    if (base >=c) and (base<=candidates[c]): 
                        #it's in a TF! 
                        tf_bind=1
                        break 
    if chrom in dnase_dict: 
        for b in bins: 
            if (dnase_peak==0) and (b in dnase_dict[chrom]): 
                candidates=dnase_dict[chrom][b]
                for c in candidates: 
                    if (base>=c) and (base<=candidates[c]): 
                        #it's in a DNAse peak! 
                        dnase_peak=1
                        break 
    outf.write(str(dnase_peak)+'\t'+str(tf_bind)+'\t'+line+'\n')

        
