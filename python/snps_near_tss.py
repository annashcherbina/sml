from Parameters import * 
import sys 
def main(): 
    dist=100000
    binsize=50000
    effect_data=open(effect,'r').read().split('\n') 
    print "loaded effect data" 
    #vcf_data=open(vcf,'r').read().split('\n') 
    vcf_data=open(sys.argv[1],'r').read().split('\n') 
    outf=open(sys.argv[2],'w') 
    print "loaded VCF data" 
    #BUILD DICTIONARY OF TSS POSITIONS 
    tss_pos=dict()
    started=False 
    for line in effect_data: 
        if started==False: 
            if line.__contains__('XLOC'): 
                started=True 
        else: 
            line=line.split('\t')
            if len(line)<5: 
                continue 
            #print str(line) 
            chrom=line[2] 
            pos=int(line[3]) 
            cur_bin=pos/binsize
            if chrom not in tss_pos: 
                tss_pos[chrom]=dict() 
            if cur_bin not in tss_pos[chrom]: 
                tss_pos[chrom][cur_bin]=dict() 
            tss_pos[chrom][cur_bin][pos]=[line[1],line[5],line[7],line[9]]
    #print "tss_pos:"+str(tss_pos) 
    print "Built dictionary of transcription start sites" 
    #curl=0 
    #total=str(len(vcf_data)) 
    for line in vcf_data: 
        #curl+=1
        #if curl%100==0: 
        #    print str(curl)+":"+total 
        if line.startswith('#'): 
            outf.write(line+'\n') 
        else: 
            tokens=line.split('\t') 
            if len(line)<8: 
                continue 
            chrom=tokens[0]
            pos=int(tokens[1]) 
            cur_bin=pos/binsize 
            all_bins=[cur_bin-1,cur_bin,cur_bin+1] 
            if chrom not in tss_pos: 
                continue 
            best_dist=dist+1 
            closest_gene=None 
            best_meta=None 
            for b in all_bins: 
                if b not in tss_pos[chrom]: 
                    continue 
                candidates=tss_pos[chrom][b] 
                for c in candidates: 
                    if abs(c-pos)< best_dist: 
                        best_dist=abs(c-pos) 
                        closest_gene=c 
                        best_meta=candidates[closest_gene] 
            if best_dist < dist: 
                outf.write('\t'.join(best_meta)+'\t'+str(best_dist)+'\t'+line+'\n') 

if __name__=="__main__": 
    main() 
