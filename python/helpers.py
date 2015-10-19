def get_maf_bin(maf): 
    if maf=="NA": 
        return 4 
    else:
        maf=float(maf) 
        maf=min(maf,1-maf)#sometimes the minor allele is reference, and sometimes alternate
        if maf  < 0.01: 
            return 0
        elif maf < 0.05: 
            return 1 
        elif maf < 0.2: 
            return 2
        else: 
            return 3 

def build_categories(terms):
    term_dict=dict()
    lt=len(terms)
    for t1 in range(lt):
        cur_key=[terms[i] for i in [t1]]
        term_dict[tuple(cur_key)]=dict() 
        for t2 in range(t1+1,lt):
            cur_key=[terms[i] for i in [t1,t2]]
            term_dict[tuple(cur_key)]=dict() 
            for t3 in range(t2+1,lt) :
                cur_key=[terms[i] for i in [t1,t2,t3]]
                term_dict[tuple(cur_key)]=dict() 
                for t4 in range(t3+1,lt):
                    cur_key=[terms[i] for i in [t1,t2,t3,t4]]
                    term_dict[tuple(cur_key)]=dict() 
                    for t5 in range(t4+1,lt):
                        cur_key=[terms[i] for i in [t1,t2,t3,t4,t5]]
                        term_dict[tuple(cur_key)]=dict() 
                        for t6 in range(t5+1,lt):
                            cur_key=[terms[i] for i in [t1,t2,t3,t4,t5,t6]]
                            term_dict[tuple(cur_key)]=dict() 
                            for t7 in range(t6+1,lt):
                                cur_key=[terms[i] for i in [t1,t2,t3,t4,t5,t6,t7]]
                                term_dict[tuple(cur_key)]=dict() 
                                for t8 in range(t7+1,lt):
                                    cur_key=[terms[i] for i in [t1,t2,t3,t4,t5,t6,t7,t8]]
                                    term_dict[tuple(cur_key)]=dict() 
                                    for t9 in range(t8+1,lt):
                                        cur_key=[terms[i] for i in [t1,t2,t3,t4,t5,t6,t7,t8,t9]]
                                        term_dict[tuple(cur_key)]=dict()
    print str(len(term_dict))
    return term_dict
 
        
