import sys 
data=open(sys.argv[1],'r').read().split('\n') 
while '' in data: 
    data.remove('')  
outf=open(sys.argv[2],'w') 
for line in data: 
    tokens=line.split('\t') 
    try:
        pos=int(tokens[1]) 
    except: 
        print str(line) 
        continue 
    pos2=pos+1 
    outf.write("chr"+line+'\t'+str(pos2)+'\n') 
