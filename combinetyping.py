import sys
from mymod import openpickle

enterobacdict=openpickle('%s/typing/pickle.incdict_enterobacteriaceae'%sys.argv[1])
gramposdict=openpickle('%s/typing/pickle.incdict_gram_positive'%sys.argv[1])

f2=open('%s/typing/inctypes.tsv'%sys.argv[1],'w')
with open('./%s/accessions_final.txt'%(sys.argv[1])) as f:
    for line in f:
        sample=line.strip()
        enterobac=enterobacdict[sample]
        grampos=gramposdict[sample]
        newlist=[]
        for indx, (a,b) in enumerate(zip(enterobac,grampos)):
            if indx==3:
                newlist.append(a)
                newlist.append(b)
                if a=='-' and b=='-':
                    newlist.append('-')
                elif a=='-':
                    newlist.append(b)
                elif b=='-':
                    newlist.append(a)
                else:
                    newlist.append(a+b)
                continue
            if a=='-' and b=='-':
                newlist.append('-')
                continue
            if a=='-':
                newlist.append(b)
                continue
            if b=='-':
                newlist.append(a)
                continue
            newlist.append(a+','+b)
        f2.write('%s\t%s\n'%(sample, '\t'.join(newlist)))
f2.close()
