import sys
from Bio import Entrez, SeqIO

plasmidfinderdatabases=sys.argv[2:]
outdir=sys.argv[1]

#get rmlstonly / rmlstrep accessions
repaccessions=[]
for database in plasmidfinderdatabases:
    with open('%s/plasmidfinder/BLASTtablebesthits_%s.tsv'%(outdir,database)) as f:
        for line in f:
            data=line.strip().split('\t')
            repaccessions.append(data[0])

repaccessions=set(repaccessions)

rmlstaccessions=[]
rmlstaccessionsdict={}
with open('%s/rmlst/BLASTtablebesthits.tsv'%outdir) as f:
    for line in f:
        data=line.strip().split('\t')
        accession=data[0]
        allele=data[1]
        rmlstaccessions.append(accession)
        if accession in rmlstaccessionsdict:
            rmlstaccessionsdict[accession].append(allele)
        else:
            rmlstaccessionsdict[accession]=[]
            rmlstaccessionsdict[accession].append(allele)

rmlstaccessions=set(rmlstaccessions)

#print(rmlstaccessions.intersection(repaccessions), 'rmlst accessions that have rep type')
#print(rmlstaccessions.difference(repaccessions), 'rmlst accessions without rep type')

rmlstrepaccessions=sorted(list(rmlstaccessions.intersection(repaccessions)))
rmlstonlyaccessions=sorted(list(rmlstaccessions.difference(repaccessions)))

#write accessions to file along with info from accessions_filtered.tsv
f2=open('%s/rmlstrepaccessions.tsv'%outdir,'w')
f3=open('%s/rmlstonlyaccessions.tsv'%outdir,'w')

with open('%s/accessions_filtered.tsv'%outdir) as f:
    for line in f:
        data=line.strip().split('\t')
        accession=data[0]
        if accession in rmlstrepaccessions:
            alleles='|'.join(rmlstaccessionsdict[accession])
            data.append(alleles)
            f2.write('%s\n'%'\t'.join(data))
        if accession in rmlstonlyaccessions:
            alleles='|'.join(rmlstaccessionsdict[accession])
            data.append(alleles)
            f3.write('%s\n'%'\t'.join(data))
            
        
f2.close()
f3.close()



#filter all rmlst accessions - these are likely to be either misannotated chromosomal sequences or chromids (see diCenzo 2017; Harrison 2010 papers on chromids)

f2=open('%s/accessions_final.fa'%outdir,'w')
f3=open('%s/accessions_final.txt'%outdir,'w')
accessions=[]
with open('%s/accessions_filtered.fa'%outdir) as f:
    for indx, seq_record in enumerate(SeqIO.parse(f,"fasta")):
        fastaheader=str(seq_record.id)
        accession=fastaheader.split(' ')[0].lstrip('>')
        if accession in rmlstaccessions:
            continue
        accessions.append(accession)
        SeqIO.write(seq_record, f2, "fasta")
        f3.write('%s\n'%accession)
f2.close()
f3.close()

f2=open('%s/accessions_final.tsv'%outdir,'w')
with open('%s/accessions_filtered.tsv'%outdir) as f:
    for line in f:
        data=line.strip().split('\t')
        accession=data[0]
        if accession in accessions:
            f2.write(line)
f2.close()

#OLD
#don't filter at this stage - first manually check rmlstrep accessions then filter rmlstonly + any dodgy rmlstrep

