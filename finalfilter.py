import sys
from Bio import SeqIO
from pythonmods import inctyping, rmlstprofile,rmlsttypingalleles

plasmidfinderdatabases=sys.argv[4:]
rmlstprofilepath=sys.argv[1]
outdir=sys.argv[2]
sequenceorigin=sys.argv[3]

#also need to include whether using inhouse sequences or ncbi (is there accessions_filtered file?) ? should be _deduplicated?

#get replicon types
enterobacaccessions=[]
gram_posaccessions=[]
enterobacaccessionsdict={}
gram_posaccessionsdict={}
for database in plasmidfinderdatabases:
    with open('%s/plasmidfinder/BLASTtablebesthits_%s.tsv'%(outdir,database)) as f:
        for line in f:
            data=line.strip().split('\t')
            accession=data[0]
            allele=data[1]
            if database=="enterobacteriaceae":
                enterobacaccessions.append(accession)
                if accession in enterobacaccessionsdict:
                    enterobacaccessionsdict[accession].append(allele)
                else:
                    enterobacaccessionsdict[accession]=[]
                    enterobacaccessionsdict[accession].append(allele)
            else:
                gram_posaccessions.append(accession)
                if accession in gram_posaccessionsdict:
                    gram_posaccessionsdict[accession].append(allele)
                else:
                    gram_posaccessionsdict[accession]=[]
                    gram_posaccessionsdict[accession].append(allele)
                
enterobacaccessions=list(set(enterobacaccessions))
gram_posaccessions=list(set(gram_posaccessions))
repaccessions=list(set(enterobacaccessions+gram_posaccessions))

#do replicon typing
reptypedict={}
for accession in repaccessions:
    reptypedict[accession]=[]
    if accession in enterobacaccessionsdict:
        types,families,probes,length=inctyping(accession, enterobacaccessionsdict,db="enterobacteriaceae")
    else:
        types,families,probes,length=['-','-','-','-']
    reptypedict[accession].extend([types,probes,length])
    if accession in gram_posaccessionsdict:
        types,families,probes,length=inctyping(accession, enterobacaccessionsdict,db="gram_positive")
    else:
        types,families,probes,length=['-','-','-','-']
    reptypedict[accession].extend([types,probes,length])  #not including families
    

#get rmlst types
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

rmlstaccessions=list(set(rmlstaccessions))


#do rmlst typing
rmlstprofileoutput=rmlstprofile('%s/profiles.txt'%rmlstprofilepath)
rmlsttypedict={}
for accession in rmlstaccessions:
    alleles=rmlstaccessionsdict[accession]
    rmlsttype=rmlsttypingalleles(alleles,rmlstprofileoutput)
    rmlsttypedict[accession]=rmlsttype #a tuple of species, rST top match, num_matches, num_mismatches, num_missing loci, rmlstalleles



###

#output final files (dependent on sequenceorigin)

if sequenceorigin=='ncbi':
    ###non plasmid accessions; output tsv with rmlst and replicon typing info
    f2=open('%s/nonplasmidaccessions.tsv'%outdir,'w')
    f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\tSpecies\tribosomalST\tNum_Matches\tNum_Mismatches\tNum_Missing_Loci\trMLST_alleles\tEnterobacteriaceae_Type\tProbe_Hits\tNum_Probe_hits\tGram_positive_Type\tProbe_Hits\tNum_Probe_hits\n')
    with open('%s/accessions_filtered.tsv'%outdir) as f:
        for indx, line in enumerate(f):
            if indx==0:
                continue
            data=line.strip().split('\t')
            accession=data[0]
            if accession in rmlstaccessions: #non plasmid accession
                #get rmlst type
                #species,rst,toprst,num_matches,num_mismatches,num_missingloci,rmlstalleles=rmlsttypedict[accession]
                rmlsttype=rmlsttypedict[accession]
                #get replicon type
                if accession in reptypedict:
                    reptype=reptypedict[accession]
                else:
                    reptype=('-','-','-','-','-','-')
                #write to file
                f2.write('%s\t%s\t%s\n'%('\t'.join(data),'\t'.join(rmlsttype),'\t'.join(reptype)))            
    f2.close()
    ###plasmid accessions; output fasta and tsv with replicon typing info
    f2=open('%s/plasmids.fa'%outdir,'w')
    plasmidaccessions=[]
    with open('%s/accessions_filtered_deduplicated.fa'%outdir) as f:
        for indx, seq_record in enumerate(SeqIO.parse(f,"fasta")):
            fastaheader=str(seq_record.id)
            accession=fastaheader.split(' ')[0].lstrip('>')
            if accession in rmlstaccessions:
                continue
            plasmidaccessions.append(accession)
            SeqIO.write(seq_record, f2, "fasta")
    f2.close()
    
    f2=open('%s/plasmids.tsv'%outdir,'w')
    f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\tEnterobacteriaceae_Type\tProbe_Hits\tNum_Probe_hits\tGram_positive_Type\tProbe_Hits\tNum_Probe_hits\n')
    with open('%s/accessions_filtered.tsv'%outdir) as f:
        for indx,line in enumerate(f):
            if indx==0:
                continue
            data=line.strip().split('\t')
            accession=data[0]
            if accession in plasmidaccessions:
                if accession in reptypedict:
                    reptype=reptypedict[accession]
                else:
                    reptype=('-','-','-','-','-','-')
                f2.write('%s\t%s\n'%('\t'.join(data),'\t'.join(reptype)))
    f2.close()


else:
    ###output tsv with rmlst and replicon typing info
    #get seqlengths
    seqlengthdict={}
    accessions=[]
    with open('%s/seqlengths.tsv'%outdir) as f:
        for line in f:
            data=line.strip().split('\t')
            accession=data[0]
            length=data[1]
            seqlengthdict[accession]=length
            accessions.append(accession)
    #write tsv
    f2=open('%s/typing.tsv'%outdir,'w')
    f2.write('Accession\tLength\t\tSpecies\tribosomalST\tNum_Matches\tNum_Mismatches\tNum_Missing_Loci\trMLST_alleles\tEnterobacteriaceae_Type\tProbe_Hits\tNum_Probe_hits\tGram_positive_Type\tProbe_Hits\tNum_Probe_hits\n')
    for accession in accessions:
        #get replicon types
        if accession in reptypedict:
            reptype=reptypedict[accession]
        else:
            reptype=('-','-','-','-','-','-')
        #get rmlst type
        if accession in rmlsttypedict:
            rmlsttype=rmlsttypedict[accession]
        else:
            rmlsttype=('-','-','-','-','-','-')
        f2.write('%s\t%s\t%s\t%s\n'%(accession,seqlengthdict[accession],'\t'.join(rmlsttype),'\t'.join(reptype)))
    f2.close()

#OLD CODE
# #get rmlstonly / rmlstrep accessions

# #print(rmlstaccessions.intersection(repaccessions), 'rmlst accessions that have rep type')
# #print(rmlstaccessions.difference(repaccessions), 'rmlst accessions without rep type')

# rmlstrepaccessions=sorted(list(rmlstaccessions.intersection(repaccessions)))
# rmlstonlyaccessions=sorted(list(rmlstaccessions.difference(repaccessions)))

# #write accessions to file along with info from accessions_filtered.tsv
# f2=open('%s/rmlstrepaccessions.tsv'%outdir,'w')
# f3=open('%s/rmlstonlyaccessions.tsv'%outdir,'w')
# f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\trMLST_alleles\n')
# f3.write('Accession\tTopology\tLength\tTitle\tCompleteness\trMLST_alleles\n')

# with open('%s/accessions_filtered.tsv'%outdir) as f:
#     for indx, line in enumerate(f):
#         if indx==0:
#             continue
#         data=line.strip().split('\t')
#         accession=data[0]
#         if accession in rmlstrepaccessions:
#             alleles='|'.join(rmlstaccessionsdict[accession])
#             data.append(alleles)
#             f2.write('%s\n'%'\t'.join(data))
#         if accession in rmlstonlyaccessions:
#             alleles='|'.join(rmlstaccessionsdict[accession])
#             data.append(alleles)
#             f3.write('%s\n'%'\t'.join(data))
            
        
# f2.close()
# f3.close()



# #filter all rmlst accessions - these are likely to be either misannotated chromosomal sequences or chromids (see diCenzo 2017; Harrison 2010 papers on chromids)

# f2=open('%s/plasmids.fa'%outdir,'w')
# #f3=open('%s/plasmids.txt'%outdir,'w')
# accessions=[]
# with open('%s/accessions_filtered_deduplicated.fa'%outdir) as f:
#     for indx, seq_record in enumerate(SeqIO.parse(f,"fasta")):
#         fastaheader=str(seq_record.id)
#         accession=fastaheader.split(' ')[0].lstrip('>')
#         if accession in rmlstaccessions:
#             continue
#         accessions.append(accession)
#         SeqIO.write(seq_record, f2, "fasta")
# f2.close()


# f2=open('%s/plasmids.tsv'%outdir,'w')
# f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\n')
# with open('%s/accessions_filtered.tsv'%outdir) as f:
#     for indx,line in enumerate(f):
#         if indx==0:
#             continue
#         data=line.strip().split('\t')
#         accession=data[0]
#         if accession in accessions:
#             f2.write(line)
# f2.close()

# #OLD
# #don't filter at this stage - first manually check rmlstrep accessions then filter rmlstonly + any dodgy rmlstrep

