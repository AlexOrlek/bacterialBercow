import sys
import numpy as np
from Bio import SeqIO
from pythonmods import inctyping, rmlstprofile,rmlsttypingalleles


rmlstprofilepath=sys.argv[1]
outdir=sys.argv[2]
sequenceorigin=sys.argv[3]
if sequenceorigin=='ncbi':
    typing='both' #both replicon and rMLST typing are performed in order to curate NCBI sequences
    plasmidfinderdatabases=sys.argv[4:6]
else: #in-house
    typing=sys.argv[4]
    plasmidfinderdatabases=sys.argv[5:7]
    contigcompletenessfile=sys.argv[7]
    contigsamplesfile=sys.argv[8]

#filtering depends on whether using ncbi or inhouse sequences, and if the latter, whether both rmlst and replicon typing have been conducted

###replicon typing
if typing=='both' or typing=='replicon':
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
            types,families,probes,length=inctyping(accession, gram_posaccessionsdict,db="gram_positive")
        else:
            types,families,probes,length=['-','-','-','-']
        reptypedict[accession].extend([types,probes,length])  #not including families


###rmlst typing
if typing=='both' or typing=='rmlst':
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
        rmlsttypedict[accession]=rmlsttype #a tuple of top species, top rST, num_matches, num_mismatches, num_missing loci, num_multiallelic loci, rmlstalleles



###

#output final files (dependent on sequenceorigin)

if sequenceorigin=='ncbi':
    ###non plasmid accessions; output tsv with rmlst and replicon typing info + contig_description ('chromosome' or 'chromosomal')
    plasmidaccessions=[]
    f2=open('%s/nonplasmids.tsv'%outdir,'w')
    f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\tSpecies\tribosomalST\tNum_Matches\tNum_Mismatches\tNum_Missing_Loci\tNum_MultiallelicLoci\trMLST_alleles\tEnterobacteriaceae_Type\tProbe_Hits\tNum_Probe_hits\tGram_positive_Type\tProbe_Hits\tNum_Probe_hits\tContig_Classification\n')
    with open('%s/accessions_filtered_deduplicated.tsv'%outdir) as f:
        for indx, line in enumerate(f):
            if indx==0:
                continue
            data=line.strip().split('\t')
            accession=data[0]
            length=data[2]
            if accession in rmlstaccessions: #non plasmid accession
                #get rmlst type, and typing stats
                rmlsttype=rmlsttypedict[accession]
                topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
                #get replicon type
                if accession in reptypedict:
                    reptype=reptypedict[accession]
                else:
                    reptype=('-','-','-','-','-','-')
                ###classify contig
                if int(num_missingloci)<3 and int(length)>500000: #smallest bacterial genome is ~580Kb
                    contigclass='chromosome'
                else:
                    contigclass='chromosomal'
                #write to file
                f2.write('%s\t%s\t%s\t%s\n'%('\t'.join(data),'\t'.join(rmlsttype),'\t'.join(reptype),contigclass))
            else:
                plasmidaccessions.append(accession)
    f2.close()
    ###plasmid accessions; output fasta and tsv with replicon typing info
    #output fasta
    f2=open('%s/plasmids.fa'%outdir,'w')
    with open('%s/accessions_filtered_deduplicated.fa'%outdir) as f:
        for indx, seq_record in enumerate(SeqIO.parse(f,"fasta")):
            fastaheader=str(seq_record.id)
            accession=fastaheader.split(' ')[0].lstrip('>')
            if accession in rmlstaccessions:
                continue
            SeqIO.write(seq_record, f2, "fasta")
    f2.close()
    #output tsv
    f2=open('%s/plasmids.tsv'%outdir,'w')
    f2.write('Accession\tTopology\tLength\tTitle\tCompleteness\tEnterobacteriaceae_Type\tProbe_Hits\tNum_Probe_hits\tGram_positive_Type\tProbe_Hits\tNum_Probe_hits\n')
    with open('%s/accessions_filtered_deduplicated.tsv'%outdir) as f:
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
    ###output tsv with rmlst and/or replicon typing info + contig classification if rmlst is conducted
    #get seqlengths
    seqlengthdict={}
    contigs=[]
    with open('%s/seqlengths.tsv'%outdir) as f:
        for line in f:
            data=line.strip().split('\t')
            contig=data[0]
            length=int(data[1])
            seqlengthdict[contig]=length
            contigs.append(contig)
    #write tsv
    f2=open('%s/typing.tsv'%outdir,'w')
    #first write column header
    header=['Contig','Length']
    if typing=='both' or typing=='rmlst':
        header.extend(['Species','ribosomalST','Num_Matches','Num_Mismatches','Num_Missing_Loci','Num_MultiallelicLoci','rMLST_alleles'])
    if typing=='both' or typing=='replicon':
        header.extend(['Enterobacteriaceae_Type','Probe_Hits','Num_Probe_hits','Gram_positive_Type','Probe_Hits','Num_Probe_hits'])
    if typing=='both' or typing=='rmlst':
        header.extend(['Contig_Classification'])
    header='\t'.join(header)
    f2.write('%s\n'%header)
    #if only replicon typing is specified, write rep typing info to file, but don't attempt to classify contig
    if typing=='replicon':
        for contig in sorted(contigs):
            #get replicon types
            if contig in reptypedict:
                reptype=reptypedict[contig]
            else:
                reptype=('-','-','-','-','-','-')
            #write to file
            f2.write('%s\t%s\t%s\n'%(contig,seqlengthdict[contig],'\t'.join(reptype)))
        f2.close()
    else:
        ###get reptype if both rmlst and reptyping are specified
        if typing=='both':
            for contig in contigs:
                #get replicon types
                if contig in reptypedict:
                    reptype=reptypedict[contig]
                else:
                    reptype=('-','-','-','-','-','-')
        ###extract contig completeness info if available
        contigcompletenessdict={}
        if contigcompletenessfile!='None':
            with open(contigcompletenessfile) as f:
                for indx,line in enumerate(f):
                    data=line.strip().split('\t')
                    contig=data[0].strip()
                    completeness=data[1].strip()
                    if contig not in contigs and indx==0:
                        continue #assume header line
                    if contig not in contigs:
                        sys.exit('Error: contig name: %s does not match any fasta header names'%contig)
                    if completeness.lower() in ['complete','circular','complete_linear']:
                        completeness='complete'
                    elif completeness.lower() in ['incomplete','linear','unknown']:
                        completeness='incomplete'
                    else:
                        sys.exit('Error: un-recognised contig completeness annotation: %s'%completeness)
                    contigcompletenessdict[contig]=completeness
        ###iterate through samples if sample groupings are available (longest contig within sample likely to be chromosome)
        if contigsamplesfile!='None':
            samplecontigdict={}
            with open(contigsamplesfile) as f:
                for indxa,line in enumerate(f):
                    data=line.strip().split('\t')
                    contig=data[0].strip()
                    sample=data[1].strip()
                    if contig not in contigs and indxa==0:
                        continue #assume header line
                    if contig not in contigs:
                        sys.exit('Error: contig name: %s does not match any fasta header names'%contig)
                    if sample in samplecontigdict:
                        samplecontigdict[sample].append(contig)
                    else:
                        samplecontigdict[sample]=[]
                        samplecontigdict[sample].append(contig)
            #iterate through samples
            for sample in sorted(samplecontigdict.keys()):
                contigs=sorted(samplecontigdict[sample])
                lengths=[]
                for contig in contigs:
                    lengths.append(seqlengthdict[contig])
                indices=list(np.argsort(lengths))
                indices.reverse()
                contigs=[contigs[i] for i in indices] #contigs are now sorted in length order, starting with longest
                lengths=[lengths[i] for i in indices]
                #iterate through contigs, sorted by length
                for indxb,(contig,length) in enumerate(zip(contigs,lengths)):
                    if contig in rmlsttypedict:
                        rmlsttype=rmlsttypedict[accession]
                    else:
                        rmlsttype=('-','-','-','-','-','-','-')
                    topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
                    if contig in contigcompletenessdict:
                        contigcompleteness=contigcompletenessdict[contig]
                    else:
                        contigcompleteness='unknown'
                    #check for chromosome contig: longest contig, rmlst typed with max 3 missing loci, longer than 500kb
                    chromosomeclass=True
                    if indxb==0 and topspecies!='-' and int(length)>500000:
                        if int(num_missingloci)<=3:
                            if contigcompleteness=='complete':
                                contigclass='complete chromosome'
                            else:
                                contigclass='chromosome'
                        else:
                            chromosomeclass=False
                    else:
                        chromosomeclass=False
                    #classify non-chromosome contigs
                    if chromosomeclass==False:
                        if topspecies!='-':
                            if contigcompleteness=='complete':
                                contigclass='chromid'
                            else:
                                contigclass='chromosomal'
                        else:
                            if contigcompleteness=='complete':
                                contigclass='plasmid'
                            else:
                                contigclass='putative plasmid'
                                
                    ##write to file##
                    if typing=='both':
                        if reptype in reptypedict:
                            reptype=reptypedict[contig]
                        else:
                            reptype=('-','-','-','-','-','-')
                        f2.write('%s\t%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),'\t'.join(reptype),contigclass))
                    if typing=='rmlst':
                        f2.write('%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),contigclass))
                        
        ###iterate through contigs if sample groupings are not available
        else:
            for contig in sorted(contigs):
                length=seqlengthdict[contig]
                if contig in rmlsttypedict:
                    rmlsttype=rmlsttypedict[accession]
                else:
                    rmlsttype=('-','-','-','-','-','-','-')
                topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
                if contig in contigcompletenessdict:
                    contigcompleteness=contigcompletenessdict[contig]
                else:
                    contigcompleteness='unknown'
                #check for chromosome contig: rmlst typed with max 3 missing loci, longer than 500kb
                chromosomeclass=True
                if int(length)>500000 and topspecies!='-':
                    if int(num_missingloci)<=3:
                        if contigcompleteness=='complete':
                            contigclass='complete chromosome'
                        else:
                            contigclass='chromosome'
                    else:
                        chromosomeclass=False
                else:
                    chromosomeclass=False
                #classify non-chromosome contigs
                if chromosomeclass==False:
                    if topspecies!='-':
                        if contigcompleteness=='complete':
                            contigclass='chromid'
                        else:
                            contigclass='chromosomal'
                    else:
                        if contigcompleteness=='complete':
                            contigclass='plasmid'
                        else:
                            contigclass='putative plasmid'
                
                ##write to file##
                if typing=='both':
                    if contig in reptypedict:
                        reptype=reptypedict[contig]
                    else:
                        reptype=('-','-','-','-','-','-')
                    f2.write('%s\t%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),'\t'.join(reptype),contigclass))
                if typing=='rmlst':
                    f2.write('%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),contigclass))


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

