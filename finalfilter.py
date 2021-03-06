import sys
import numpy as np
from Bio import SeqIO
from pythonmods import inctyping, rmlstprofile,rmlsttypingalleles
from collections import defaultdict

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
    sampleoutput=sys.argv[9]
    typedcontigsonly=sys.argv[10]


def defaultreptype():
    return ['-','-','-','-','-','-'] #enterobacteriaceae types, probes, number of probes, gram-positive ...
def defaultrmlst():
    return ['-','-','-','-','-','-','-'] #topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles
def defaultcompleteness():
    return 'unknown'

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
                        enterobacaccessionsdict[accession]=[allele]
                else:
                    gram_posaccessions.append(accession)
                    if accession in gram_posaccessionsdict:
                        gram_posaccessionsdict[accession].append(allele)
                    else:
                        gram_posaccessionsdict[accession]=[allele]

    enterobacaccessions=list(set(enterobacaccessions))
    gram_posaccessions=list(set(gram_posaccessions))
    repaccessions=list(set(enterobacaccessions+gram_posaccessions))

    #do replicon typing
    reptypedict=defaultdict(defaultreptype)
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
                rmlstaccessionsdict[accession]=[allele]

    rmlstaccessions=list(set(rmlstaccessions))

    #do rmlst typing
    rmlstprofileoutput=rmlstprofile('%s/profiles.txt'%rmlstprofilepath)
    rmlsttypedict=defaultdict(defaultrmlst)
    for accession in rmlstaccessions:
        alleles=rmlstaccessionsdict[accession]
        rmlsttype=rmlsttypingalleles(alleles,rmlstprofileoutput)
        rmlsttypedict[accession]=rmlsttype #a tuple of top species, top rST, num_matches, num_mismatches, num_missing loci, num_multiallelic loci, rmlstalleles


###for inhouse, extract sequence lengths and contig completeness info if available
if sequenceorigin!='ncbi':
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
            
    contigcompletenessdict=defaultdict(defaultcompleteness)
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



###########output final files; output depends on whether using ncbi or inhouse sequences, and if the latter, whether rmlst has been conducted

#N.B nomenclature used below for sequences: "accession" if it's an NCBI sequence; "contig" if it's an inhouse sequence
                
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
                reptype=reptypedict[accession]
                ###classify accession
                if int(num_missingloci)<3 and int(length)>500000: #smallest bacterial genome is ~580Kb
                    accessionclass='chromosome'
                else:
                    accessionclass='chromosomal'
                #write to file
                f2.write('%s\t%s\t%s\t%s\n'%('\t'.join(data),'\t'.join(rmlsttype),'\t'.join(reptype),accessionclass))
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
                reptype=reptypedict[accession]
                f2.write('%s\t%s\n'%('\t'.join(data),'\t'.join(reptype)))
    f2.close()

else:
    ###output tsv with rmlst and/or replicon typing info + contig classification if rmlst is conducted
    #write column header for contig-level file
    f2=open('%s/contigtyping.tsv'%outdir,'w')
    header=['Contig','Length']
    if typing=='both' or typing=='rmlst':
        header.extend(['Species','ribosomalST','Num_Matches','Num_Mismatches','Num_Missing_Loci','Num_MultiallelicLoci','rMLST_alleles'])
    if typing=='both' or typing=='replicon':
        header.extend(['Enterobacteriaceae_Type','Probe_Hits','Num_Probe_hits','Gram_positive_Type','Probe_Hits','Num_Probe_hits'])
    if typing=='both' or typing=='rmlst':
        header.extend(['Contig_Classification'])
    header='\t'.join(header)
    f2.write('%s\n'%header)
    ###if only replicon typing is specified, write rep typing info to file, but don't attempt to classify contig; if sample groupings are specified, iterate through samples (and sample contigs); otherwise just iterate through contigs
    if contigsamplesfile!='None':
        #write column header for sample-level file
        if sampleoutput=='True':
            f3=open('%s/sampletyping.tsv'%outdir,'w')
            header=['Sample','Length']
            if typing=='both' or typing=='rmlst':
                header.extend(['Species','ribosomalST','Num_Matches','Num_Mismatches','Num_Missing_Loci','Num_MultiallelicLoci','rMLST_alleles'])
            if typing=='both' or typing=='replicon':
                header.extend(['Enterobacteriaceae_Type','Probe_Hits','Num_Probe_hits','Gram_positive_Type','Probe_Hits','Num_Probe_hits'])
            header='\t'.join(header)
            f3.write('%s\n'%header)
    
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
                    samplecontigdict[sample]=[contig]
        #iterate through samples (first pass - annotate chromosome, chromosomal, chromid contigs; leave other contigs as unknown for now)
        contigclassesdict={}
        for sample in sorted(samplecontigdict.keys()):           
            samplecontigs=samplecontigdict[sample]
            lengths=[]
            for contig in samplecontigs:
                lengths.append(seqlengthdict[contig])
            indices=list(np.argsort(lengths))
            indices.reverse()
            samplecontigs=[samplecontigs[i] for i in indices] #sample contigs are now sorted in length order, starting with longest
            lengths=[lengths[i] for i in indices]
            #iterate through sample contigs, sorted by length
            for indxb,(contig,length) in enumerate(zip(samplecontigs,lengths)):
                #get typing info and contig completeness
                reptype=reptypedict[contig]
                enterobactypes,enterobacprobes,enterobacnum,grampostypes,gramposprobes,gramposnum=reptype
                rmlsttype=rmlsttypedict[contig]
                topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
                contigcompleteness=contigcompletenessdict[contig]
                ##classify contig
                #check for chromosome contig: longest contig, rmlst typed with max 3 missing loci, longer than 500kb
                contigclass='unknown'
                if indxb==0 and topspecies!='-' and int(length)>500000:
                    if int(num_missingloci)<=3:
                        if contigcompleteness=='complete':
                            contigclass='complete chromosome'
                        else:
                            contigclass='chromosome'
                #classify non-chromosome sample contigs
                if contigclass=='unknown':
                    if topspecies!='-':  #at least one rMLST locus
                        if contigcompleteness=='complete':
                            contigclass='complete chromid'
                        else:
                            contigclass='chromosomal'
                    #else:  #no rMLST (so no topspecies); leave contig classified as unknown; after iterating through all contigs again, reassign according to decision tree
                #store contig class
                contigclassesdict[contig]=contigclass

        #iterate through samples (second pass - assign contigs without rMLST locus)
        for sample in sorted(samplecontigdict.keys()):
            f2.write('>%s\n'%sample)  
            samplecontigs=samplecontigdict[sample]
            contigclasses=set()
            for contig in samplecontigs:
                contigclasses.add(contigclassesdict[contig])
            lengths=[]
            for contig in samplecontigs:
                lengths.append(seqlengthdict[contig])
            indices=list(np.argsort(lengths))
            indices.reverse()
            samplecontigs=[samplecontigs[i] for i in indices] #sample contigs are now sorted in length order, starting with longest
            lengths=[lengths[i] for i in indices]
            #initialise sample-level lists
            if sampleoutput=='True':  #sampleoutput=True by default in lastest version, so this if statement is now uneccessary
                samplermlstalleles=[]
                sampleenterobactypes=[]
                sampleenterobacprobes=[]
                samplegrampostypes=[]
                samplegramposprobes=[]
                samplelength=str(sum(lengths))
            #iterate through sample contigs, sorted by length
            for contig,length in zip(samplecontigs,lengths):
                #get typing info and contig completeness
                contigclass=contigclassesdict[contig]
                reptype=reptypedict[contig]
                enterobactypes,enterobacprobes,enterobacnum,grampostypes,gramposprobes,gramposnum=reptype
                rmlsttype=rmlsttypedict[contig]
                topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
                contigcompleteness=contigcompletenessdict[contig]
                #append typing info to sample-level lists
                if sampleoutput=='True':
                    if rmlstalleles!='-':
                        samplermlstalleles.extend(rmlstalleles.split('|'))
                    if enterobactypes!='-':
                        sampleenterobactypes.extend(enterobactypes.split(','))
                        sampleenterobacprobes.extend(enterobacprobes.split(','))
                    if grampostypes!='-':
                        samplegrampostypes.extend(grampostypes.split(','))
                        samplegramposprobes.extend(gramposprobes.split(','))
                ##classify contig
                if contigclass=='unknown':
                    if enterobactypes!='-' or grampostypes!='-':
                        contigclass='putative plasmid'
                        if contigcompleteness=='complete':
                            contigclass='complete plasmid'
                        else:
                            if ('chromosomal' not in contigclasses) and ('chromosome' not in contigclasses): #if other rMLST contigs are either complete chromosome or chromid
                                contigclass='plasmid'
                    else:
                        if ('chromosomal' not in contigclasses) and ('chromosome' not in contigclasses):
                            contigclass='putative plasmid'
                                
                ##write contig-level typing to file##
                if typedcontigsonly=='True':  #skip untyped contigs
                    if topspecies=='-' and enterobactypes=='-' and grampostypes=='-':
                        continue
                if typing=='both':
                    f2.write('%s\t%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),'\t'.join(reptype),contigclass))
                elif typing=='rmlst':
                    f2.write('%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),contigclass))
                else: #replicon typing only
                    f2.write('%s\t%s\t%s\n'%(contig,length,'\t'.join(reptype)))
                    
            ###write sample-level typing to file
            if sampleoutput=='True':
                if len(sampleenterobactypes)==0:
                    sampleenterobactypes='-'
                    sampleenterobacprobes='-'
                    sampleenterobacnum='-'
                else:
                    indices=list(np.argsort(sampleenterobactypes))
                    sampleenterobactypes=[sampleenterobactypes[i] for i in indices]
                    sampleenterobacprobes=[sampleenterobacprobes[i] for i in indices]
                    sampleenterobacnum=str(len(sampleenterobactypes))
                    sampleenterobactypes=','.join(sampleenterobactypes)
                    sampleenterobacprobes=','.join(sampleenterobacprobes)
                if len(samplegrampostypes)==0:
                    samplegrampostypes='-'
                    samplegramposprobes='-'
                    samplegramposnum='-'
                else:
                    indices=list(np.argsort(samplegrampostypes))
                    samplegrampostypes=[samplegrampostypes[i] for i in indices]
                    samplegramposprobes=[samplegramposprobes[i] for i in indices]
                    samplegramposnum=str(len(samplegrampostypes))
                    samplegrampostypes=','.join(samplegrampostypes)
                    samplegramposprobes=','.join(samplegramposprobes)
                samplereptype=[sampleenterobactypes,sampleenterobacprobes,sampleenterobacnum,samplegrampostypes,samplegramposprobes,samplegramposnum]
                if len(samplermlstalleles)==0:
                    samplermlsttype=('-','-','-','-','-','-','-')
                else:
                    samplermlsttype=rmlsttypingalleles(samplermlstalleles,rmlstprofileoutput)
                    
                if typing=='both':
                    f3.write('%s\t%s\t%s\t%s\n'%(sample,samplelength,'\t'.join(samplermlsttype),'\t'.join(samplereptype)))
                elif typing=='rmlst':
                    f3.write('%s\t%s\t%s\n'%(sample,samplelength,'\t'.join(samplermlsttype)))
                else: #replicon typing only
                    f3.write('%s\t%s\t%s\n'%(sample,samplelength,'\t'.join(samplereptype)))
        if sampleoutput=='True':
            f3.close()
        
    else: #no sample groupings, iterate through contigs
        for contig in sorted(contigs):
            length=seqlengthdict[contig]
            reptype=reptypedict[contig]
            enterobactypes,enterobacprobes,enterobacnum,grampostypes,gramposprobes,gramposnum=reptype
            rmlsttype=rmlsttypedict[contig]
            topspecies,toprst,num_matches,num_mismatches,num_missingloci,num_multiallelicloci,rmlstalleles=rmlsttype
            contigcompleteness=contigcompletenessdict[contig]
            #check for chromosome contig: rmlst typed with max 3 missing loci, longer than 500kb                                        
            contigclass='unknown'
            if int(length)>500000 and topspecies!='-':
                if int(num_missingloci)<=3:
                    if contigcompleteness=='complete':
                        contigclass='complete chromosome'
                    else:
                        contigclass='chromosome'
            #classify non-chromosome contigs                                                                                            
            if contigclass=='unknown':
                if topspecies!='-':
                    if contigcompleteness=='complete':
                        contigclass='complete chromid'
                    else:
                        contigclass='chromosomal'
                else:
                    if enterobactypes!='-' or grampostypes!='-':
                        contigclass='putative plasmid'
                        if contigcompleteness=='complete':
                            contigclass='complete plasmid'

            ##write contig-level typing to file##
            if typedcontigsonly=='True':  #skip untyped contigs
                if topspecies=='-' and enterobactypes=='-' and grampostypes=='-':
                    continue
                
            if typing=='both':
                f2.write('%s\t%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),'\t'.join(reptype),contigclass))
            elif typing=='rmlst':
                f2.write('%s\t%s\t%s\t%s\n'%(contig,length,'\t'.join(rmlsttype),contigclass))
            else: #replicon typing only
                f2.write('%s\t%s\t%s\n'%(contig,length,'\t'.join(reptype)))            

    f2.close()




