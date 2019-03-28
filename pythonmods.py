################Python modules

###wrapper for the subprocess command
def runsubprocess(args,stderrpath=None, stdoutpath=None, writefile=None,shell=False,verbose=True):
    import subprocess,sys #os
    try:
        import thread
    except:
        import _thread
    """takes a subprocess argument list and runs Popen/communicate(); by default, both output and error are printed to screen; stderrpath and stdoutpath for saving output can be optionally set; a redirect can be optionally set (writefile argument); errors are handled at multiple levels i.e. subthread error handling; can set shell=True; the function can be used 'fruitfully' since stdout is returned"""
    if shell==True: #e.g. args=['ls *.txt]
        processname=args[0] #ls *.txt
        processname=processname.split()#['ls', '*.txt'] #list argument syntax 
    else:
        processname=args
    processname=(" ".join(a for a in args))
    if stderrpath==None:
        pass
    else:
        if stderrpath.endswith('stderr.txt'): #want to make sure file ends with non-duplicated 'stderr.txt'
            stderrpath=str(stderrpath[:-10]).strip()
        stderrstrip=stderrpath.split('/')[-1]
        if stderrstrip=='': #there was nothing to strip after / i.e. was just /stderr.txt or stderr.txt
            pass
        else:
            stderrpath=stderrpath[:-(len(stderrstrip))]
        stderrpath=stderrpath+processname+'_'+stderrstrip+'stderr.txt'
    if stdoutpath==None:
        pass
    else:
        if stdoutpath.endswith('stdout.txt'): 
            stdoutpath=str(stdoutpath[:-10]).strip()
        stdoutstrip=stdoutpath.split('/')[-1]
        if stdoutstrip=='': 
            pass
        else:
            stdoutpath=stdoutpath[:-(len(stdoutstrip))]
        stdoutpath=stdoutpath+processname+'_'+stdoutstrip+'stdout.txt'
    if verbose==True:
        print('{} {}'.format(processname, 'processname'))
    try:
        if writefile==None:
            if shell==False:
                p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if shell==True:
                p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr= p.communicate()
            if verbose==True:
                try:
                    print('{} {}'.format(stdout.decode(), 'stdout'))
                except:
                    pass
                try:
                    print('{} {}'.format(stderr.decode(), 'stderr'))
                except:
                    pass
            if stdoutpath==None:
                pass
            else:
                with open(stdoutpath,'w') as stdoutfile:
                    stdoutfile.write(stdout)
            if stderrpath==None:
                pass
            else:
                with open(stderrpath,'w') as stderrfile:
                    stderrfile.write(stderr)
        else:
            with open(writefile,'w') as stdout:
                if shell==False:
                    p=subprocess.Popen(args,stdout=stdout, stderr=subprocess.PIPE)
                if shell==True:
                    p=subprocess.Popen(args,stdout=stdout, stderr=subprocess.PIPE, shell=True)
                stdout, stderr= p.communicate()
                if verbose==True:
                    try:
                        print('{} {}'.format(stdout.decode(), 'stdout'))
                    except:
                        pass
                    try:
                        print('{} {}'.format(stderr.decode(), 'stderr'))
                    except:
                        pass
                #n.b stdout is None - can't write to file
                if stderrpath==None:
                    pass
                else:
                    with open(stderrpath,'w') as stderrfile:
                        stderrfile.write(stderr)
        if p.returncode==0:
            if verbose==True:
                print('{} {}'.format(processname, 'code has run successfully'))
        else:
            if verbose==False:
                print('{} {}'.format(processname, 'processname'))
            print('{} {}'.format('source code fail'))
    except:
        if verbose==False:
            print('{} {}'.format(processname, 'processname'))
        print('runsubprocess code fail')
        sys.exit()
        
    if p.returncode!=0:
        if 'driverscript' in str(sys.argv[0]):
            sys.exit()
        else:
            #os._exit
            try:
                thread.interrupt_main()
            except:
                _thread.interrupt_main()
            #keyboard interrupt is more drastic than sys.exit - only use for subscript subprocess errors which wouldn't otherwise trigger driverscript abort
    else:
        return stdout



###wrapper for blast command
def runblastn(query, database, blastoutput, evalue=str(10), outfmt='custom qcov', task='blastn', num_threads=str(1), max_target_seqs=str(500), max_hsps=False, perc_identity=False, qcov_hsp_perc=False, culling_limit=False, word_size=False): #default evalue is 10; default wordsize is 11 for blastn - just use -tasks parameter which also changes gap settings); default task with no -task parameter set is megablast; N.b use cutstom qcov as standard for blast-based typing and plasmidpipeline, otherwise use outfmt 6 or specify custom outfmt 6                                                                      
    import subprocess
    evalue=str(evalue) #in case I forget to enter as string                                                                                                                                                 
    if outfmt=='custom qcov': #allow alias for custom output flags                                                                                                                                          
        outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen'
    cmdArgs=['blastn', '-query', query, '-db', database, '-out', blastoutput, '-evalue', evalue, '-outfmt', outfmt, '-task', task, '-num_threads', num_threads, '-max_target_seqs', max_target_seqs]
    if max_hsps!=False:
        cmdArgs.extend(['-max_hsps', max_hsps])
    if perc_identity!=False:
        cmdArgs.extend(['-perc_identity', perc_identity])
    if qcov_hsp_perc!=False:
        cmdArgs.extend(['-qcov_hsp_perc', qcov_hsp_perc])
    if culling_limit!=False:
        cmdArgs.extend(['-culling_limit', culling_limit])
    if word_size!=False:
        cmdArgs.extend(['-word_size', word_size])
    subprocess.call(cmdArgs)



def unlist(listed, d=','):
    """takes a list and converts to delimiter-separated string; comma delimited by default"""
    unlisted=(d.join(a for a in listed))
    return unlisted


###selecting best blast hits
def blastfilter(blastoutput,finalfile,sortedfile, idtable=False, pidthresh=90,coveragethresh=0.9,overlapthresh=0.5, longerqueries=True, keepblastoutput=False, keepsortedfile=False):
    """takes a blast output table and filters for best non-overalpping hits; if provided, idtable filepath will be used to rename subject sequence ids post filtering"""
    from pythonmods import runsubprocess, unlist

    #use blastfilter.sh to filter and sort blastoutput (outputs sortedfile)                                                                                                                                 
    formatcontigcol=False #hardcoded since this option only applies to mlst filtering                                                                                                                       
    args=['bash', 'blastfilter.sh','%s'%longerqueries,'%s'%formatcontigcol,'%s'%blastoutput,'%s'%sortedfile,'%s'%pidthresh,'%s'%coveragethresh]
    runsubprocess(args)

    #make dictionary based on score-sorted file of pid/coverage-filtered output, with query sequences as keys and score-ordered hits as nested list of values                                               
    uniquequerynames=[]
    querydict={}
    fileObj=open(sortedfile)
    for indx, line in enumerate(fileObj):
        data=line.strip().split('\t')
        queryname=data[0] #filtered qseqid                                                                                                                                                                  
        if queryname in querydict:
            querydict[queryname].append(data) #!previously did data[1:] but this is confusing                                                                                                               
        else:
            querydict[queryname]=[]
            querydict[queryname].append(data)
            uniquequerynames.append(queryname)
    fileObj.close()
    print("finished making querydict")

    #if idtable path is set, make oldname:newname iddict so oldname can be reassigned while writing to file                                                                                                 
    if idtable!=False:
        iddict={} #oldname : newname                                                                                                                                                                        
        with open(idtable) as f:
            for line in f:
                data=line.strip().split('\t')
                oldname=data[2]
                newname=data[1]
                if newname=='None':
                    iddict[oldname]=oldname
                else:
                    iddict[oldname]=newname

    #Reminder of file columns for outfmt custom qcov: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen                                          
    if longerqueries==True: #if the queries are the longer sequence then the 'gene' length is slen and hit overlaps should be considered on the query                                                       
        genelenindex=int(15) #slen                                                                                                                                                                          
        startindex=int(6) #qstart                                                                                                                                                                           
        endindex=int(7) #qend                                                                                                                                                                               
    else:
        genelenindex=int(14) #qlen                                                                                                                                                                          
        startindex=int(8) #sstart                                                                                                                                                                           
        endindex=int(9) #ssend

    #write best hits to file                                                                                                                                                                                
    fileObj=open(finalfile,'w')
    for indxa, query in enumerate(uniquequerynames):
        print('{} {}'.format(indxa, 'query indx, blastfilter'))
        hitrangegenelengths=[] #list of tuples of inlcuded hits and the associated gene lengths                                                                                                             
        for indxb, data in enumerate(querydict[query]): #running through all hits associated with a given query                                                                                             
            if indxb==0:  #for the highest scoring hit, include, irrespecitve of overlaps                                                                                                                   
                print('{} {} {} {} {}'.format(data, 'data = row of hit info', startindex, endindex, data[0]))
                ranges=[int(data[startindex]),int(data[endindex])]
                genelength=int(data[genelenindex])
                minrange=ranges[ranges.index(min(ranges))]
                maxrange=ranges[ranges.index(max(ranges))]
                hitrange=range(minrange,maxrange+int(1))   #this is a list of cosecutive number spanning the hit range                                                                                      
                hitrangegenelengths.append((hitrange, genelength))
                if idtable!=False:
                    oldname=data[1]
                    data[1]=iddict[oldname]
                fileObj.write('%s\n'%('\t'.join(data)))
            else:
                ranges=[int(data[startindex]),int(data[endindex])]
                genelength=int(data[genelenindex])
                minrange=ranges[ranges.index(min(ranges))]
                maxrange=ranges[ranges.index(max(ranges))]
                hitrange=range(minrange,maxrange+int(1))

                #if hit range intersects dont include otherwise include and add hitrange to ranges.                                                                                                         
                for indxc, includedhit in enumerate(hitrangegenelengths):  #checking non-best hit against all inlcuded hits; there will be at least one hit (the best hit) in hit ranges                    
                    includedhitrange=includedhit[0]
                    includedgenelength=includedhit[1]
                    genelengths=[genelength, includedgenelength]
                    mingenelength=genelengths[genelengths.index(min(genelengths))]
                    intersectlength=len(set(hitrange).intersection(set(includedhitrange)))
                    pairwiseoverlap=float(float(intersectlength)/float(mingenelength))  #floats are essential float(int/int) is 0.0                                                                         
                    if pairwiseoverlap<float(0.5):
                        add=True
                        continue

                    else:
                        add=False
                        break
                if add==True:
                    hitrangegenelengths.append((hitrange, genelength))
                    if idtable!=False:
                        oldname=data[1]
                        data[1]=iddict[oldname]
                    fileObj.write('%s\n'%('\t'.join(data)))
                else:
                    pass #if add is false i.e. there was an overlap with a better hit, don't write to file, don't append to hitrangegenelengths                                                             
    fileObj.close()
    #delete intermediate files                                                                                                                                                                              
    if keepblastoutput==False:
        runsubprocess(['rm %s'%blastoutput],shell=True)
    if keepsortedfile==False:
        runsubprocess(['rm %s'%sortedfile],shell=True)



###selecting best mlst blast hits

#N.B mlstfilter will select only the best allele per locus for each query; so even if there are multiple non-overlapping hits per locus only the best hit will be returned. This simplifies filtering/typing but multiple high-scoring non-overlapping hits could be used as indicator of contamination (on the other hand, perhaps sometimes chromosomes do have multiple rmlst/mlst loci with mosaic phylogenetic origins - so not sure whether multiple hits per locus would be a good marker of contamination except in extreme cases where there are multiple complete chromosomes of different species in a single sample [which would probably be picked up as contaminant anyway]).
def mlstfilter(blastoutput,finalfile,sortedfile, idtable=False, pidthresh=90,coveragethresh=0.9,overlapthresh=0.5, incF=False, longerqueries=True, formatcontigcol=True, keepblastoutput=False, keepsortedfile=False):
    """takes a blast output table of mlst alleles and filters for best alleles; similar code to blast filter (differences include: formatcontigcol can be set during filtering); formatcontigcol=True means do filtering on a per-sample basis (hence requiring the contig query column to be prefixed with sample-level column); formatcontigcol=False means do filtering on a per-contig level"""
    from mymod import runsubprocess, unlist

    #use blastfilter.sh to filter and sort blastoutput (outputs sortedfile)                                                                                                                                 
    args=['bash', 'blastfilter.sh','%s'%longerqueries,'%s'%formatcontigcol,'%s'%blastoutput,'%s'%sortedfile,'%s'%pidthresh,'%s'%coveragethresh]
    runsubprocess(args)

    #if idtable path is set, make oldname:newname iddict so oldname can be reassigned while writing to file                                                                                                 
    if idtable!=False:
        iddict={} #oldname : newname                                                                                                                                                                        
        with open(idtable) as f:
            for line in f:
                data=line.strip().split('\t')
                oldname=data[2]
                newname=data[1]
                if newname=='None':
                    iddict[oldname]=oldname
                else:
                    iddict[oldname]=newname

    #Reminder of file columns for outfmt custom qcov: (qseqid sample-level [if formatcontigcol=True]) qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen                                                                                                                                                                                                     
    #write to file the top hit allele for each locus                                                                                                                                                        

    if formatcontigcol==True:
        sseqindex=int(2)
    else:
        sseqindex=int(1)

    accessions=[]
    fileObj=open(sortedfile)
    fileObj2=open(finalfile,'w')
    for indx, line in enumerate(fileObj):
        data=line.strip().split('\t')
        accession=data[0]
        if formatcontigcol==True: ###added to handle contig output                                                                                                                                          
            allele=str(data[2].split("-")[0]) #for mlst allele in format allele-n (checked for ecoli/salmonella)                                                                                            
            allele=str(allele.rsplit('_',1)[0]) #for pmlst alleles in format allele_n or allele_x_n                                                                                                         
        else:
            allele=str(data[1].split("-")[0]) #for mlst allele in format allele-n (checked for ecoli/salmonella)                                                                                            
            allele=str(allele.rsplit('_',1)[0]) #for pmlst alleles in format allele_n or allele_x_n                                                                                                         
        if incF!=False:
            if allele.startswith('FIC') or allele.startswith('FII'):
                allele='F'

        if accession in accessions: #if the row represents an accession for which at least one locus is already represented by a top hit allele, only write to file if the row represents a locus that is not yet represented                                                                                                                                                                                          
            if allele in alleles: #i.e. top hit allele already written to file                                                                                                                              
                continue
            else:
                alleles.append(allele)
                if idtable!=False:
                    oldname=data[sseqindex]
                    data[sseqindex]=iddict[oldname]
                fileObj2.write('%s\n'%('\t'.join(data)))
        else: #if the row represents the top scoring hit for a given accession, write to file                                                                                                               
            alleles=[]
            accessions.append(accession)
            alleles.append(allele)
            if idtable!=False:
                oldname=data[sseqindex]
                data[sseqindex]=iddict[oldname]
            fileObj2.write('%s\n'%('\t'.join(data)))
    fileObj.close()
    fileObj2.close()
    #delete intermediate files                                                                                                                                                                              
    if keepblastoutput==False:
        runsubprocess(['rm %s'%blastoutput],shell=True)
    if keepsortedfile==False:
        runsubprocess(['rm %s'%sortedfile],shell=True)



###assigning replicon types

def inctypingprobes(inctype_probes,db='enterobacteriaceae'):
    """takes a list of probetypes; outputs inctypes, incfamilies, inclength"""
    from mymod import unlist
    import sys
    inclength=str(len(inctype_probes))
    nestinctypes=[] #inc type
    nestinctypes_concise=[] #inc family
    if db=='enterobacteriaceae':
        for idx, probe in enumerate(inctype_probes):
            probe=str(probe).strip()
            probesplit=probe.split('_')
            if probe.startswith('Col'):
                nestinctypes.append(probesplit[0]) #just including all col type plasmids as probe_ids
                nestinctypes_concise.append('Col')
            elif probe.startswith('IncA/C'):
                if probe.startswith('IncA/C_1'):
                    nestinctypes.append('IncA/C1')
                else:
                    nestinctypes.append(probesplit[0]) #IncA/C2
                nestinctypes_concise.append('IncA/C')
            elif probe.startswith('IncB/O/K/Z'):
                nestinctypes.append(str(probesplit[0])+str(probesplit[1]))
                nestinctypes_concise.append('IncB/O/K/Z')
            elif probe.startswith('IncF'):
                if probe.startswith('IncFII') or probe.startswith('IncFII('):
                    nestinctypes.append('IncFII')
                elif probe.startswith('IncFIA'):
                    nestinctypes.append('IncFIA')
                elif probe.startswith('IncFIB'):
                    nestinctypes.append('IncFIB')
                elif probe.startswith('IncFIC'):
                    nestinctypes.append('IncFIC')
                else:
                    print('{} {}'.format('IncF probe handling error', probe)); sys.exit()
                nestinctypes_concise.append('IncF')
            elif probe.startswith('IncHI'):
                if probe.startswith('IncHI1'):
                    nestinctypes.append('IncHI1')
                elif probe.startswith('IncHI2'):
                    nestinctypes.append('IncHI2')
                else:
                    print('{} {}'.format('IncHI probe handling error', probe)); sys.exit()
                nestinctypes_concise.append('IncH')
            elif probe.startswith('IncI'):
                if probe.startswith('IncI1'):
                    nestinctypes.append('IncI1')
                elif probe.startswith('IncI2'):
                    nestinctypes.append('IncI2')
                else:
                    print('{} {}'.format('IncI probe handling error', probe)); sys.exit()
                nestinctypes_concise.append('IncI')
            elif probe.startswith('IncL/M'):
                nestinctypes.append('IncL/M')
                nestinctypes_concise.append('IncL/M')
            elif probe.startswith('IncN'):
                if probe.startswith('IncN_1'):
                    nestinctypes.append('IncN1')
                else:
                    nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncN')
            elif probe.startswith('IncP') or probe.startswith('P1_alpha'):
                if probe.startswith('P1_alpha'):
                    nestinctypes.append('IncP-1alpha')
                elif probe.startswith('IncP(Beta)'):
                    nestinctypes.append('IncP-1beta')
                elif probe.startswith('IncP(6)'):
                    nestinctypes.append('IncP6')
                elif probe.startswith('IncP1'):
                    nestinctypes.append('IncP1')
                elif probe.startswith('IncP6'):
                    nestinctypes.append('IncP6')
                else:
                    print('{} {}'.format('IncP probe handling error', probe)); sys.exit()
                nestinctypes_concise.append('IncP')
            elif probe.startswith('IncQ'):
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncQ')
            elif probe.startswith('IncR'):
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncR')
            elif probe.startswith('IncT'):
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncT')
            elif probe.startswith('IncU'):
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncU')
            elif probe.startswith('IncW'):
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncW')
            elif probe.startswith('IncX'):
                if probe.startswith('IncX3('):
                    nestinctypes.append('IncX3')
                else:
                    nestinctypes.append(probesplit[0])
                nestinctypes_concise.append('IncX')
            elif probe.startswith('repA') or probe.startswith('Rep'): #these are recently added rep genes which don't fit with previous inc types                                                 
                if probe.startswith('Rep_1_pKPC-2_CP011573'):
                    nestinctypes.append('rep_CP011573')
                    nestinctypes_concise.append('rep_CP011573')
                elif probe.startswith('repA_1_pKPC-2_CP013325'):
                    nestinctypes.append('repA_CP013325')
                    nestinctypes_concise.append('repA_CP013325')
                elif probe.startswith('repA_2_pKPC-2_JX397875'):
                    nestinctypes.append('repA_JX397875')
                    nestinctypes_concise.append('repA_JX397875')
                elif probe.startswith('RepA_1_pKPC-CAV1321_CP011611'):
                    nestinctypes.append('repA_CP011611')
                    nestinctypes_concise.append('repA_CP011611')
                else:
                    print('{} {}'.format('rep/repA probe handling error', probe)); sys.exit()
            elif probe.startswith('FIA') or probe.startswith('FIA('):
                nestinctypes.append('IncFIA')
                nestinctypes_concise.append('IncF')
            elif probe.startswith('FII') or probe.startswith('FII('):
                nestinctypes.append('IncFII')
                nestinctypes_concise.append('IncF')
            else:
                nestinctypes.append(probesplit[0])
                nestinctypes_concise.append(probesplit[0])

        inctype_probes=sorted(inctype_probes) #sort in alphabetical order                                                                                                                         
        inctypes=sorted(nestinctypes)
        inctypes_concise=sorted(nestinctypes_concise)

        inctype_probes=str(unlist(inctype_probes))
        inctypes=str(unlist(inctypes))
        inctypes_concise=str(unlist(inctypes_concise))

        return (inctypes, inctypes_concise, inctype_probes, inclength)
    elif db=='gram_positive':
        for idx, probe in enumerate(inctype_probes):
            probe=str(probe).strip()
            probesplit=probe.split('_')
            nestinctypes.append(probesplit[0])
            nestinctypes_concise.append(probesplit[0])
            
        inctype_probes=sorted(inctype_probes) #sort in alphabetical order                                                                                                                         
        inctypes=sorted(nestinctypes)
        inctypes_concise=sorted(nestinctypes_concise)

        inctype_probes=str(unlist(inctype_probes))
        inctypes=str(unlist(inctypes))
        inctypes_concise=str(unlist(inctypes_concise))

        return (inctypes, inctypes_concise, inctype_probes, inclength)
    else:
        print('Error: unrecognised database')
        sys.exit()


        
###assigning rmlst types

def rmlstprofile():
    """extracts rmlst profile information; returns a tuple of lists allelenames (the allele column headers),rSTs,alleles (listed list of allele rows), genuses, species,clonal complexes; if cell is empty, returns 'unspecified' for that cell"""
    sts=[]
    alleles=[]
    genuses=[]
    species=[]
    ccs=[]
    with open('/well/bag/orlek/datasets/rmlstalleles/profiles.tsv') as f:
        for indx, line in enumerate(f):
            line=line.strip()
            data=line.split('\t')
            if indx==0:
                allelenames=data[1:-7]
            else:
                sts.append(data[0]) #STs       
                alleles.append(data[1:54]) #alleles list of lists
                try:
                    assert data[54]!=''
                    assert not data[54].isspace()
                    genuses.append(data[54])
                except:
                    genuses.append('unspecified')
                try:
                    assert data[55]!=''
                    assert not data[55].isspace()
                    species.append(data[55])
                except:
                    species.append('unspecified')
                try:
                    assert data[59]!=''
                    assert not data[59].isspace()
                    ccs.append(data[59])
                except:
                    ccs.append('unspecified')

    return (allelenames,sts,alleles,genuses,species,ccs)



def occurrencelist(mylist):
    """take a list and returns a string of element occurrences"""
    from collections import Counter
    occurrencetuple=Counter(mylist).most_common() #if nothing specified inside parentheses most_common returns element,occurrence tuple of all elements, ordered by most frequent first
    s=''
    for a,b in occurrencetuple:
        s+=str(a)+','+str(b)+'; '
    s=s.rstrip('; ')
    return(s)


def rmlsttypingalleles(rmlstalleles, missinglociallowed=True, mismatchesallowed=0, nonexactmismatchesallowed=5, numallelesfornonexactmatching=30):
    """takes a list of rmlst alleles and returns rST, species, rmlstalleles, clonal complex; if set of detected alleles do not match a recognised rST, outputs 'unknown'; but with missinglociallowed=True,
 N-based mismatches are ignored; if there are multiple different genuses/species/clonal complexes, an occurence list is returned; mismatchesallowed=0/missinglociallowed=True corresponds to the definition
 of an exact match on the rmlst webserver; if nonexactmismatchesallowed !=False and there are sufficient alleles (as specified by numallelesfornonexactmatching) then do second round of typing with relaxe
d mismatch threshold; non-exact output is prefixed with 'non-exact'"""
    from pythonmods import rmlstprofile,unlist,occurrencelist
    ###first extract profile information
    allelenames,sts,alleles,genuses,species,ccs=rmlstprofile()
    myalleledict={}
    for myallele in rmlstalleles:
        myallele=myallele.split('_')
        assert len(myallele)==2
        myallelename=myallele[0]
        myalleletype=myallele[1]
        myalleledict[myallelename]=myalleletype
    myalleles=[]
    for allelename in allelenames:
        try:
            myalleles.append(myalleledict[allelename])
        except:
            myalleles.append('N')
    assert len(myalleles)==53
    ###try to extract rST, genus, species, clonal complex without excluding N-based mismatching
    try:
        alleleindex=alleles.index(myalleles) #if not in list, this will raise error
        myrst=sts[alleleindex]
        mygenus=genuses[alleleindex]
        myspecies=species[alleleindex]
        mycc=ccs[alleleindex]
    except:
        #print('rmlstalleles do not match any set of alleles in the profile table exactly. Alleles: %s' %myalleles)
        myrst='unknown'
        mygenus='unknown'
        myspecies='unknown'
        mycc='unknown'
        ###but can you type if you exclude 'N' loci?...
        if missinglociallowed==True: #either an incomplete query genome or a false-positive (i.e. loci detected in query but recorded as N in a matching profile)
            print('checking taxonomy with Ns excluded')
            specieslist=[]
            genuslist=[]
            stlist=[]
            cclist=[]
            for indx, profilealleles in enumerate(alleles):
                mismatches=0
                for a,b in zip(myalleles,profilealleles):
                    if a=='N' or b=='N':
                        continue
                    if a==b:
                        continue
                    else:
                        mismatches=mismatches+1
                if mismatches<=int(mismatchesallowed):
                    specieslist.append(species[indx])
                    genuslist.append(genuses[indx])
                    stlist.append(sts[indx])
                    cclist.append(ccs[indx])
            specieslistset=sorted(list(set(specieslist)))
            genuslistset=sorted(list(set(genuslist)))
            cclistset=sorted(list(set(cclist)))
            if len(stlist)==0:
                print('no profile allele sets exactly match my alleles, even after excluding missing alleles')
                if nonexactmismatchesallowed!=False and len(myalleles)>numallelesfornonexactmatching:
                    ###recursive code to check for non-exact matching - redefining lists and retyping
                    print('checking for non-exact matching')
                    specieslist=[]
                    genuslist=[]
                    stlist=[]
                    cclist=[]
                    for indx, profilealleles in enumerate(alleles):
                        mismatches=0
                        for a,b in zip(myalleles,profilealleles):
                            if a=='N' or b=='N':
                                continue
                            if a==b:
                                continue
                            else:
                                mismatches=mismatches+1
                        if mismatches<=int(nonexactmismatchesallowed):
                            specieslist.append(species[indx])
                            genuslist.append(genuses[indx])
                            stlist.append(sts[indx])
                            cclist.append(ccs[indx])
                    specieslistset=sorted(list(set(specieslist)))
                    genuslistset=sorted(list(set(genuslist)))
                    cclistset=sorted(list(set(cclist)))
                    if len(stlist)==0:
                        print('no profile allele sets exactly match my alleles, even after allowing non-exact matching')
                    elif len(stlist)==1: #one profile allele set matches my alleles
                        mygenus='non-exact_'+genuslist[0]
                        myspecies='non-exact_'+specieslist[0]
                        myrst='non-exact_'+stlist[0]
                        mycc='non-exact_'+cclist[0]
                    else:  #multiple rSTs match my alleles; if there is a single genus/species/clonal complex output this, otherwise output occurrencelist
                        myrst='non-exact_'+unlist(sorted(stlist))
                        if len(genuslistset)>1:
                            mygenus='non-exact_'+occurrencelist(genuslist)
                        else:
                            mygenus='non-exact_'+genuslistset[0]
                        if len(specieslistset)>1:
                            myspecies='non-exact_'+occurrencelist(specieslist)
                        else:
                            myspecies='non-exact_'+specieslistset[0]
                        if len(cclistset)>1:
                            mycc='non-exact_'+occurrencelist(cclist)
                        else:
                            mycc='non-exact_'+cclistset[0]
                        ###end of recursive code
            elif len(stlist)==1: #one profile allele set matches my alleles
                mygenus=genuslist[0]
                myspecies=specieslist[0]
                myrst=stlist[0]
                mycc=cclist[0]
            else:  #multiple rSTs match my alleles; if there is a single genus/species/clonal complex output this, otherwise output occurrencelist 
                myrst=unlist(sorted(stlist))
                if len(genuslistset)>1:
                    mygenus=occurrencelist(genuslist)
                else:
                    mygenus=genuslistset[0]
                if len(specieslistset)>1:
                    myspecies=occurrencelist(specieslist)
                else:
                    myspecies=specieslistset[0]
                if len(cclistset)>1:
                    mycc=occurrencelist(cclist)
                else:
                    mycc=cclistset[0]
    return(myrst, mygenus, myspecies, mycc, unlist(sorted(rmlstalleles)))


