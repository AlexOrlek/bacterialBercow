################Python modules

###miscellaneous
def unique(mylist):
    seen = set()
    return [x for x in mylist if not (x in seen or seen.add(x))]


def unlist(listed, d=','):
    """takes a list and converts to delimiter-separated string; comma delimited by default"""
    unlisted=(d.join(a for a in listed))
    return unlisted


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



###selecting best blast hits
def blastfilter(blastoutput,finalfile,sortedfile,sourcedir,idtable=False, pidthresh=90,coveragethresh=0.9,overlapthresh=0.5,longerqueries=True, keepblastoutput=False, keepsortedfile=False):
    """takes a blast output table and filters for best non-overalpping hits; if provided, idtable filepath will be used to rename subject sequence ids post filtering"""
    #use blastfilter.sh to filter and sort blastoutput (outputs sortedfile)
    formatcontigcol=False #hardcoded since this option only applies to mlst filtering
    args=['bash', '%s/blastfilter.sh'%sourcedir,'%s'%longerqueries,'%s'%formatcontigcol,'%s'%blastoutput,'%s'%sortedfile,'%s'%pidthresh,'%s'%coveragethresh]
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
        #print('{} {}'.format(indxa, 'query indx, blastfilter'))
        hitrangegenelengths=[] #list of tuples of inlcuded hits and the associated gene lengths     
        for indxb, data in enumerate(querydict[query]): #running through all hits associated with a given query
            if indxb==0:  #for the highest scoring hit, include, irrespecitve of overlaps                                                                                             
                #print('{} {} {} {} {}'.format(data, 'data = row of hit info', startindex, endindex, data[0]))
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
def mlstfilter(blastoutput,finalfile,sortedfile,sourcedir,idtable=False, pidthresh=90,coveragethresh=0.9,overlapthresh=0.5, incF=False, longerqueries=True, formatcontigcol=True, keepblastoutput=False, keepsortedfile=False):
    """takes a blast output table of mlst alleles and filters for best alleles; similar code to blast filter (differences include: formatcontigcol can be set during filtering); formatcontigcol=True means do filtering on a per-sample basis (hence requiring the contig query column to be prefixed with sample-level column); formatcontigcol=False means do filtering on a per-contig level"""
    #from pythonmods import runsubprocess, unlist
    #use blastfilter.sh to filter and sort blastoutput (outputs sortedfile)
    args=['bash', '%s/blastfilter.sh'%sourcedir,'%s'%longerqueries,'%s'%formatcontigcol,'%s'%blastoutput,'%s'%sortedfile,'%s'%pidthresh,'%s'%coveragethresh]
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
    import sys,re
    inclength=str(len(inctype_probes))
    nestinctypes=[] #inc type
    nestinctypes_concise=[] #inc family
    if db=='enterobacteriaceae':
        for idx, probe in enumerate(inctype_probes):
            probe=str(probe).strip()
            probesplit=probe.split('_')
            if probe.startswith('Col'):
                nestinctypes.append(probesplit[0]) #just including all col type plasmids as probe_id prefixes
                nestinctypes_concise.append('Col')
            elif probe.startswith('IncA/C'):
                if probe.startswith('IncA/C_1'):
                    nestinctypes.append('IncA/C1')
                elif probe.startswith('IncA/C2_1'):
                    nestinctypes.append('IncA/C2')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncA/C probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncA/C')
            elif probe.startswith('IncB/O/K/Z'):
                if re.match('IncB/O/K/Z_[0-9]+_',probe):               
                    nestinctypes.append(str(probesplit[0])+str(probesplit[1]))
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncIncB/O/K/Z probe could not be assigned a known replicon type', probe))
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
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncF probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncF')
            elif probe.startswith('IncHI'):
                if probe.startswith('IncHI1'):
                    nestinctypes.append('IncHI1')
                elif probe.startswith('IncHI2'):
                    nestinctypes.append('IncHI2')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncHI probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncH')
            elif probe.startswith('IncI'):
                if probe.startswith('IncI1'):
                    nestinctypes.append('IncI1')
                elif probe.startswith('IncI2'):
                    nestinctypes.append('IncI2')
                elif probe.startswith('IncI_Gamma_1') or probe.startswith('IncI_1_Gamma') or probe.startswith('IncI_gamma_1') or probe.startswith('IncI_1_gamma'):
                    nestinctypes.append('IncI1-gamma')
                elif probe.startswith('IncI_Alpha_1') or probe.startswith('IncI_1_Alpha') or probe.startswith('IncI_alpha_1') or probe.startswith('IncI_1_alpha'):
                    nestinctypes.append('IncI1-alpha')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncI probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncI')
            elif probe.startswith('IncL/M'):
                nestinctypes.append('IncL/M') #the IncL/M family isn't subdivided into replicon types in plasmidfinder 
                nestinctypes_concise.append('IncL/M')
            elif probe.startswith('IncN'):
                if probe.startswith('IncN_1') or probe.startswith('IncN1'):
                    nestinctypes.append('IncN1')
                elif probe.startswith('IncN2'):
                    nestinctypes.append('IncN2')
                elif probe.startswith('IncN3'):
                    nestinctypes.append('IncN3')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncN probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncN')
            elif probe.startswith('IncP') or probe.startswith('P1_alpha'):
                if probe.startswith('P1_alpha'):
                    nestinctypes.append('IncP1-alpha')
                elif probe.startswith('IncP(Beta)'):
                    nestinctypes.append('IncP1-beta')
                elif probe.startswith('IncP(6)') or probe.startswith('IncP6'):
                    nestinctypes.append('IncP6')
                elif probe.startswith('IncP1'):
                    nestinctypes.append('IncP1')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncP probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncP')
            elif probe.startswith('IncQ'):
                if probe.startswith('IncQ1'):
                    nestinctypes.append('IncQ1')
                elif probe.startswith('IncQ2'):
                    nestinctypes.append('IncQ2')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncQ probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncQ')
            elif probe.startswith('IncR'):
                nestinctypes.append(probesplit[0]) #not subdivided into different replicon types in plasmidfinder
                nestinctypes_concise.append('IncR')
            elif probe.startswith('IncT'):
                nestinctypes.append(probesplit[0]) #not subdivided into different replicon types in plasmidfinder
                nestinctypes_concise.append('IncT')
            elif probe.startswith('IncU'):
                nestinctypes.append(probesplit[0]) #not subdivided into different replicon types in plasmidfinder
                nestinctypes_concise.append('IncU')
            elif probe.startswith('IncW'):
                nestinctypes.append(probesplit[0]) #not subdivided into different replicon types in plasmidfinder
                nestinctypes_concise.append('IncW')
            elif probe.startswith('IncX'):
                if probe.startswith('IncX3('):
                    nestinctypes.append('IncX3')
                elif re.match('IncX[0-9]+_',probe):
                    nestinctypes.append(probesplit[0])
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncX probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncX')
            # elif probe.startswith('repA') or probe.startswith('Rep'): #these are recently added rep genes which don't fit with previous inc types                                                 
            #     if probe.startswith('Rep_1_pKPC-2_CP011573'):
            #         nestinctypes.append('rep_CP011573')
            #         nestinctypes_concise.append('rep_CP011573')
            #     elif probe.startswith('repA_1_pKPC-2_CP013325'):
            #         nestinctypes.append('repA_CP013325')
            #         nestinctypes_concise.append('repA_CP013325')
            #     elif probe.startswith('repA_2_pKPC-2_JX397875'):
            #         nestinctypes.append('repA_JX397875')
            #         nestinctypes_concise.append('repA_JX397875')
            #     elif probe.startswith('RepA_1_pKPC-CAV1321_CP011611'):
            #         nestinctypes.append('repA_CP011611')
            #         nestinctypes_concise.append('repA_CP011611')
            #     else:
            #         nestinctypes.append(probesplit[0])
            #         print('{} {}'.format('Rep/repA probe could not be assigned a known replicon type', probe))
            elif probe.startswith('FIA') or probe.startswith('FIA('):
                nestinctypes.append('IncFIA')
                nestinctypes_concise.append('IncF')
            elif probe.startswith('FII') or probe.startswith('FII('):
                nestinctypes.append('IncFII')
                nestinctypes_concise.append('IncF')
            elif probe.startswith('IncY'):
                if probe.startswith('IncY_1'):
                    nestinctypes.append('IncY') #currently only 1 IncY type IncY_1
                else:
                    nestinctypes.append(probesplit[0])
                    print('{} {}'.format('IncY probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncY')
            else:
                print('{} {}'.format('unknown Inc probe could not be assigned a known replicon type/family', probe))
                nestinctypes.append(probe)
                nestinctypes_concise.append(probe)
                
        #sort probes,types,families in alphabetical order, according to probes
        inctype_probes_list=[]
        inctypes_list=[]
        inctypes_concise_list=[]
        for p,t,f in sorted(zip(inctype_probes,nestinctypes,nestinctypes_concise)):
            inctype_probes_list.append(p)
            inctypes_list.append(t)
            inctypes_concise_list.append(f)
        #convert sorted lists to strings
        inctype_probes=str(unlist(inctype_probes_list))
        inctypes=str(unlist(inctypes_list))
        inctypes_concise=str(unlist(inctypes_concise_list))

        return (inctypes, inctypes_concise, inctype_probes, inclength)
    elif db=='gram_positive':
        import re
        gramposregex=re.compile(r'((?:^.*\|)?)(.*)') #allows probes of syntax Familyprefix|reptype_probeinfo OR reptype_probeinfo, by capturing both parts separately 
        for idx, probe in enumerate(inctype_probes):
            probe=str(probe).strip()
            probegroup1=gramposregex.search(probe).group(1)
            probegroup2=gramposregex.search(probe).group(2)
            probesplit=probegroup2.split('_')
            nestinctypes.append('%s%s'%(probegroup1,probesplit[0]))
            nestinctypes_concise.append('%s%s'%(probegroup1,probesplit[0]))

        #sort probes,types,families in alphabetical order, according to probes
        inctype_probes_list=[]
        inctypes_list=[]
        inctypes_concise_list=[]
        for p,t,f in sorted(zip(inctype_probes,nestinctypes,nestinctypes_concise)):
            inctype_probes_list.append(p)
            inctypes_list.append(t)
            inctypes_concise_list.append(f)
        #convert sorted lists to strings
        inctype_probes=str(unlist(inctype_probes_list))
        inctypes=str(unlist(inctypes_list))
        inctypes_concise=str(unlist(inctypes_concise_list))

        return (inctypes, inctypes_concise, inctype_probes, inclength)
    else:
        print('Error: unrecognised database')
        sys.exit()




def inctyping(sample, blastdict, db='enterobacteriaceae'):
    """takes a dictionary of samples/accessions : plasmidfinder inctype probes; outputs inctypes, incfamilies, probetypes, number of replicons"""
    blastkeys=list(blastdict.keys())
    if sample in blastkeys:
        unsortedincprobes=blastdict[sample] #blastdict[sample] is a list of one or more inc types associated with a given sample            
        inctypes,incfamilies,incprobes,inclength=inctypingprobes(unsortedincprobes,db)
    else:
        incprobes="-"
        inctypes="-"
        incfamilies="-"
        inclength="-"
    return (inctypes, incfamilies, incprobes, inclength)


        
###assigning rmlst types

def rmlstprofile(profilepath):
    """extracts rmlst profile information; returns a tuple of lists: loci (the allele column headers),rSTs; and listed lists: alleles (listed list of allele rows); and lists: genuses, species, clonal complexes; if cell is empty, returns 'unspecified' for that cell"""
    sts=[]
    profiles=[]
    genuses=[]
    species=[]
    ccs=[]
    with open(profilepath) as f:  #profile.tsv file
        for indx, line in enumerate(f):
            data=line.strip().split('\t')
            if indx==0:
                lociindices=[indx for indx,i in enumerate(data) if i.startswith('BACT')]
                loci=data[lociindices[0]:lociindices[-1]+1]
            else:
                sts.append(data[0]) #STs       
                profiles.append(data[lociindices[0]:lociindices[-1]+1]) #profiles list of lists  #data[1:54]
                try:
                    #genuses.append(data[54])
                    genuses.append(data[lociindices[-1]+1])
                except:
                    genuses.append('unspecified')
                try:
                    #species.append(data[55])
                    species.append(data[lociindices[-1]+2])
                except:
                    species.append('unspecified')

    return (loci,sts,profiles,genuses,species)



def rmlsttypingalleles(rmlstalleles, rmlstprofileoutput):
    """takes a list of rmlst alleles and returns top matching species, rST, num_matches, num_mismatches, num_missingloci; information on other close matches; N.B num_matches/mismatches/missingloci is calculated after excluding 'N' alleles in the profile; but top matching species includes N-based matches (so if there are very few loci detected, the top match will be a profile with many N positions)"""
    import operator
    ###first extract profile information, and get list of my alleles in the same order as profiles, filling missing loci with 'N'
    loci,sts,profiles,genuses,species=rmlstprofileoutput
    myalleledict={}
    for myallele in rmlstalleles:
        myallele=myallele.split('_')
        assert len(myallele)==2
        mylocus=myallele[0]
        myallele=myallele[1]
        myalleledict[mylocus]=myallele
    myalleles=[] #alleles in order of locus column heading ie. corresponding to allele profile order
    for locus in loci:
        try:
            myalleles.append(myalleledict[locus])
        except:
            myalleles.append('N')
    #assert len(myalleles)==53
    
    ###get top match(es)
    nummatches=[]
    for profile in profiles:
        matches=[]
        for myallele,profileallele in zip(myalleles,profile):
            if myallele=='N':
                continue
            if myallele==profileallele:
                matches.append(myallele)
        nummatches.append(len(matches))

    topmatches=max(nummatches)
    topindices=[indx for indx,i in enumerate(nummatches) if i == topmatches]

    ###get top species if known
    toprsts=[]
    topspecies=[]
    topgenuses=[]
    for indx in topindices:
        toprsts.append(sts[indx])
        topspecies.append(species[indx])
        topgenuses.append(genuses[indx])
  
    if len(set(topspecies))==1:
        species=topspecies[0]
    else:
        if len(set(genuses))==1:
            species='unknown: possible species: %s'%';'.join(sorted(list(set(topspecies))))
        else:
            species='unknown: possible genera: %s'%';'.join(sorted(list(set(topgenuses))))
            
    ###get top rST
    if len(toprsts)==1:
        toprst=toprsts[0]
    else:
        toprst='%i top matches including %s'%(len(toprsts),toprsts[0]) #select first rST top match

    ###get rST match stats (if there are multiple rST top matches, select first rST)
    topprofile=profiles[topindices[0]]
    topmatchcounter=0
    topmismatchcounter=0
    topmissinglocicounter=0
    for x,y  in zip(myalleles, topprofile):
        if x==y:
            topmatchcounter=topmatchcounter+1  #this include N==N-based matches
        elif x=='N':
            topmissinglocicounter=topmissinglocicounter+1
        else:
            topmismatchcounter=topmismatchcounter+1
        
    return(str(species),str(toprst),str(topmatchcounter),str(topmismatchcounter),str(topmissinglocicounter),str('|'.join(rmlstalleles)))


