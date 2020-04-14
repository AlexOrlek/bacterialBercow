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

def runsubprocess(args,verbose=False,shell=False,polling=False):
    import subprocess,sys
    try:
        import thread
    except:
        import _thread
    """takes a subprocess argument list and runs Popen/communicate or Popen/poll() (if polling=True); if verbose=True, processname (string giving command call) is printed to screen (processname is always printed if a process results in error); errors are handled at multiple levels i.e. subthread error handling"""
    if shell==True:
        processname=args[0]
        processname=processname[0].split()
        processname=(" ".join(a for a in processname))
    else:
        processname=(" ".join(a for a in args))
    if verbose==True:
        print('{0} {1}'.format(processname, 'processname'))
    try:
        if polling==True:
            p=subprocess.Popen(args, stdout=subprocess.PIPE,shell=shell)
            while True:
                stdout=p.stdout.readline()
                if p.poll() is not None:
                    break
                if stdout: #if stdout not empty...
                    print('{0}'.format(stdout.decode().strip()))
        else:
            p=subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell)
            stdout, stderr= p.communicate()
            if stdout:
                print('{0}'.format(stdout.decode()))
            if stderr:
                try: #want to output to stderr stream
                    if (sys.version_info > (3, 0)):
                        print('{0}'.format(stderr.decode()),file=sys.stderr) #Python3
                    else:
                        print>>sys.stderr,stderr  #Python2
                except: #if above code block fails for some reason, print stderr (to stdout)
                    print('{0}'.format(stderr.decode()))

        if p.returncode==0:
            if verbose==True:
                print('{0} {1}'.format(processname, 'code has run successfully'))
        else:
            sys.exit() #triggers except below
    except:
        print('{0} {1}'.format(processname, '#this pipeline step produced error'))
        print('unexpected error; exiting')
        sys.exit()
        
    if p.returncode!=0:
        print('unexpected error; exiting')
        try:
            thread.interrupt_main()
        except:
            _thread.interrupt_main()



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
    for line in fileObj:
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
    for query in uniquequerynames:
        #print('{0} {1}'.format(indxa, 'query indx, blastfilter'))
        hitrangegenelengths=[] #list of tuples of inlcuded hits and the associated gene lengths     
        for indx, data in enumerate(querydict[query]): #running through all hits associated with a given query
            if indx==0:  #for the highest scoring hit, include, irrespecitve of overlaps                                                                                             
                #print('{0} {1} {2} {3} {4}'.format(data, 'data = row of hit info', startindex, endindex, data[0]))
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
                for includedhit in hitrangegenelengths:  #checking non-best hit against all inlcuded hits; there will be at least one hit (the best hit) in hit ranges
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
    for line in fileObj:
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
        for probe in inctype_probes:
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
                    print('{0} {1}'.format('IncA/C probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncA/C')
            elif probe.startswith('IncB/O/K/Z'):
                if re.match('IncB/O/K/Z_[0-9]+_',probe):               
                    nestinctypes.append(str(probesplit[0])+str(probesplit[1]))
                else:
                    nestinctypes.append(probesplit[0])
                    print('{0} {1}'.format('IncIncB/O/K/Z probe could not be assigned a known replicon type', probe))
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
                    print('{0} {1}'.format('IncF probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncF')
            elif probe.startswith('IncHI'):
                if probe.startswith('IncHI1'):
                    nestinctypes.append('IncHI1')
                elif probe.startswith('IncHI2'):
                    nestinctypes.append('IncHI2')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{0} {1}'.format('IncHI probe could not be assigned a known replicon type', probe))
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
                    print('{0} {1}'.format('IncI probe could not be assigned a known replicon type', probe))
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
                    print('{0} {1}'.format('IncN probe could not be assigned a known replicon type', probe))
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
                    print('{0} {1}'.format('IncP probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncP')
            elif probe.startswith('IncQ'):
                if probe.startswith('IncQ1'):
                    nestinctypes.append('IncQ1')
                elif probe.startswith('IncQ2'):
                    nestinctypes.append('IncQ2')
                else:
                    nestinctypes.append(probesplit[0])
                    print('{0} {1}'.format('IncQ probe could not be assigned a known replicon type', probe))
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
                    print('{0} {1}'.format('IncX probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncX')
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
                    print('{0} {1}'.format('IncY probe could not be assigned a known replicon type', probe))
                nestinctypes_concise.append('IncY')
            else:
                print('{0} {1}'.format('unknown Inc probe could not be assigned a known replicon type/family', probe)) #includes recently added rep genes which don't fit with previous inc types e.g. repA_CP011611
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
        for probe in inctype_probes:
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
    """takes a list of rmlst alleles and returns top matching species, top matching rST(s), num_matches, num_mismatches, num_missingloci; rMLST alleles detected. There may be multiple top matching species/rST (if there are equally good matches). If top species encompass multiple genera, then top genera are given. num_matches includes matches to allele 'N' in top rST profile i.e. cases where allele is missing in sample and in profile. num_mismtches doesn't include  missing loci which are recorded as N in top rST profile. num_missingloci doesn't include loci recorded as 'N' in profile"""
    import operator
    ###first extract profile information, and get list of my alleles in the same order as profiles, filling missing loci with 'N'
    loci,sts,profiles,genuses,species=rmlstprofileoutput
    myalleledict={}
    for myallele in rmlstalleles:
        myallele=myallele.split('_')
        assert len(myallele)==2
        mylocus=myallele[0]
        myallele=myallele[1]
        if mylocus in myalleledict:
            myalleledict[mylocus].append(myallele)
        else:
            myalleledict[mylocus]=[myallele]
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
            #if there is a contaminant locus (multiple alleles), record, but still see if one allele matches
            if len(myallele)>1:
                for mymultipleallele in myallele:
                    if mymultipleallele==profileallele:
                        matches.append(mymultipleallele)
            else:
                if myallele[0]==profileallele:
                    matches.append(myallele[0])
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
        if y=='N':
            continue
        if len(x)>1: #multple alleles detected, so it cannot = N, so don't need to add to missingloci counter; either at least one allele matches or none match in which case it's a mismatch
            matchfound=False
            for mymultipleallele in x:
                if mymultipleallele==y and mymultipleallele!='N' and matchfound==False: #this excludes N==N-based matches
                    topmatchcounter=topmatchcounter+1
                    matchfound=True
            if matchfound==False:
                topmismatchcounter=topmismatchcounter+1
        else:    
            if x[0]==y and x[0]!='N': #this excludes N==N-based matches
                topmatchcounter=topmatchcounter+1  
            elif x[0]=='N':
                topmissinglocicounter=topmissinglocicounter+1
            else:
                topmismatchcounter=topmismatchcounter+1

    ###get multi-allelic loci (potential contaminants)
    contaminantloci=[]
    for indx,myallele in enumerate(myalleles):
        if topprofile[indx]=='N' or myallele=='N':
            continue
        if len(myallele)>1:
            contaminantloci.append({loci[indx]:myallele})

    return(str(species),str(toprst),str(topmatchcounter),str(topmismatchcounter),str(topmissinglocicounter),str(len(contaminantloci)),str('|'.join(rmlstalleles))) #top species, top rST, num_matches, num_mismatches, num_missingloci, num_multiallelicloci,alleles


