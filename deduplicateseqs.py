import sys, operator
from Bio import SeqIO
from pythonmods import unique
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

deduplicationmethod=sys.argv[1]
outdir=sys.argv[2]


###########create list of plasmid accessions; find deduplicated accessions (choose the refseq duplicate over genbank duplciate), then write to file, excluding 'bad' entries


###find deduplicated accessions

#first download whole plasmid sequences as well as accession ids and number of features


#get accession info from accessions_filtereed_biosample.tsv
biosamplenamedict={} #biosample accession : name of submittor
with open('%s/accessions_filtered_metadata.tsv'%outdir) as f:
    for indx, line in enumerate(f):
        if indx==0:
            continue
        data=line.strip().split('\t')
        biosampleaccession=data[0]
        name=data[1]
        owner=data[2]
        biosamplenamedict[biosampleaccession]=[name,owner]  #duplicate data (due to nucleotide accessions with same biosample accession id) will be overwritten

#N.B if the nucleotide accession had no biosample elink, this will just result in no output row in accessions_filtered_metadata

#print biosamplenamedict["SAMN05717715"]
#sys.exit()
nuclaccessions=[]  #accessions_filtered_dblinks.tsv should have same number of rows (nucleotide accessions) as accessions_filtered.tsv
metadatadict={} #accession : [biosample accession, name]
bioprojectdict={}
bioprojectdict_noaccessionversion={}
with open('%s/accessions_filtered_dblinks.tsv'%outdir) as f:
    for indx, line in enumerate(f):
        if indx==0:
            continue
        data=line.strip().split('\t')
        accession=data[0]
        biosample=data[1]
        bioproject=data[2]
        if int(bioproject)==0: #projectid 0 means no bioproject link
            bioproject='-'
        nuclaccessions.append(accession)
        if biosample in biosamplenamedict:
            name,owner=biosamplenamedict[biosample]
            if name=='- -':
                name='-'
        else:
            name='-'
            owner='-'
        metadatadict[accession]=[biosample,name,owner]
        bioprojectdict[accession]=bioproject
        accession_noversion=accession.split('.')[0]
        bioprojectdict_noaccessionversion[accession_noversion]=bioproject

#print metadatadict["NZ_CP017408.1"]
#sys.exit()


#get dictionary of refseq:genbank accession ids and use to edit accession:bioprojectid dict
refseqgenbankdict={}
with open('%s/accessions_filtered_refseq_gb.tsv'%outdir) as f:
    for indx, line in enumerate(f):
        if indx==0:
            continue
        data=line.strip().split('\t')
        refseqaccession=data[0]
        genbankaccession=data[1]
        refseqgenbankdict[refseqaccession]=genbankaccession #remember that refseq accession is in "accession.version" format whereas genbank is in "accession" format



for accession in nuclaccessions:
    #if accession is refseq, edit the project id in bioprojectdict
    if '_' in accession: #refseq accession
        gbaccession=refseqgenbankdict[accession]
        if gbaccession in bioprojectdict_noaccessionversion:
            gbbioproject=bioprojectdict_noaccessionversion[gbaccession]
            bioprojectdict[accession]=gbbioproject
        else:  #no cognate genbank accession - overwriting with "-" (alternative would be to keep refseq bioproject to allow for refseq only input, but in this case probably best to use owner metadata)
            bioprojectdict[accession]='-'
            
#for accession in sorted(nuclaccessions):
#    print accession,bioprojectdict[accession]
#sys.exit()

###testing - if all 3 metadata items are changed for one accession, both accessions written to fasta
#metadatadict["NZ_CP026758.1"][0]='testSAMN'
#metadatadict["NZ_CP026758.1"][1]='testname'
#metadatadict["NZ_CP026758.1"][2]='testowner'
###

#determine identical by using sequence as dictionary key; then select single deduplicated accessions to write to file based on accession type (refseq/genbank) and biosample

plasmiddict={}

fileObj=open('%s/accessions_filtered.fa'%outdir)
for indx, seq_record in enumerate(SeqIO.parse(fileObj,"fasta")):
    accession=str(seq_record.id)
    if '_' in accession:
        accessiontype='refseq'
    else:
        accessiontype='genbank'
    sequence=str(seq_record.seq)
    biosample,name,owner=metadatadict[accession]
    bioproject=bioprojectdict[accession]
    if sequence in plasmiddict:
        plasmiddict[sequence].append([accession,accessiontype,biosample,name,owner,bioproject])
    else:
        plasmiddict[sequence]=[]
        plasmiddict[sequence].append([accession,accessiontype,biosample,name,owner,bioproject])






#write to file
deduplicatebymetadata=False #deduplicate all identical accessions
f2=open('%s/accessions_filtered_deduplicated.fa'%outdir,'w')
f3=open('%s/identicalaccessions.tsv'%outdir,'w')
f3.write('Accessions\tRefseq_Accessions\tGenbank_Accessions\tBiosamples\tSubmitter_Names\tOwners\tBioProjects\tUnique_Biosamples\tUnique_Names\tUnique_Owners\tUnique_BioProjects\tNum_Accessions\tNum_Refseq_Accessions\tNum_Genbank_Accessions\tNum_Unique_Biosamples\tNum_Unique_Names\tNum_Unique_Owners\tNum_Unique_BioProjects\tSame_Biosample\tSame_Submitter_Name\tSame_Owner\tSame_BioProject\tAccessions_Deduplicatedby_Metadata\tAccessions_Deduplicatedby_BioProject\n')
for seq, values in plasmiddict.iteritems():
    if len(values)>1: #multiple accessions for a given sequence
        accessions=map(operator.itemgetter(0), values)
        accessiontypes=map(operator.itemgetter(1), values)
        biosamples=map(operator.itemgetter(2), values)
        names=map(operator.itemgetter(3), values)
        owners=map(operator.itemgetter(4), values)
        bioprojects=map(operator.itemgetter(5), values)
        #group accessions by metadata - if they have either the same biosample accession or the same name (ignoring missing info), they are assigned the same group
        #problem: may have different biosample and name not provided - this would lead to different group but maybe should be more leniant and require both to be different? Yes - name / lab group info is the right idea - just need to make sure it's there (lab group would be better?)
        #reason for both requiring to be different: could have same lab group but copy pasting sequence with different biosample names; could have same biosample being sequenced / used by different lab groups.
        biosamplelist=[]
        namelist=[]
        ownerlist=[]
        bioprojectlist=[]
        metadatagroups=[]
        bioprojectgroups=[]
        for biosample,name,owner,bioproject in zip(biosamples,names,owners,bioprojects):
            if biosample not in biosamplelist and biosample!='-':
                biosamplelist.append(biosample)
            if name not in namelist and name!='-':
                namelist.append(name)
            if owner not in ownerlist and owner!='-':
                ownerlist.append(owner)
            if bioproject not in bioprojectlist and bioproject!='-':
                bioprojectlist.append(bioproject)
            if biosample=='-':
                biosampleindx='NA'
            else:
                biosampleindx=[indx for indx,b in enumerate(biosamplelist) if b==biosample][0]
            if name=='-':
                nameindx='NA'
            else:
                nameindx=[indx for indx,n in enumerate(namelist) if n==name][0]
            if owner=='-':
                ownerindx='NA'
            else:
                ownerindx=[indx for indx,o in enumerate(ownerlist) if o==owner][0]
            if bioproject=='-':
                bioprojectindx='NA'
            else:
                bioprojectindx=[indx for indx,p in enumerate(bioprojectlist) if p==bioproject][0]
            #print biosampleindx,nameindx
            metadataindices=[]
            bioprojectindices=[]
            if biosampleindx!='NA':
                metadataindices.append(biosampleindx)
                bioprojectindices.append(biosampleindx)
            if nameindx!='NA':
                metadataindices.append(nameindx)
            if ownerindx!='NA':
                metadataindices.append(ownerindx)
            if bioprojectindx!='NA':
                bioprojectindices.append(bioprojectindx)
            if len(metadataindices)==3: #i.e. all 3 fields present
                metadatagroup=min(metadataindices)
            else:
                metadatagroup=0
            metadatagroups.append(metadatagroup)
            if len(bioprojectindices)==2: #i.e. all 2 fields present
                bioprojectgroup=min(bioprojectindices)
            else:
                bioprojectgroup=0
            bioprojectgroups.append(bioprojectgroup)
        metadatagroupsset=set(metadatagroups)
        bioprojectgroupsset=set(bioprojectgroups)
        #print metadatagroupsset, biosamples,names,len(metadatagroupsset)
        refseqaccessions=[a for a in accessions if '_' in a]
        genbankaccessions=[a for a in accessions if '_' not in a]
        if len(refseqaccessions)==0:
            refseqaccessions='-'
        if len(genbankaccessions)==0:
            genbankaccessions='-'
        biosamplesunique=unique(biosamples)
        namesunique=unique(names)
        ownersunique=unique(owners)
        bioprojectsunique=unique(bioprojects)
        biosamplesuniquenonmissing=[i for i in biosamplesunique if i!='-']
        namesuniquenonmissing=[i for i in namesunique if i!='-']
        ownersuniquenonmissing=[i for i in ownersunique if i!='-']
        bioprojectsuniquenonmissing=[i for i in bioprojectsunique if i!='-']
        if len(biosamplesuniquenonmissing)==len(biosamplesunique):
            if len(biosamplesuniquenonmissing)==1:
                samebiosample='TRUE'
            else: #more than one unique biosample
                samebiosample='FALSE'
        else:
            samebiosample='missing data'
            
        if len(namesuniquenonmissing)==len(namesunique):
            if len(namesuniquenonmissing)==1:
                samename='TRUE'
            else:
                samename='FALSE'
        else:
            samename='missing data'

        if len(ownersuniquenonmissing)==len(ownersunique):
            if len(ownersuniquenonmissing)==1:
                sameowner='TRUE'
            else:
                sameowner='FALSE'
        else:
            sameowner='missing data'

        if len(bioprojectsuniquenonmissing)==len(bioprojectsunique):
            if len(bioprojectsuniquenonmissing)==1:
                samebioproject='TRUE'
            else:
                samebioproject='FALSE'
        else:
            samebioproject='missing data'

        #f3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('|'.join(accessions),'|'.join(refseqaccessions),'|'.join(genbankaccessions),'|'.join(biosamples),'|'.join(names),'|'.join(owners),'|'.join(bioprojects),'|'.join(biosamplesunique),'|'.join(namesunique),'|'.join(ownersunique),'|'.join(bioprojectsunique),len(accessions),len(refseqaccessions),len(genbankaccessions),len(biosamplesuniquenonmissing),len(namesuniquenonmissing),len(ownersuniquenonmissing),len(bioprojectsuniquenonmissing),samebiosample,samename,sameowner,samebioproject))
        #if len(metadatagroupsset)>1:
        #    f3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(','.join(accessions),','.join(refseqaccessions),','.join(genbankaccessions),','.join(biosamples),','.join(names),','.join(biosamplesunique),','.join(namesunique),len(accessions),len(refseqaccessions),len(genbankaccessions),len(biosamplesuniquenonmissing),len(namesuniquenonmissing),'different biosample and name'))
        #else:
        #    f3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(','.join(accessions),','.join(refseqaccessions),','.join(genbankaccessions),','.join(biosamples),','.join(names),','.join(biosamplesunique),','.join(namesunique),len(accessions),len(refseqaccessions),len(genbankaccessions),len(biosamplesuniquenonmissing),len(namesuniquenonmissing),'same biosample or name'))
        metadatadeduplicatedaccessions=[]
        bioprojectdeduplicatedaccessions=[]
        for group in metadatagroupsset:
            indices=[indx for indx,g in enumerate(metadatagroups) if g==group] #get all indices that equal group
            if len(indices)>1: #multiple accessions for a given metadata group - deduplicate by selecting refseq if available
                myrefseqindices=[indx for indx,t in enumerate(accessiontypes) if indx in indices and t=='refseq']
                if len(myrefseqindices)>0: #if there's at least 1 refseq, select the first refseq accession
                    accession=accessions[myrefseqindices[0]]
                else: #if there's only genbank, select the first genbank
                    accession=accessions[indices[0]]
            else: #1 accession for a given metadata group; select accession
                accession=accessions[indices[0]]
            metadatadeduplicatedaccessions.append(accession)
            if deduplicationmethod=="submitter": #choose one sequence per submitter metadata group
                f2.write('>%s\n'%accession)
                f2.write('%s\n'%seq)
        for group in bioprojectgroupsset:
            indices=[indx for indx,g in enumerate(bioprojectgroups) if g==group] #get all indices that equal group
            if len(indices)>1: #multiple accessions for a given bioproject group - deduplicate by selecting refseq if available
                myrefseqindices=[indx for indx,t in enumerate(accessiontypes) if indx in indices and t=='refseq']
                if len(myrefseqindices)>0: #if there's at least 1 refseq, select the first refseq accession
                    accession=accessions[myrefseqindices[0]]
                else: #if there's only genbank, select the first genbank
                    accession=accessions[indices[0]]
            else: #1 accession for a given bioproject group; select accession
                accession=accessions[indices[0]]
            bioprojectdeduplicatedaccessions.append(accession)
            if deduplicationmethod=="bioproject": #choose one sequence per bioproject group
                f2.write('>%s\n'%accession)
                f2.write('%s\n'%seq)
                
        if deduplicationmethod=="all": #choose one sequence across all duplicates irrespective of whether or not metadata/bioproject is shared
            myrefseqindices=[indx for indx,t in enumerate(accessiontypes) if t=='refseq']
            if len(myrefseqindices)>0: #if there's at least 1 refseq, select the first refseq accession
                accession=accessions[myrefseqindices[0]]
            else: #if there's only genbank, select the first genbank
                accession=accessions[0]
            f2.write('>%s\n'%accession)
            f2.write('%s\n'%seq)

        f3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('|'.join(accessions),'|'.join(refseqaccessions),'|'.join(genbankaccessions),'|'.join(biosamples),'|'.join(names),'|'.join(owners),'|'.join(bioprojects),'|'.join(biosamplesunique),'|'.join(namesunique),'|'.join(ownersunique),'|'.join(bioprojectsunique),len(accessions),len(refseqaccessions),len(genbankaccessions),len(biosamplesuniquenonmissing),len(namesuniquenonmissing),len(ownersuniquenonmissing),len(bioprojectsuniquenonmissing),samebiosample,samename,sameowner,samebioproject,'|'.join(metadatadeduplicatedaccessions),'|'.join(bioprojectdeduplicatedaccessions)))
    else: #1 accession for a given sequence
        accession=values[0][0]
        f2.write('>%s\n'%accession)
        f2.write('%s\n'%seq)

f2.close()
f3.close()


#OLD CODE - USING BIOSAMPLE AND CREATEDATE, AND TREATING "-" AS A DISTINCT BIOSAMPLE        


#metadatadict["NC_001735.4"][1]='alexorlekedit'
#metadatadict["U67194.4"][1]='alexedit'
#metadatadict["NC_001735.4"][0]='SAMNedit'
#metadatadict["U67194.4"][0]='SAMNedit2'

###testing - the following changes mean duplicate sequences are written to file for acessions LR134277.1 and NZ_LR134277.1
#metadatadict["NC_001735.4"][1]='2016/09/29'
#metadatadict["LR134277.1"][0]='testSAMN'
#metadatadict["LR134277.1"][1]='2019/12/19'
###


#write to file
# f2=open('%s/accessions_filtered_deduplicated.fa'%outdir,'w')
# f3=open('%s/duplicateaccessions.tsv'%outdir,'w')
# f3.write('Accessions\tUnique_Biosamples\tUnique_Dates\tNum_Accessions\tIs_Biosample_Same\n')
# for seq, values in plasmiddict.iteritems():
#     if len(values)>1: #multiple accessions for a given sequence
#         accessions=map(operator.itemgetter(0), values)
#         accessiontypes=map(operator.itemgetter(1), values)
#         biosamples=map(operator.itemgetter(2), values)
#         dates=map(operator.itemgetter(3), values)
#         #group accessions by metadata - if they have either the same biosample accession or the same date, they are assigned the same group
#         biosamplelist=[]
#         datelist=[]
#         groups=[]
#         for biosample,date in zip(biosamples,dates):
#             if biosample not in biosamplelist:
#                 biosamplelist.append(biosample)
#             if date not in datelist:
#                 datelist.append(date)
#             biosampleindx=[indx for indx,b in enumerate(biosamplelist) if b==biosample]
#             dateindx=[indx for indx,d in enumerate(datelist) if d==date]
#             groupindx=min(biosampleindx,dateindx)[0]
#             groups.append(groupindx)
#         groupsset=set(groups)
#         #print groupsset, biosamples,dates,len(groupsset)
#         if len(groupsset)>1:
#             f3.write('%s\t%s\t%s\t%s\t%s\n'%(','.join(accessions),','.join(unique(biosamples)),','.join(unique(dates)),len(accessions),'different biosample or date'))
#         else:
#             f3.write('%s\t%s\t%s\t%s\t%s\n'%(','.join(accessions),','.join(unique(biosamples)),','.join(unique(dates)),len(accessions),'same biosample or date'))
#         for group in groupsset:
#             indices=[indx for indx,g in enumerate(groups) if g==group] #get all indices that equal group
#             if len(indices)>1: #multiple accessions for a given metadata group - deduplicate by selecting refseq if available
#                 myrefseqindices=[indx for indx,t in enumerate(accessiontypes) if indx in indices and t=='refseq']
#                 if len(myrefseqindices)>0: #if there's at least 1 refseq, select the first refseq accession
#                     accession=accessions[myrefseqindices[0]]
#                 else: #if there's only genbank, select the first genbank
#                     accession=accessions[indices[0]]
#             else: #1 accession for a given metadata group; select accession
#                 accession=accessions[indices[0]]
#             f2.write('>%s\n'%accession)
#             f2.write('%s\n'%seq)
#     else: #1 accession for a given sequence
#         accession=values[0][0]
#         f2.write('>%s\n'%accession)
#         f2.write('%s\n'%seq)

# f2.close()
# f3.close()






#OLD CODE

# f2=open('%s/accessions_filtered_deduplicated.fa'%outdir,'w')
# f3=open('%s/accessions_filtered_duplicate.tsv'%outdir,'w')
# f3.write('Accessions\tNum_Accessions\tIs_Biosample_Same\n')
# for seq, values in plasmiddict.iteritems():
#     if len(values)>1: #multiple accessions for a given sequence
#         accessions=map(operator.itemgetter(0), values)
#         accessiontypes=map(operator.itemgetter(1), values)
#         biosamples=map(operator.itemgetter(2), values)
#         dates=map(operator.itemgetter(3), values)
#         biosampleset=set(biosamples)
#         if len(biosampleset)>1:
#             f3.write('%s\t%s\t%s\n'%(','.join(accessions),len(accessions),'sample biosample'))
#         else:
#             f3.write('%s\t%s\t%s\n'%(','.join(accessions),len(accessions),'different biosample'))
#         for biosample in biosampleset:
#             indices=[indx for indx,b in enumerate(biosamples) if b==biosample] #get all indices that equal biosample set item
#             if len(indices)>1: #multiple accessions for a given biosample - deduplicate by selecting refseq if available
#                 myrefseqindices=[indx for indx,t in enumerate(accessiontypes) if indx in indices and t=='refseq']
#                 if len(myrefseqindices)>0: #if there's at least 1 refseq, select the first refseq accession
#                     accession=accessions[myrefseqindices[0]]
#                 else: #if there's only genbank, select the first genbank
#                     accession=accessions[indices[0]]
#             else: #1 accession for a given biosample; select accession
#                 accession=accessions[indices[0]]
#     else: #1 accession for a given sequence
#         accession=values[0][0]
#     f2.write('>%s\n'%accession)
#     f2.write('%s\n'%seq)

# f2.close()
# f3.close()




#OLD CODE
        #currentaccession,currentaccessiontype,currentbiosample=plasmiddict[sequence]
        #print(currentaccession,accession)
        #print(currentaccessiontype,accessiontype)
        #print(currentbiosample,biosample)
        #if currentbiosample!=biosample:
        #    print('interesting!')
        #    sys.exit()


# #find deduplicated accessions - based on favouring refseqs, if no refseq then favour most richly annotated accessions; see old file for previous code (favouring most richly annotated)

# #create a dictionary with unique sequences as keys and accessions, no. features, and date as values

# plasmid_dict={}
# plasmid_file = open("%s/2nd_task/output/%s/%splasmidseqdownload.tsv" %(filepath,sys.argv[1],sys.argv[1]))

# for i, line in enumerate(plasmid_file):
#     data=line.split('\t')
#     key=data[6].strip()
#     if key in plasmid_dict:
#         plasmid_dict[key].append([data[0].strip(),data[5].strip(), data[7].strip()])   #added date to dictionary
#     else:
#         plasmid_dict[key]=[]
#         plasmid_dict[key].append([data[0].strip(),data[5].strip(), data[7].strip()])

# print len(plasmid_dict), "len plasmid dict"



# #find best unique accessions, and the earliest date for the unique sequence across duplicate acccessions

# #approach: if there is only one accesions, append to list; else: make list of duplicate accessions/no. features/date; find earliest date; find refseq accession, if there is a single refseq, append to list; elif there are no refseqs: find most feature rich accession; elif there are multiple refseqs: pick first one

# #UPDATE: no longer getting 'earliest dates' 


# fileObj = open("%s/2nd_task/output/%s/%sdeduplicatesequencecheck.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")  #shows details of accessions that have duplicate seuqences
# from datetime import datetime

# #earliestdates={}
# deduplicatedaccessions=[]
# for indx, (key, values) in enumerate(plasmid_dict.iteritems()):
#     if len(values)>1:
#         fileObj.write('%s\n' % key[0:50])
#         features=[]
#         accessions=[]
#         bestaccessionindex=[]
#         dates=[]
#         for value in values:
#             fileObj.write('%s\n' %value)
#             features.append(value[1])
#             accessions.append(value[0])
#             date=str(value[2])
#             dates.append(datetime.strptime(date, '%d-%b-%Y'))
            
#         #find earliest date
#         #earliestdate=str(min(dates).strftime('%d-%b-%Y'))

#         #find best accession
#         refseqaccessionindex=[(a,b) for a, b in zip(accessions,features) if '_' in a]
#         if len(refseqaccessionindex)==1:
#             deduplicatedaccessions.append(refseqaccessionindex[0][0]) ###if there is one refseq accession
#             #earliestdates[str(refseqaccessionindex[0][0])]=[earliestdate] 
#         else:
#             bestaccessionindex=[(a,b) for a, b in zip(accessions,features) if b == max(features)] #if there are no refseq accessions (or multiple refseqs), select from genbank accessions
#             if len(refseqaccessionindex)==0:
#                 deduplicatedaccessions.append(bestaccessionindex[0][0]) ###select genbank accession with most number of features, or if there are tied max feature accessions, select the first accession
#                 #earliestdates[str(bestaccessionindex[0][0])]=[earliestdate] 
#             else:
#                 deduplicatedaccessions.append(refseqaccessionindex[0][0]); print "duplicate refseqs"; print refseqaccessionindex ###if there are multiple refseq accessions, select only first refseq accession
#                 #earliestdates[str(refseqaccessionindex[0][0])]=[earliestdate] ####
#     else:
#         deduplicatedaccessions.append(values[0][0])  ###if there is only one accession
#         #earliestdates[str(values[0][0])]=[values[0][2]]


# print len(deduplicatedaccessions), "len deduplicated accessions (after removing accessions with 100% sequence identity"
# #print len(earliestdates), "len earliest dates"
# fileObj.close()




# #WRITE TO FILE

# #n.b. changed to v1 (except for excluded duplicateaccessions) since these files are futher filtered by text mining (flterbadaccessions.py)

# excludedids=[]
# fileObj = open("%s/2nd_task/downloads/%s/%sdownload_%s.gb" %(filepath,sys.argv[1],sys.argv[1],downloaddate))
# fileObj2 = open("%s/2nd_task/output/%s/%sfeaturesdownload_v1.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w") #contains plasmid column
# fileObj3 = open("%s/2nd_task/output/%s/%snucleotideseqdownload_v1.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w") #also contains protein sequences so no need for separate code in extractproteins directory writing enterobactranslationdownload.tsv (as done previously)
# fileObj4 = open("%s/2nd_task/output/%s/%srecorddownload_v1.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# fileObj5 = open("%s/2nd_task/output/%s/%sdeduplicateplasmidsdownload_v1.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# fileObj6 = open("%s/2nd_task/output/%s/%sduplicateaccessions.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w") 


# for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
#     print j
#     if rec.id not in list(notcdsplasmids): 
#         if rec.id in list(deduplicatedaccessions): #deduplicated accessions, based on sequence id
#             record_id = rec.id
#             #name = rec.name #can tell from rec.id
#             plasmid=plasmids[j] 
#             record_description = rec.description
#             record_organism = rec.annotations['organism']
#             record_length = len(rec.seq)
#             record_sequence = rec.seq
#             record_date= rec.annotations['date']  #added so can look at trends by time
#             #record_earliestdate=unlist(earliestdates[str(record_id)])  #not useful - just the earliest update date
#             record_features = len(rec.features)
#             record_dbxrefs=rec.dbxrefs    #"The dbxrefs list gets populated from any PROJECT or DBLINK lines" i.e. biosample and bioproject accesions form the DBLINK line are retrieved
#             if not record_dbxrefs: #empty list is false
#                 record_dbxrefs="-"
#             else:
#                 record_dbxrefs=unlist(record_dbxrefs)
#             featurecounter=0
#             for feature in rec.features:
#                 if feature.type == "CDS": #there will be some features since featureless excluded
#                     featurecounter=featurecounter+1
#                     strand = feature.location.strand  #1 is positive strand, -1 is negative strand
#                     startlocation = feature.location.start
#                     endlocation = feature.location.end
#                     #try:  #db_xref is GI id which is being phased out
#                     #    db_xref=unlist(feature.qualifiers["db_xref"])
#                     #except KeyError:
#                     #    db_xref="-"
#                     try:
#                         note=unlist(feature.qualifiers["note"])
#                     except KeyError:
#                         note="-"
#                     try:
#                         gene=unlist(feature.qualifiers["gene"])
#                     except KeyError:
#                         gene="-"
#                     try:
#                         locus_tag=unlist(feature.qualifiers["locus_tag"])
#                     except KeyError:
#                         locus_tag="-"
#                     try:
#                         product=unlist(feature.qualifiers["product"])
#                     except KeyError:
#                         product="-"
#                     try:
#                         protein_id=unlist(feature.qualifiers["protein_id"])
#                     except KeyError:
#                         protein_id="-"
#                     try: 
#                         translation=unlist(feature.qualifiers["translation"]) #this is important in order to exclude nucleotide sequences that don't code for protein
#                     except KeyError:
#                         translation="-"
#                     try:
#                         nucleotideseq=feature.extract(rec.seq) #feature.location.extract(rec).seq *also works; rec.seq[startlocation:endlocation] *works but may need reverse complementation  
#                     except KeyError:
#                         nucleotideseq="-" 

#                     fileObj2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_organism, plasmid, strand, startlocation, endlocation, gene, product, protein_id, locus_tag, record_dbxrefs, record_length, note))  #removed db_xref (GI ID - being phased out); added plasmid length - quick dirty way to flag dodgy start-end locations
#                     #fileObj3.write('%s|%s|%s\t%s\t%s\n' % (protein_id, product, record_organism, nucleotideseq, translation))#translation
#                     fileObj3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, protein_id, gene, product, record_organism, nucleotideseq, translation)) #May 2017: added gene
#                 record_features_CDS=featurecounter
#             fileObj4.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_description, record_organism, plasmid, record_length, record_features, record_features_CDS, record_date, record_dbxrefs))
#             fileObj5.write('%s\t%s\n' % (record_id, record_sequence))
#         else:

#             record_id = rec.id
#             record_description = rec.description
#             record_organism = rec.annotations['organism']

#             for i, feature in enumerate(rec.features):
#                 if i<1:
#                     record_plasmid=unlist(feature.qualifiers['plasmid'])
#             record_length = len(rec.seq)
#             record_features= len(rec.features)
#             featurecounter=0
#             for feature in rec.features:
#                 if feature.type == "CDS":
#                     featurecounter=featurecounter+1
#             record_features_CDS=featurecounter
#             record_date= rec.annotations['date']
#             fileObj6.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_description, record_organism, record_plasmid, record_length, record_features, record_features_CDS, record_date))  #duplicate accessions not selected as best duplicate
#     else:
#         pass


# fileObj.close()
# fileObj2.close()
# fileObj3.close()
# fileObj4.close()
# fileObj5.close()
# fileObj6.close()



# #make enterobacplasmidseqdownload.tsv readable by removing plasmid sequence column

# fileObj = open("%s/2nd_task/output/%s/%splasmidseqdownload.tsv" %(filepath,sys.argv[1],sys.argv[1]))
# fileObj2 = open("%s/2nd_task/output/%s/%splasmidseqdownload_readable.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")

# for i, line in enumerate(fileObj):   
#     data=line.split('\t')  
#     fileObj2.write('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (data[0], data[1], data[2], data[3], data[4], data[5]))
    


# #record_id, record_description, record_organism, record_length, record_features, record_features_CDS, record_sequence



# #############################

# #DESCRIPTIVE STATISTICS - PRE-FILTERING

# #post-deduplication descriptive statistics

# #count the number of occurences of each species
# fileObj=open("%s/2nd_task/output/%s/%srecorddownload_v1.tsv" %(filepath,sys.argv[1],sys.argv[1]))
# organisms=[]

# for line in fileObj:
#     data = line.split('\t')
#     organisms.append(data[2])
# fileObj.close()

# species=[]
# for organism in organisms:
#     data=organism.split()
#     data2=data[0]+' '+data[1]
#     species.append(data2)

# #count occurences
# from collections import Counter
# organismsdict = Counter(organisms)
# organismslist=organismsdict.items()

# fileObj2=open("%s/2nd_task/output/%s/%sorganismcounts.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# for i in range(0, len(organismslist)):
#     organismnames=organismslist[i][0]
#     organismcounts=organismslist[i][1]
#     fileObj2.write('%s\t%s\n' % (organismnames, organismcounts))

# #count occurences
# from collections import Counter
# speciesdict = Counter(species)
# specieslist=speciesdict.items()

# fileObj3=open("%s/2nd_task/output/%s/%sspeciescounts.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# for i in range(0, len(specieslist)):
#     speciesnames=specieslist[i][0]
#     speciescounts=specieslist[i][1]
#     fileObj3.write('%s\t%s\n' % (speciesnames, speciescounts))



# #pre-deduplication descriptive statistics

# #count the number of occurences of each species
# fileObj=open("%s/2nd_task/output/%s/%splasmidseqdownload.tsv" %(filepath,sys.argv[1],sys.argv[1]))
# organisms=[]

# for line in fileObj:
#     data = line.split('\t')
#     organisms.append(data[2])
# fileObj.close()

# species=[]
# for organism in organisms:
#     data=organism.split()
#     data2=data[0]+' '+data[1]
#     species.append(data2)

# #count occurences
# from collections import Counter
# organismsdict = Counter(organisms)
# organismslist=organismsdict.items()

# fileObj2=open("%s/2nd_task/output/%s/%sorganismcounts_prededuplication.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# for i in range(0, len(organismslist)):
#     organismnames=organismslist[i][0]
#     organismcounts=organismslist[i][1]
#     fileObj2.write('%s\t%s\n' % (organismnames, organismcounts))

# #count occurences
# from collections import Counter
# speciesdict = Counter(species)
# specieslist=speciesdict.items()

# fileObj3=open("%s/2nd_task/output/%s/%sspeciescounts_prededuplication.tsv" %(filepath,sys.argv[1],sys.argv[1]), "w")
# for i in range(0, len(specieslist)):
#     speciesnames=specieslist[i][0]
#     speciescounts=specieslist[i][1]
#     fileObj3.write('%s\t%s\n' % (speciesnames, speciescounts))








    
########################OLD CODE


#check whether record is plasmid by looking at feature source - code prints all that aren't referenced to plasmid

#!!not really necessary - this should be the same as the plasmid filter (except in rare case X78316.1 which would be excluded anyway)
"""

notplasmidids=[]
fileObj = open("%s/2nd_task/downloads/%s/enterobacdownload_%s.gb" %(filepath,downloaddate))

for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    for i, feature in enumerate(rec.features):
        if i<1:
            try:
                print feature.qualifiers['plasmid']
            except:
                notplasmidids.append(rec.id)
        else:
            break #only need to check for plasmid annotation on one feature


print notplasmidids, "notplasmidids"
fileObj.close()

savepickle(notplasmidids, "%s/2nd_task/pickle.notplasmidids" %filepath)
"""




#write excluded accessions to file
#JUST DONE FOR FEATURELESS

"""
picklefileObj=open("%s/2nd_task/pickle.featurelessids" %filepath, "rb")
featurelessids=pickle.load(picklefileObj)
picklefileObj2=open("%s/2nd_task/pickle.notplasmidids" %filepath, "rb")
notplasmidids=pickle.load(picklefileObj2)
combined=list(featurelessids)+list(notplasmidids)
notcdsplasmids=set(combined)

picklefileObj.close()
picklefileObj2.close()


fileObj = open("%s/2nd_task/downloads/%s/enterobacdownload_%s.gb" %(filepath,downloaddate))
fileObj1 = open("%s/2nd_task/enterobacfeatureless.tsv" %filepath, "w") 
fileObj2 = open("%s/2nd_task/enterobacnotplasmids.tsv" %filepath, "w") 
fileObj3 = open("%s/2nd_task/enterobacnotcdsplasmids.tsv" %filepath, "w") #changed from featurelessexcluded 

for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    print j
    if rec.id in list(featurelessids):
        record_id1 = rec.id
        record_description1 = rec.description
        record_organism1 = rec.annotations['organism']
        record_length1 = len(rec.seq)
        fileObj1.write('%s\t%s\t%s\t%s\n' % (record_id1, record_organism1, record_length1, record_description1)) 
    if rec.id in list(notplasmidids):
        record_id2 = rec.id
        record_description2 = rec.description
        record_organism2 = rec.annotations['organism']
        record_length2 = len(rec.seq)
        fileObj2.write('%s\t%s\t%s\t%s\n' % (record_id2, record_organism2, record_length2, record_description2))
    if rec.id in list(notcdsplasmids):
        record_id3 = rec.id
        record_description3 = rec.description
        record_organism3 = rec.annotations['organism']
        for i, feature in enumerate(rec.features):
            if i<1:
                try:
                    record_plasmid3=unlist(feature.qualifiers['plasmid'])
                except: #this doesn't seem to be necessary - unnamed seems to be the default
                    record_plasmid3=" "
        record_length3 = len(rec.seq)
        record_features3= len(rec.features)
        featurecounter3=0
        for feature in rec.features:
            if feature.type == "CDS":
                featurecounter3=featurecounter3+1
        record_features_CDS3=featurecounter3
        record_date3= rec.annotations['date']
        fileObj3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id3, record_description3, record_organism3, record_plasmid3, record_length3, record_features3, record_features_CDS3, record_date3))  
    else:
        continue



fileObj.close()
fileObj1.close()
fileObj2.close()
fileObj3.close()
"""

####

"""
picklefileObj2=open("%s/2nd_task/pickle.notplasmidids" %filepath, "rb")
notplasmidids=pickle.load(picklefileObj2)
combined=list(featurelessids)+list(notplasmidids)
notcdsplasmids=set(combined) #changed from deduplicated

picklefileObj2.close()

print "set of ids to exclude created"
"""



#PICKLING - unecessary
"""
picklefileObj=open("%s/2nd_task/output/%s/pickle.%sfeaturelessids" %(filepath,sys.argv[1],sys.argv[1]), "rb")
featurelessids=pickle.load(picklefileObj)
picklefileObj.close()
"""

#savepickle(featurelessids, "%s/2nd_task/output/%s/pickle.%sfeaturelessids" %(filepath,sys.argv[1],sys.argv[1]))

#savepickle(plasmids, "%s/2nd_task/pickle.plasmids" %filepath)
#savepickle(plasmid_dict, "%s/2nd_task/output/%s/pickle.%splasmid_dict" %(filepath,sys.argv[1],sys.argv[1]))

#picklefileObj=open("%s/2nd_task/output/%s/pickle.%splasmid_dict" %(filepath,sys.argv[1],sys.argv[1]), "rb")
#plasmid_dict=pickle.load(picklefileObj)

#savepickle(deduplicatedaccessions, "%s/2nd_task/pickle.deduplicatedaccessions" %filepath)
#savepickle(earliestdates, "%s/2nd_task/pickle.earliestdates_dict" %filepath)

#picklefileObj=open("%s/2nd_task/pickle.plasmids" %filepath, "rb")
#plasmids=pickle.load(picklefileObj)
#print "opened pickle plasmid"

#picklefileObj1=open("%s/2nd_task/pickle.earliestdates_dict" %filepath, "rb")
#earliestdates=pickle.load(picklefileObj1)
#print "opened pickle deduplicatedaccessions"

#picklefileObj2=open("%s/2nd_task/pickle.deduplicatedaccessions" %filepath, "rb")
#deduplicatedaccessions=pickle.load(picklefileObj2)
#print "opened pickle deduplicatedaccessions"

#picklefileObj.close()
#picklefileObj1.close()
#picklefileObj2.close()

"""
picklefileObj3=open("%s/2nd_task/pickle.featurelessids" %filepath, "rb")
featurelessids=pickle.load(picklefileObj3)
picklefileObj4=open("%s/2nd_task/pickle.notplasmidids" %filepath, "rb")
notplasmidids=pickle.load(picklefileObj4)
combined=list(featurelessids)+list(notplasmidids)
notcdsplasmids=set(combined)
picklefileObj3.close()
picklefileObj4.close()
"""






















########################################################################################################OLD CODE (PRE-NOV 2016 RE-WRITE)

###############TESTING



"""
#TESTING DATES
#testing - accession NZ_CP0123323.1 - megablasted on ncbi and the earliest date of matching seuqence didn't correspond to my earliestdate value (which was earlier). This code confirms that the right date was selected and re-running OLDdeduplication folder code to extract date showed that CP013323.1 has been updated which explains the discrepancy.  


picklefileObj=open("%s/2nd_task/pickle.plasmid_dict" %filepath, "rb")
plasmid_dict=pickle.load(picklefileObj)

plasmid_file = open("%s/2nd_task/enterobacplasmidseqdownload.tsv" %filepath)
for line in plasmid_file:
    data=line.split('\t')
    if str(data[0])=="NZ_CP013323.1":
        print data[0], "accession"
        sequence=data[6]

print plasmid_dict[sequence]

"""

###############TESTING deduplication

"""
#find deduplicated accessions

picklefileObj=open("%s/2nd_task/pickle.plasmid_dict" %filepath, "rb")
plasmid_dict=pickle.load(picklefileObj)     


#find best unique accessions
counter1=0
counter2=0
counter3=0
counter4=0
counter5=0

deduplicatedaccessions=[]
for indx, (key, values) in enumerate(plasmid_dict.iteritems()):
    if len(values)>1:  
        features=[]
        accessions=[]
        bestaccessionindex=[]
        for value in values:
            features.append(value[1])
            accessions.append(value[0])            
        bestaccessionindex=[(a,b) for a, b in zip(accessions,features) if b == max(features)]#first select based on most features
        if len(bestaccessionindex)==1:
            if counter1==1:
                continue
            deduplicatedaccessions.append(bestaccessionindex[0][0]) ###if there is one accession with most number of features
            counter1=1
            print values, "one accession with most features", bestaccessionindex[0][0]
        else:
            refseqaccessionindex=[i for i in bestaccessionindex if '_' in i[0]]
            if len(refseqaccessionindex)==1:
                if counter2==1:
                    continue
                deduplicatedaccessions.append(refseqaccessionindex[0][0]) ###if there are multiple accessions with same no. features and one refseq accession
                counter2=1
                print values, "one refseq accession", refseqaccessionindex[0][0]
            elif len(refseqaccessionindex)==0:
                if counter3==1:
                    continue
                deduplicatedaccessions.append(bestaccessionindex[0][0]) ###if there are multiple accessions with same no. features but no refseqs
                counter3=1
                print values, "no refseqs", bestaccessionindex[0][0]

            else:
                if counter4==1:
                    continue
                deduplicatedaccessions.append(refseqaccessionindex[0][0]); print "error - duplicate refseqs"  ###if there are multiple accessions with same no. features and multiple refseqs
                counter4=1
                print values, "multiple refseqs", refseqaccessionindex[0][0], indx

    else:
        if counter5==1:
            continue
        deduplicatedaccessions.append(values[0][0])  ###if there is only one accession
        counter5=1
        print values, "one accessions", values[0][0]
        
print counter4


"""





####################OLD CODE

##############code moved to separate scripts:

#####convert nucleotide sequences to FASTA file for cd-hit-est
"""
fileObj=open("%s/2nd_task/enterobacnucleotideseqdownload.tsv" %filepath)
#fileObj2=open("%s/2nd_task/enterobacnucleotideseqexcluded_noprotein.tsv" %filepath, "w")  #previously excluded nucleotides with no associated protein for consistency with cdhit, but this is no longer necessary.
comments=[]
sequences=[]

for idx, line in enumerate(fileObj):
    print idx
    data = line.split('\t')
    comment, nucleotideseq, translation = data
    #if translation.count('-')>0:
    #    fileObj2.write('%s\t%s\t%s' %(comment,nucleotideseq,translation))   #remove non-coding DNA; n.b. there is a space after the "-" so == "-" doesn't work and don't need \n
    #else:
        #comments.append(comment)
        #sequences.append(nucleotideseq)
    comments.append(comment)
    sequences.append(nucleotideseq)
    
fileObj.close()

writeFastaSeqs(comments, sequences, "%s/2nd_task/enterobacnucleotide.fasta" %filepath)

print "finished writing fasta"
"""




#cd-hit-est code

#example: cd-hit-est -i est_human -o est_human95 -c 0.95 -n 10 -d 0 -M 16000 - T 8   #d is len description; M is memory limit; T is number of threads
"""
cmdArgs=['cd-hit-est', '-i', 'enterobacnucleotide.fasta', '-o', './cdhitest_output/output_80covthreshold/enterobac_accurate_coverage80_%s_%s.zip' %(sys.argv[1], sys.argv[2]), '-c', '%s' %sys.argv[1], '-n', '%s' %sys.argv[2], '-g', '1', '-aL', '0.80']  #-T 6 can be added for more time consuming thresholds; , '-T', '6'
subprocess.call(cmdArgs)
"""







#####convert protein sequences to FASTA file for cd-hit
"""
fileObj=open("%s/2nd_task/enterobacnucleotideseqdownload.tsv" %filepath)

comments=[]
sequences=[]

for idx, line in enumerate(fileObj):
    print idx
    data = line.split('\t')
    comment, nucleotideseq, translation = data
    if translation.count('-')>0:
        continue   #CDS with no protein annotation not included in cdhit (would be excluded automatically anyway)
    else:
        comments.append(comment)
        sequences.append(translation)
    
fileObj.close()

writeFastaSeqs(comments, sequences, "%s/2nd_task/enterobactranslation.fasta" %filepath)

print "finished writing fasta"

"""


#cdhit code
"""
cmdArgs=['cd-hit', '-i', './enterobactranslation.fasta', '-o', './../cdhit_output/output_80covthreshold/enterobac_accurate_coverage80_%s.zip' %sys.argv[1], '-c', '%s' %sys.argv[1], '-n', '5', '-g', '1', '-aL', '0.80']  #-n is always 5 at the pids I'm using
subprocess.call(cmdArgs)
"""





############OTHER OLD CODE

#getting deduplicated accessions  -correct code but no earliest date extraction

"""
deduplicatedaccessions=[]
for indx, (key, values) in enumerate(plasmid_dict.iteritems()):
    if len(values)>1:  
        features=[]
        accessions=[]
        bestaccessionindex=[]
        for value in values:
            features.append(value[1])
            accessions.append(value[0])            
        bestaccessionindex=[(a,b) for a, b in zip(accessions,features) if b == max(features)]#first select based on most features
        if len(bestaccessionindex)==1:
            deduplicatedaccessions.append(bestaccessionindex[0][0]) ###if there is one accession with most number of features
        else:
            refseqaccessionindex=[i for i in bestaccessionindex if '_' in i[0]]
            if len(refseqaccessionindex)==1:
                deduplicatedaccessions.append(refseqaccessionindex[0][0]) ###if there are multiple accessions with same no. features and one refseq accession
            elif len(refseqaccessionindex)==0:
                deduplicatedaccessions.append(bestaccessionindex[0][0]) ###if there are multiple accessions with same no. features but no refseqs
            else:
                deduplicatedaccessions.append(refseqaccessionindex[0][0]); print "error - duplicate refseqs"  ###if there are multiple accessions with same no. features and multiple refseqs
    else:
        deduplicatedaccessions.append(values[0][0])  ###if there is only one accession

print len(deduplicatedaccessions), "len deduplciated accessions"
savepickle(deduplicatedaccessions, "%s/2nd_task/pickle.deduplicatedaccessions" %filepath)
"""




#programmatic download
#first need to do esearch and then use efetch to get more details using the ids
"""
#handle = Entrez.esearch(db="nucleotide", term="plasmid AND complete sequence AND Enterobacteriaceae AND pTC2") #trial
#handle = Entrez.esearch(db="nucleotide", term="plasmid AND complete sequence AND Enterobacteriaceae") #this doesn't work: Annotation too long for DBLINK line  
handle = Entrez.esearch(db="nucleotide", term="plasmid AND complete sequence AND enterobacteria[porgn] AND biomol_genomic[PROP] AND plasmid[filter]")
recordObj = Entrez.read(handle, "genbank") #SeqRecord
handle.close()

fileObj2 = open("%s/2nd_task/enterobacdownloadTRIAL.gb" %filepath, "w")
for indx, record in enumerate(SeqIO.parse(fileObj, "genbank")):
    print indx
    SeqIO.write(record, fileObj2, "genbank")   #must use SeqIO.write not fileObj2.write('%s\t' % record.id) to maintain format


fileObj.close()
fileObj2.close()

#create genbank file
#fileObj = Entrez.efetch(db="nucleotide", id=recordObj["IdList"], rettype="gb", retmode="text")
#print recordObj["IdList"]
"""

#example 2
"""
fileObj = Entrez.esearch(db="nucleotide", term="HG941719.1")
record= Entrez.read(fileObj)
p_handle=Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta")
print p_handle.read()

"""










#cd-hit code
"""
#c is sequence identity threshold, n is wordlength, g 1 specifies accurate clustering, aL is the alignment coverage (alignement length/representative length) 
#cd-hit -i enterobactranslation.fasta -o enterobacoutput.zip -c 0.9 -n 5
#cd-hit -i enterobactranslation.fasta -o enterobacoutput_accurate_v2.zip -c 0.9 -n 5 -g 1
#cd-hit -i enterobactranslation.fasta -o enterobacoutput_accurate_coverage90.zip -c 0.9 -n 5 -g 1 -aL 0.90
#cmdArgs=['cd-hit', '-i', 'enterobactranslation.fasta', '-o', 'enterobac_accurate_coverage90_%s.zip' %sys.argv[1], '-c', '%s' %sys.argv[1], '-n', '5', '-g', '1', '-aL', '0.90']
cmdArgs=['cd-hit', '-i', 'enterobactranslation.fasta', '-o', './cdhit_output/output_80covthreshold/enterobac_accurate_coverage80_%s.zip' %sys.argv[1], '-c', '%s' %sys.argv[1], '-n', '5', '-g', '1', '-aL', '0.80']
subprocess.call(cmdArgs)
"""









"""
for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    if j<1000000:
        print j, "j"
        if rec.id not in list(deduplicated):
            record_id = rec.id
            for i, feature in enumerate(rec.features):
                if i<1:
                    plasmids.append(unlist(feature.qualifiers['plasmid']))
                else:
                    break
        else:
            continue
    else:
        break

print plasmids

fileObj.close()
"""

"""
fileObj = open("/mnt/nfs/home/alexo/2nd_task/Enterobac_3725hits_sequence.gb")#Enterobac_3725hits_sequence.gb #download.gb

for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    print j#total: 3725

print len(deduplicated), "deduplicated excluded" #362


picklefileObj=open("/mnt/nfs/home/alexo/2nd_task/pickle.plasmids", "rb")
plasmids=pickle.load(picklefileObj)
print len(plasmids), "plasmid len" #3363

"""





"""




for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    if j<10:
        print j, "j"
        if rec.id not in list(deduplicated):
            record_id = rec.id
            for i, feature in enumerate(rec.features):
                if i<1:
                    print feature.qualifiers['plasmid'], "plasmid"
                else:
                    break
        else:
            break
    else:
        break



"""







"""
for j, rec in enumerate(SeqIO.parse(fileObj, "genbank")):
    if j<10:
        print j, "j"
        if rec.id not in list(deduplicated):
            record_id = rec.id
            record_description = rec.description
            record_organism = rec.annotations['organism']
            record_length = len(rec.seq)
            record_features = len(rec.features)
            #featurecounter=0
            for i, feature in enumerate(rec.features):
                #print i
                if feature.type == "CDS":
                    #strand = feature.location.strand  #1 is positive strand, -1 is negative strand
                    #print i, "first"
                    plasmid = feature.qualifiers['plasmid']
                    #fileObj2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_organism, strand, startlocation, endlocation))
                    #fileObj3.write('%s|%s|%s\t%s\n' % (protein_id, product, record_organism, translation))
                    #fileObj3.write('%s\t%s\n' % (protein_id, translation)) #this doesn't make extracting protein id easier
                #if feature.type == "CDS" and i>0:
                    #print i, "after first"
                    #strand = feature.location.strand  #1 is positive strand, -1 is negative strand
                    #fileObj2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_organism, strand, startlocation, endlocation, gene, product, protein_id, db_xref, locus_tag, note))
                    #fileObj3.write('%s|%s|%s\t%s\n' % (protein_id, product, record_organism, translation))
                    #fileObj3.write('%s\t%s\n' % (protein_id, translation)) #this doesn't make extracting protein id easier
            print plasmid, j, i
        else:
            excludedids.append(rec.id)
    else:
        break

fileObj.close()
fileObj2.close()
#fileObj3.close()
#fileObj4.close()

"""












#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC

#fileObj2 = open("trialoutput.fasta", "w")
#for i, proteinSeq in enumerate(sequences):
#    seqObj = Seq(proteinSeq, IUPAC.protein)
#    proteinObj = SeqRecord(seqObj, id="TEST")
#    SeqIO.write([proteinObj,], fileObj2, 'fasta')





#enumerate and if reaches end without writing to file record as 'no features'
#could filter at this stage - if taxa != enteobac...

#in outer nest, add rec.name, rec.description etc. (understand)



#code before I wrote code to make notplasmid and featureless pickle objects:
"""
featurelessids=[]

fileObj = open("/mnt/nfs/home/alexo/2nd_task/Enterobac_3725hits_sequence.gb")#Enterobac_3725hits_sequence.gb #download.gb
fileObj2 = open("/mnt/nfs/home/alexo/2nd_task/enterobacfeaturesdownload.tsv", "w") #enterobacdownload.tsv
fileObj3 = open("/mnt/nfs/home/alexo/2nd_task/enterobactranslationdownload.tsv", "w") #enterobacdownload.fasta
fileObj4 = open("/mnt/nfs/home/alexo/2nd_task/enterobacrecorddownload.tsv", "w")

for rec in SeqIO.parse(fileObj, "genbank"):
    if rec.features:
        record_id = rec.id
        #name = rec.name #can tell from rec.id
        record_description = rec.description
        record_organism = rec.annotations['organism']
        record_length = len(rec.seq)
        record_features = len(rec.features)
        featurecounter=0
        for feature in rec.features:
            if feature.type == "CDS":
                featurecounter=featurecounter+1
                strand = feature.location.strand  #1 is positive strand, -1 is negative strand
                startlocation = feature.location.start
                endlocation = feature.location.end
                try:
                    db_xref=unlist(feature.qualifiers["db_xref"])
                except KeyError:
                    db_xref="-"
                try:
                    note=unlist(feature.qualifiers["note"])
                except KeyError:
                    note="-"
                try:
                    gene=unlist(feature.qualifiers["gene"])
                except KeyError:
                    gene="-"
                try:
                    locus_tag=unlist(feature.qualifiers["locus_tag"])
                except KeyError:
                    locus_tag="-"
                try:
                    product=unlist(feature.qualifiers["product"])
                except KeyError:
                    product="-"
                try:
                    protein_id=unlist(feature.qualifiers["protein_id"])
                except KeyError:
                    protein_id="-"
                try:
                    translation=unlist(feature.qualifiers["translation"])
                except KeyError:
                    translation="-"
                fileObj2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (record_id, record_organism, strand, startlocation, endlocation, gene, product, protein_id, db_xref, locus_tag, note))
                fileObj3.write('%s|%s|%s\t%s\n' % (protein_id, product, record_organism, translation))
        record_features_CDS=featurecounter
        fileObj4.write('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (record_id, record_description, record_organism, record_length, record_features, record_features_CDS)) 
    else:
        if not rec.features:
            featurelessids.append(rec.id)

print len(featurelessids)


fileObj.close()
fileObj2.close()
fileObj3.close()
fileObj4.close()

"""





"""
fileObj = Entrez.efetch(db="nucleotide", id=recordObj["IdList"], rettype="gb", retmode="text")
fileObj2 = open('/mnt/nfs/home/alexo/2nd_task/download.gb', "w")

for record in SeqIO.parse(fileObj, "genbank"):  #when using SeqIO.read() ...the handle can't be assigned multiple records (ValueError: more than one record found in handle)
    fileObj2.write('%s\t' % record.id)   #write a specific SeqRecord object to file
    print record.features
fileObj.close()
fileObj2.close()

"""

"""
fileObj2 = open("/mnt/nfs/home/alexo/2nd_task/download.gb", "w")
SeqIO.write(records, fileObj2, "genbank")
fileObj.close()
fileObj2.close()
"""


"""
#GFF parsing


import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
 
in_file = "/mnt/nfs/home/alexo/2nd_task/download.gb"
out_file = "/mnt/nfs/home/alexo/2nd_task/download.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")
 
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
 
in_handle.close()
out_handle.close()
"""


"""
fileObj = open("/mnt/nfs/home/alexo/2nd_task/download.gb")
featureslist=[]
for record in SeqIO.parse(fileObj, "genbank"):
    featureslist.append(record.features.type)
    
fileObj.close()


print featureslist[0]



#use genbank file
#old code#######

#create FASTA files
fileObj = open("/mnt/nfs/home/alexo/2nd_task/download.gb")
fileObj2 = open("/mnt/nfs/home/alexo/2nd_task/download.fasta", "w")
for record in SeqIO.parse(fileObj, "genbank"):
    SeqIO.write(record, fileObj2, "fasta")
    
fileObj.close()
fileObj2.close()

#look at fasta files:

fileObj = open("/mnt/nfs/home/alexo/2nd_task/download.fasta")
for trial in SeqIO.parse(fileObj, "fasta"):
    print trial.id
    print repr(trial.seq)
"""


"""
fileObj = open("/mnt/nfs/home/alexo/2nd_task/download.gb")
for rec in SeqIO.parse(fileObj, "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
                print feature.location
                print feature.qualifiers["protein_id"]
                print feature.location.extract(rec).seq

#############



qualifiers=["db_xref", "note"]
for i in qualifiers:
list[i]=[]
                for i in qualifiers:
                list[i]=[]
                    try:
                        list[i]=feature.qualifiers[qualifier]
                    except KeyError:
                        list[i]="-"




                print nestedlist[i][0]
                print nestedlist[i][1] 









filtering cd-hit 
# parse through the cluster file and store the cluster name + sequences in the dictionary
from itertools import *
cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
for cluster in cluster_groups:
    name = cluster.next().strip()
    seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.next()]
    cluster_dic[name] = seqs
print cluster_dic

"""

#descriptive statistics
"""
fileObj=open("/mnt/nfs/home/alexo/2nd_task/enterobacrecorddownload.tsv")
accessions=[]
taxas=[]
plasmid_sizes=[]
features=[]
CDSfeatures=[]

for idx, line in enumerate(fileObj):
    data = line.split('\t')
    accession, taxa, plasmid_size, feature, CDSfeature = data
    accessions.append(accession)
    taxas.append(taxa)
    plasmid_sizes.append(plasmid_size)
    features.append(feature)
    CDSfeatures.append(CDSfeature)

fileObj.close()
 

import matplotlib.pyplot as plt
import numpy as np

import plotly.plotly as py
# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api


plt.hist(plasmid_sizes)
plt.title("Plasmid Sizes Histogram")
plt.xlabel("Length of plasmid genome (bp)")
plt.ylabel("Frequency")

fig = plt.gcf()

plot_url = py.plot_mpl(fig, filename='mpl-basic-histogram')

"""


#count organism occurences




"""
data = '''\
ashwin programmer india
amith programmer india'''

c = Counter()
for line in data.splitlines():
    c.update(line.split())
print(c)
"""


#removing tails of histogram


"""
fileObj = open("/mnt/nfs/home/alexo/2nd_task/plasmidsequence.gb")
trial=SeqIO.parse(fileObj, "genbank")
for i, feature in enumerate(trial.features):
    while i<10:
        print feature.qualifiers['plasmid']
    #print rec.features

"""


"""
print "next plasid"

fileObj2 = open("/mnt/nfs/home/alexo/2nd_task/chromosomalsequence.gb")
trial2=SeqIO.parse(fileObj2, "genbank")
for rec in trial2:
    for i, feature in enumerate(rec.features):
        while i<10:
            print feature.qualifiers['plasmid']

"""
    #print rec.features






#can't determine plasmid source at seqrecord object level:
"""
fileObj = open("/mnt/nfs/home/alexo/2nd_task/plasmidsequence.gb")
trial=SeqIO.parse(fileObj, "genbank")
print #dir(trial)
for rec in trial:
    print rec.annotations
    #print rec.features
"""







#fileObj.close()
#fileObj2.close()
#fileObj3.close()
#fileObj4.close()
