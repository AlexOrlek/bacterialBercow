import sys
from pythonmods import runblastn, runsubprocess, mlstfilter

#do blastn of accessions_filtered.fa file against rmlst database

rmlstdbpath=sys.argv[1]
threads=sys.argv[2]
outdir=sys.argv[3]
sequenceorigin=sys.argv[4]
sourcedir=sys.argv[5]
if sequenceorigin=='ncbi':
    query='%s/accessions_filtered_deduplicated.fa'%outdir
else:
    query=sys.argv[6]

cmdArgs='find %s'%rmlstdbpath+' '+'-maxdepth 1 -mindepth 1 -name "*.n??" -printf "%f\n" | rev | cut -d"." -f2- | rev | sort | uniq >'+' '+'%s/rmlstloci.txt'%rmlstdbpath
runsubprocess([cmdArgs],shell=True)
    
with open('%s/rmlstloci.txt'%rmlstdbpath) as f:
    for line in f:
        locus=line.strip()
        database='%s/%s'%(rmlstdbpath,locus)
        blastoutput='%s/rmlst/BLASTtable_%s.tsv' %(outdir,locus)
        runblastn(query, database, blastoutput, num_threads='%s'%threads, max_hsps='1', perc_identity='95', evalue='1e-10',word_size='28') #since running per-locus blasts, can be stringent with max_hsps***; also using stingent pid, e-value and word size
        print('runblastn finshed for locus %s'%locus)
#***"Maximum number of HSPs (alignments) to keep for any single query-subject pair. The HSPs shown will be the best as judged by expect value. This number should be an integer that is one or greater. If this option is not set, BLAST shows all HSPs meeting the expect value criteria. Setting it to one will show only the best HSP for every query-subject pair"

#combining per-locus outputs prior to filtering
args=['find %s/rmlst/ -maxdepth 1 -mindepth 1 -type f -name "BLASTtable_*.tsv" ! -name "BLASTtable_combined*.tsv" -exec cat {} \; > %s/rmlst/BLASTtable_combined.tsv'%(outdir,outdir)]
runsubprocess(args,shell=True)

#filtering
blastoutput='%s/rmlst/BLASTtable_combined.tsv'%outdir
finalfile='%s/rmlst/BLASTtablebesthits.tsv'%outdir
sortedfile='%s/rmlst/BLASTtablesorted.tsv'%outdir
mlstfilter(blastoutput, finalfile, sortedfile, sourcedir, pidthresh=95, coveragethresh=0.95, formatcontigcol=False) #95% coverage/pid are the thresholds used by Larsen et al. 2014 "Benchmarking methods for genomic taxonomy"; and a high coverage threhsold is appropriate for complete (or near complete) genomes
print('mlstfilter finished')




#OLD CODE - !PROBLEM: there are 1239737 sequences in the combined database - would need to set max_target_seqs to infeasibly large cutoff; best to search databases individually (each database comprises ~ 30000 alleles)
#idtable='%s/idtable.tsv'%sys.argv[3] #don't need id tables


# database='%s/allallelesdb'%sys.argv[3]
# idtable='%s/idtable.tsv'%sys.argv[3]

# query='%s/output/%s_assemblies.fasta' %(sys.argv[1],sys.argv[2])
# blastoutput='%s/rmlst/output/%s_BLASTtable.tsv' %(sys.argv[1],sys.argv[2])
# runblastn(query, database, blastoutput, num_threads='%s'%sys.argv[4]) #running with default evalue and word size   
# print('runblastn finshed')
# finalfile='%s/rmlst/output/%s_BLASTtablebesthits.tsv'%(sys.argv[1],sys.argv[2])
# sortedfile='%s/rmlst/output/%s_BLASTtablesorted.tsv'%(sys.argv[1],sys.argv[2])
# blastfilter(blastoutput, finalfile, sortedfile, idtable, pidthresh=95, coveragethresh=0.95) #this is the threshold use in larsen 2014, and a high coverage threhsold is appropriate for complete (or near complete) genomes !? - should also apply to mlst , resfinder etc.
# print('blastfilter finished')
