import sys, os
from pythonmods import runblastn,blastfilter

database=sys.argv[1] #either enterobacteriaceae or gram_positive
dbpath=sys.argv[2] #plasmidfinder database path, either enterobacteriaceae or gram_positive
threads=sys.argv[3]
outdir=sys.argv[4]
sequenceorigin=sys.argv[5]
sourcedir=sys.argv[6]
if sequenceorigin=='ncbi':
    query='%s/accessions_filtered_deduplicated.fa'%outdir
else:
    query=sys.argv[7]

#do blastn of accessions_filtered.fa file against plasmidfinder database

#run database search
blastoutput='%s/plasmidfinder/BLASTtable_%s.tsv' %(outdir,database)
runblastn(query, dbpath, blastoutput, num_threads='%s'%threads) #running with default evalue and word size   
print('runblastn finshed; database: %s'%database)
finalfile='%s/plasmidfinder/BLASTtablebesthits_%s.tsv'%(outdir,database)
sortedfile='%s/plasmidfinder/BLASTtablesorted_%s.tsv'%(outdir,database)
blastfilter(blastoutput, finalfile, sortedfile, sourcedir, pidthresh=80, coveragethresh=0.60)
print('blastfilter finished; database: %s'%database)

