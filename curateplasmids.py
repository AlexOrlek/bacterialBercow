#!/usr/bin/env python
import argparse, os, sys, subprocess, signal
sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess

def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def positiveint(x):
    x = int(x)
    if x < 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" %x)
    return x

def batchsizeint(x):
    x = int(x)
    if x < 2:
         raise argparse.ArgumentTypeError("%s is too small; batch size must be greater than 1" %x)
    if x > 500:
         raise argparse.ArgumentTypeError("%s is too large; batch size must not exceed 500" %x)
    return x

parser = argparse.ArgumentParser(description='Run pipeline scripts')
parser.add_argument('-q','--taxonomyquery', help='Taxonomy search query term to be supplied to the edirect eseach -query argument (default: bacteria[porgn:__txid2])', default="bacteria[porgn:__txid2])", type=str)
parser.add_argument('-d','--datequery', help='Date search query term to be supplied to the edirect eseach -query argument (e.g. "2017/01/01"[PDAT] : "3000"[PDAT] would retrieve records since 2017) (not required)', required=False, type=str)
parser.add_argument('-b','--batchsize', help='Number of accession nucleotide records to retrieve per edirect query (default: 500; min: 2; max: 500)', default=500, type=batchsizeint)
parser.add_argument('-e','--emailaddress', help="User's email address which will be provided as an argument to edirect econtact -email (required if retrieving data from NCBI i.e. required except if --fasta file is provided)", required=False, type=str)
parser.add_argument('-s','--dbsource', help='Database source; refseq or refseq_genbank (default: refseq_genbank)', default="refseq_genbank", choices=["refseq","refseq_genbank"],type=str)
parser.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=positiveint)
parser.add_argument('-o','--out', help='Output directory (required)', required=True, type=str)
parser.add_argument('--rmlstdbpath', help='Path to the directory used to store the rmlst database (default: databases/rmlstalleles)',required=False)
parser.add_argument('--enterobacdbpath', help='Path to the "enterobacteriaceae" plasmidfinder BLAST database (default: databases/plasmidfinder/enterobacteriaceae/enterobacteriaceaedb)',required=False)
parser.add_argument('--gramposdbpath', help='Path to the "gram_positive" plasmidfinder BLAST database (default: databases/plasmidfinder/gram_positive/gram_positivedb)',required=False)
parser.add_argument('--accessions', help='A text file containing NCBI plasmid accessions in the first column; if provided, these accessions will be retrieved, rather than retrieving plasmid accessions using a query term (default: retrieve accessions using a query term)',required=False)
parser.add_argument('--sequences', help='A fasta file containing uncharacterised bacterial contig nucleotide sequences; if provided, these contigs will be typed using rmlst and replicon loci to determine whether they are likely to be plasmids or chromosomal (default: retrieve sequences from NCBI)',required=False)
parser.add_argument('--retrieveaccessionsonly', action='store_true',help='If flag is provided, stop after retrieving and filtering accessions (default: do not stop)',required=False)
parser.add_argument('--retrievesequencesonly', action='store_true',help='If flag is provided, stop after retrieving sequences from filtered accessions (default: do not stop)',required=False)

args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)

cmdArgs=['mkdir -p %s'%outputpath]
runsubprocess(cmdArgs,shell=True)

###retrieve accessions and sequences from NCBI
if args.sequences==None:
    if args.accessions==None:
        if args.datequery==None:
            datepresent="absent"
        else:
            datepresent=="present"
        runsubprocess(['bash','%s/downloadaccessions.sh'%sourcedir,datepresent,str(args.taxonomyquery),str(args.datequery),str(args.dbsource),outputpath])
        runsubprocess(['python','%s/filteraccessions.py'%sourcedir,outputpath])
    else:
        runsubprocess(['bash','%s/downloaduseraccessions.sh'%sourcedir,str(args.accessions),outputpath])
        runsubprocess(['python','%s/filteraccessions.py'%sourcedir,outputpath])
    ###retrieve sequences if args.retrieveaccessionsonly is false    
    if args.retrieveaccessionsonly==True:
        sys.exit()
    else:
        runsubprocess(['bash','%s/downloadsequences.sh'%sourcedir,str(args.batchsize),str(args.emailaddress),outputpath])
        runsubprocess(['python','%s/deduplicateseqs.py'%sourcedir,outputpath])
        


if args.retrieveaccessionsonly==True:
    sys.exit()

if args.retrievesequencesonly==True:
    sys.exit()


###characterise sequences to identify plasmids
cmdArgs=['mkdir -p %s/plasmidfinder'%outputpath]
runsubprocess(cmdArgs,shell=True)
cmdArgs=['mkdir -p %s/rmlst'%outputpath]
runsubprocess(cmdArgs,shell=True)

if args.enterobacdbpath==None:
    enterobacteriaceaedbpath='%s/databases/plasmidfinder/enterobacteriaceae/enterobacteriaceaedb'%sourcedir
else:
    enterobacteriaceaedbpath=str(args.enterobacdbpath)
    
if args.gramposdbpath==None:
    gram_positivedbpath='%s/databases/plasmidfinder/gram_positive/gram_positivedb'%sourcedir
else:
    gram_positivedbpath=str(args.gramposdbpath)
      

if args.rmlstdbpath==None:
    rmlstdbpath='%s/databases/rmlstalleles'%sourcedir
else:
    rmlstdbpath=str(args.rmlstdbpath)


if args.sequences==None:
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'enterobacteriaceae',enterobacteriaceaedbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'gram_positive',gram_positivedbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/rmlst.py'%sourcedir,rmlstdbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/finalfilter.py'%sourcedir, outputpath, 'enterobacteriaceae', 'gram_positive'])
else:
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'enterobacteriaceae',enterobacteriaceaedbpath,str(args.threads),outputpath,'user',sourcedir,str(args.sequences)])
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'gram_positive',gram_positivedbpath,str(args.threads),outputpath,'user',sourcedir,str(args.sequences)])
    runsubprocess(['python', '%s/rmlst.py'%sourcedir,rmlstdbpath,str(args.threads),outputpath,'user',sourcedir,str(args.sequences)])
    runsubprocess(['python', '%s/finalfilter.py'%sourcedir, outputpath, 'enterobacteriaceae', 'gram_positive'])


