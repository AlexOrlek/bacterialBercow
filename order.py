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

parser = argparse.ArgumentParser(description='bacterialBercow: bringing order to bacterial sequences',add_help=False)
#Help options
help_group = parser.add_argument_group('Help')
help_group.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
#General options                                                                                                                                                         
general_group = parser.add_argument_group('General options')
general_group.add_argument('-o','--out', help='Output directory (required)', required=True, type=str)
#NCBI query and retrieval options
ncbi_group = parser.add_argument_group('NCBI query and retrieval options')
ncbi_group.add_argument('-e','--emailaddress', help="User's email address which will be provided as an argument to edirect econtact -email (required if retrieving data from NCBI)", required=False, type=str)
ncbi_group.add_argument('--taxonomyquery', help='Taxonomy search query term to be supplied to the edirect eseach -query argument (default: bacteria[porgn:__txid2])', default="bacteria[porgn:__txid2])", type=str)
ncbi_group.add_argument('--datequery', help='Date search query term to be supplied to the edirect eseach -query argument (e.g. "2017/01/01"[PDAT] : "3000"[PDAT] would retrieve records since 2017) (not required)', required=False, type=str)
ncbi_group.add_argument('-s','--dbsource', help='Database source; refseq or refseq_genbank (default: refseq_genbank)', default="refseq_genbank", choices=["refseq","refseq_genbank"],type=str)
ncbi_group.add_argument('--deduplicationmethod', help='Specify how identical sequences should be deduplicated; either "all" duplicates are removed; otherwise, duplicates are removed if they share biosample accession id + "submitter" metadata; or "bioproject" accession id; or "both" submitter metadata and bioproject accession id (default: "both")', default="both", choices=["both","submitter","bioproject","all"],type=str)
ncbi_group.add_argument('-b','--batchsize', help='Number of accession nucleotide records to retrieve per edirect query (default: 200; min: 2; max: 500)', default=200, type=batchsizeint)


#Replicon and rMLST typing options
typing_group = parser.add_argument_group('Replicon and rMLST typing options')
typing_group.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=positiveint)
typing_group.add_argument('--typing', help='Specifies what sequence typing to perform (only applicable if in-house sequences are provided using --inhousesequences flag); either "replicon", "rmlst" typing or "both" (default: both)',default="both",choices=["both","replicon","rmlst"],required=False)


#Pipeline step customisation (specifying starting and stopping points)
steps_group = parser.add_argument_group('Customising pipeline steps (specifying starting / stopping points)')
steps_group.add_argument('--accessions', help='A text file containing NCBI plasmid accessions in the first column; if provided, these accessions will be retrieved, rather than retrieving plasmid accessions using a query term (default: retrieve accessions using a query term)',required=False)
steps_group.add_argument('--retrieveaccessionsonly', action='store_true',help='If flag is provided, stop after retrieving and filtering accessions (default: do not stop)',required=False)
steps_group.add_argument('--retrievesequencesonly', action='store_true',help='If flag is provided, stop after retrieving deduplicated sequences from filtered accessions (default: do not stop)',required=False)
steps_group.add_argument('--restartwithsequences', action='store_true',help='If flag is provided, re-start the pipeline using sequences retrieved from NCBI',required=False)
steps_group.add_argument('--inhousesequences', help='A fasta file containing uncharacterised bacterial contig nucleotide sequences; if provided, these contigs will be typed using rmlst and replicon loci to determine whether they are likely to be plasmids or chromosomal (default: retrieve sequences from NCBI)',required=False)

args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)

cmdArgs=['mkdir -p %s'%outputpath]
runsubprocess(cmdArgs,shell=True)

###retrieve accessions and sequences from NCBI
if args.inhousesequences==None and args.restartwithsequences==False:
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
        runsubprocess(['python','%s/deduplicateseqs.py'%sourcedir,str(args.deduplicationmethod),outputpath])
        


if args.retrieveaccessionsonly==True:
    sys.exit()

if args.retrievesequencesonly==True:
    sys.exit()


###characterise sequences to identify plasmids
cmdArgs=['mkdir -p %s/plasmidfinder'%outputpath]
runsubprocess(cmdArgs,shell=True)
cmdArgs=['mkdir -p %s/rmlst'%outputpath]
runsubprocess(cmdArgs,shell=True)

enterobacteriaceaedbpath='%s/databases/plasmidfinder_db/blastdbs/enterobacteriaceaedb'%sourcedir
gram_positivedbpath='%s/databases/plasmidfinder_db/blastdbs/gram_positivedb'%sourcedir
rmlstdbpath='%s/databases/rmlstalleles/blastdbs'%sourcedir
rmlstprofilepath='%s/databases/rmlstalleles'%sourcedir

if args.inhousesequences==None:
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'enterobacteriaceae',enterobacteriaceaedbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'gram_positive',gram_positivedbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/rmlst.py'%sourcedir,rmlstdbpath,str(args.threads),outputpath,'ncbi',sourcedir])
    runsubprocess(['python', '%s/finalfilter.py'%sourcedir, rmlstprofilepath,outputpath, 'ncbi','enterobacteriaceae', 'gram_positive'])
else:
    cmdArgs=["cat %s | bioawk -c fastx '{print $name,length($seq)}' > %s/seqlengths.tsv"%(str(args.inhousesequences),outputpath)]
    runsubprocess(cmdArgs,shell=True)
    if args.typing=='replicon' or args.typing=='both':
        runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'enterobacteriaceae',enterobacteriaceaedbpath,str(args.threads),outputpath,'user',sourcedir,str(args.inhousesequences)])
        runsubprocess(['python', '%s/plasmidfinder.py'%sourcedir,'gram_positive',gram_positivedbpath,str(args.threads),outputpath,'user',sourcedir,str(args.inhousesequences)])
    if args.typing=='rmlst' or args.typing=='both':
        runsubprocess(['python', '%s/rmlst.py'%sourcedir,rmlstdbpath,str(args.threads),outputpath,'user',sourcedir,str(args.inhousesequences)])
    runsubprocess(['python', '%s/finalfilter.py'%sourcedir, rmlstprofilepath,outputpath,'user',str(args.typing),'enterobacteriaceae', 'gram_positive'])





###OLD CODE


#typing_group.add_argument('--enterobacdbpath', help='Path to the "enterobacteriaceae" plasmidfinder BLAST database (default: databases/plasmidfinder/enterobacteriaceae/enterobacteriaceaedb)',required=False)
#typing_group.add_argument('--gramposdbpath', help='Path to the "gram_positive" plasmidfinder BLAST database (default: databases/plasmidfinder/gram_positive/gram_positivedb)',required=False)
#typing_group.add_argument('--rmlstdbpath', help='Path to the directory used to store the rmlst blast database files (default: databases/rmlstalleles/blastdbs)',required=False)
#typing_group.add_argument('--rmlstprofilepath', help='Path to the directory used to store the rmlst profile file (default: databases/rmlstalleles)',required=False)


# if args.enterobacdbpath==None:
#     enterobacteriaceaedbpath='%s/databases/plasmidfinder/enterobacteriaceae/enterobacteriaceaedb'%sourcedir
# else:
#     enterobacteriaceaedbpath=str(args.enterobacdbpath)
    
# if args.gramposdbpath==None:
#     gram_positivedbpath='%s/databases/plasmidfinder/gram_positive/gram_positivedb'%sourcedir
# else:
#     gram_positivedbpath=str(args.gramposdbpath)
      

# if args.rmlstdbpath==None:
#     rmlstdbpath='%s/databases/rmlstalleles/blastdbs'%sourcedir
# else:
#     rmlstdbpath=str(args.rmlstdbpath)

# if args.rmlstprofilepath==None:
#     rmlstprofilepath='%s/databases/rmlstalleles'%sourcedir
# else:
#     rmlstprofilepath=str(args.rmlstprofilepath)

