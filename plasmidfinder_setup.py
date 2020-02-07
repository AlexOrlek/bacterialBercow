#!/usr/bin/env python
import os, datetime, subprocess
from Bio import SeqIO
from pythonmods import runsubprocess

sourcedir=os.path.dirname(os.path.abspath(__file__))

output_folder='./databases/plasmidfinder_db'
cmdArgs=['mkdir -p %s'%output_folder]
runsubprocess(cmdArgs,shell=True)

cmdArgs=['git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git ./databases/plasmidfinder_db']
runsubprocess(cmdArgs,shell=True)

gramposfastas=[]
for file in os.listdir('./databases/plasmidfinder_db'):
    if file.endswith('.fsa') and file!='enterobacteriaceae.fsa':
        gramposfastas.append(os.path.join('./databases/plasmidfinder_db', file))

#combine gram-positive replicons into single gram-positive database
f2=open('./databases/plasmidfinder_db/gram_positive.fsa','w')
for filename in gramposfastas:
    with open(filename) as f:
        for indx,seq_record in enumerate(SeqIO.parse(f,'fasta')):
            fastaheader=str(seq_record.id)
            newfastaheader='%s|%s'%(filename,fastaheader)
            seq_record.id=newfastaheader
            seq_record.description=''
            SeqIO.write(seq_record,f2,'fasta')

#make blast databases
cmdArgs=['bash', '%s/makeblastdbs.sh'%sourcedir, '%s'%output_folder, 'enterobacteriaceae']
subprocess.call(cmdArgs)
cmdArgs=['bash', '%s/makeblastdbs.sh'%sourcedir, '%s'%output_folder, 'gram_positive']
subprocess.call(cmdArgs)


current_year = datetime.datetime.utcnow().year
current_month = datetime.datetime.utcnow().month
current_day = datetime.datetime.utcnow().day
with open(os.path.join(output_folder, 'download_date.txt'), 'w') as f:
    f.write('{}-{}-{}'.format(current_year, current_month, current_day))
            
print('Finished downloading PlasmidFinder database')
