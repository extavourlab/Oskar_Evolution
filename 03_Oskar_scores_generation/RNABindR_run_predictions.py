'''
Author: Savandara Besse
- Creation: 07-27-2017
- Last modification: 09-04-2019

Description :
Python script able to create automatic requests for RNA-binding sites predictions
using MSA fasta file as input (HERE: ./FASTA/OSKAR_sequences.fasta)
'''

import progressbar, requests, time, sys
from multiprocessing import Process
from Bio import SeqIO
import pandas as pd


fastaFile = sys.argv[1]
records = list(SeqIO.parse(fastaFile, 'fasta'))
ID = [record.id for record in records]
seq = [str(record.seq) for record in records]
description = [record.description for record in records]
description = [descr.split('|') for descr in description]
description = [descr[len(descr)-1] for descr in description]
seqInfos = [ID, description, seq]
df = pd.DataFrame(data=seqInfos).T
df.columns = ['ID', 'Description','Sequence']

processes = []
for index in df.index:
    payload = {'email': 'savandara.besse@gmail.com', 'JobTitle': 'JOB_{}_{}'.format(index,df['Description'][index]), 'QuerySeq':'>{}_{}\n{}'.format(df['ID'][index],df['Description'][index],df['Sequence'][index])}
    print(payload)
    p = Process(target=requests.post, args=('http://ailab1.ist.psu.edu/RNABindR/cgi-bin/predict.cgi', payload))
    processes.append(p)
for p in processes :
    p.start()
    time.sleep(180)
    try:
        for p in processes :
            p.join()
    except:
        [p.terminate() for p in processes]
