#Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re

#Data
deg_database = pd.read_csv('essential_genes/deg.csv')
df = pd.read_excel('data/ifo0880/elife-32110-supp1-v2.xlsx', sheet_name = 'Essential Genes')
eg_list = list(map(int, list(df.loc[df['Essential'] == 'Yes', 'Protein ID']))) #Protein IDs of essential genes

#Functions
def read_eg_fasta(name, check_list):
    fasta_seqs = SeqIO.parse(open(name),'fasta')
    data = []
    for fasta in fasta_seqs:
        if int(fasta.description.split('|')[2]) in check_list:
            if str(fasta.seq).strip().upper()[-1] == '*':
                data.append([fasta.description.split('|')[2], str(fasta.seq).strip().upper()[:-1], 'Rhodosporidium toruloides IFO0880'])
            else:
                data.append([fasta.description.split('|')[2], str(fasta.seq).strip().upper(), 'Rhodosporidium toruloides IFO0880'])
    return data

#Essential genes in IFO0880
proteins = read_eg_fasta('data/ifo0880/s888LsZ4f6T.faa', eg_list)
print(len(proteins))

#adding to new deg.csv
deg_database = pd.concat([deg_database, pd.DataFrame(proteins, columns=deg_database.columns)], ignore_index=True)
pd.DataFrame(deg_database).to_csv('essential_genes/deg_upd.csv', index = False)





