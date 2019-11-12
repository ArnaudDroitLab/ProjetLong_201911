"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""

#Import des modules
import os
import pandas
from Bio import SeqIO

def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f



## MAIN ##
fasta_list = listdir_nohidden('./mRNA_fasta')

#récupération séquences sous forme de dictionnaire - clefs = ids -
path = './mRNA_fasta/'
dico={}
mRNA_sequence = []
start = []
stop = []
gene = []

for file_name in fasta_list:
    gene_name = (file_name.split('_')[0])
    dico[gene_name+'_mRNA'] = []
    dico[gene_name+'_mRNA'].append(SeqIO.to_dict(SeqIO.parse(path+file_name, 'fasta')))
    
    for id in dico[gene_name+'_mRNA']:
        for value in id.values():
            stop_list = []
            mydna = value.seq
            mRNA_sequence.append(str(value.seq))
            start.append(mydna.find('ATG'))
            stop_list.append(mydna.find('TGA'))
            stop_list.append(mydna.find('TAG'))
            stop_list.append(mydna.find('TAA'))   
            stop.append(stop_list)
            gene.append(gene_name)


df = pandas.DataFrame(columns = ['gene', 'mRNA_sequence','start','stop','ORF_sequence'])
df['gene'] = pandas.Series(gene)
df['mRNA_sequence'] = pandas.Series(mRNA_sequence)
df['start'] = pandas.Series(start)
df['stop'] = pandas.Series(stop)

print(df)
