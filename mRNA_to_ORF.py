"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""

#Import des modules
import os
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
for file_name in fasta_list:
    gene_name = (file_name.split('_')[0])
    dico[gene_name+'_mRNA'] = []
    dico[gene_name+'_mRNA'].append(SeqIO.to_dict(SeqIO.parse(path+file_name, 'fasta')))
    for id in dico[gene_name+'_mRNA']:
        for value in id.values():
            mydna = value.seq
            print(mydna.find('ATG'))