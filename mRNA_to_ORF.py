"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""

#Import des modules
import os
import pandas
from Bio import SeqIO


#Fonctions
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f



## MAIN ##
fasta_list = listdir_nohidden('./mRNA_fasta')
path = './mRNA_fasta/'
dico={}
mRNA_name = []
mRNA_sequence = []
start = []
stop = []
gene = []
longueur_mRNA = []
orf = []

for file_name in fasta_list:
    gene_name = (file_name.split('_')[0])
    dico[gene_name+'_mRNA'] = []
    dico[gene_name+'_mRNA'].append(SeqIO.to_dict(SeqIO.parse(path+file_name, 'fasta')))
    
    for id in dico[gene_name+'_mRNA']:
        for key in id.keys():
            mRNA_name.append(key)
        for value in id.values():
            stop_list = []
            orf_list = []
            mydna = value.seq
            mRNA_sequence.append(str(value.seq))
            ATG_indice = mydna.find('ATG')
            start.append(ATG_indice)
            TGA_indice = mydna.find('TGA')
            TAG_indice = mydna.find('TAG')
            TAA_indice = mydna.find('TAA')
            stop_list.append(TGA_indice)
            stop_list.append(TAG_indice)
            stop_list.append(TAA_indice)   
            stop.append(stop_list)
            gene.append(gene_name)
            longueur_mRNA.append(len(str(value.seq)))
            orf_TGA = value.seq[ATG_indice:TGA_indice]
            orf_TAG = value.seq[ATG_indice:TAG_indice]
            orf_TAA = value.seq[ATG_indice:TAA_indice]
            orf_list.append(str(orf_TGA))
            orf_list.append(str(orf_TAG))
            orf_list.append(str(orf_TAA))
            orf.append(orf_list)


#Creation dataframe
df = pandas.DataFrame(columns = ['gene', 'mRNA_name','mRNA_sequence','start','stop','longueur_mRNA','ORF_sequences'])
df['gene'] = pandas.Series(gene)
df['mRNA_name'] = pandas.Series(mRNA_name)
df['mRNA_sequence'] = pandas.Series(mRNA_sequence)
df['start'] = pandas.Series(start)
df['stop'] = pandas.Series(stop)
df['longueur_mRNA'] = pandas.Series(longueur_mRNA)
df['ORF_sequences'] = pandas.Series(orf)

#Export vers csv
pandas.DataFrame.to_csv(df, 'ORF.csv')