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
gene = []
longueur_mRNA = []
nombre_ORF = []
orf = []

for file_name in fasta_list:
    gene_name = (file_name.split('_')[0])
    dico[gene_name+'_mRNA'] = []
    dico[gene_name+'_mRNA'].append(SeqIO.to_dict(SeqIO.parse(path+file_name, 'fasta')))
    
    for id in dico[gene_name+'_mRNA']:
        for key in id.keys():
            mRNA_name.append(key)
        
        for value in id.values():
            orf_list = []
            mRNA_sequence.append(str(value.seq))
            longueur_mRNA.append(len(str(value.seq)))
            gene.append(gene_name)
            
            for cadre in range(3):
                
                #Sens
                sens = value.seq[cadre:]
                ATG_indice_sens = sens.find('ATG')
                TGA_indice_sens = sens.find('TGA')
                TAG_indice_sens = sens.find('TAG')
                TAA_indice_sens = sens.find('TAA')
                orf_TGA_sens = sens[ATG_indice_sens:TGA_indice_sens]
                orf_TAG_sens = sens[ATG_indice_sens:TAG_indice_sens]
                orf_TAA_sens = sens[ATG_indice_sens:TAA_indice_sens]
                if len(str(orf_TGA_sens)) >= 30 and len(str(orf_TGA_sens))%3==0:
                    orf_list.append(str(orf_TGA_sens))
                if len(str(orf_TAG_sens)) >= 30 and len(str(orf_TAG_sens))%3==0:
                    orf_list.append(str(orf_TAG_sens))
                if len(str(orf_TAA_sens)) >= 30 and len(str(orf_TAA_sens))%3==0:
                    orf_list.append(str(orf_TAA_sens))
            
                #Antisens
                antisens = sens.complement()
                ATG_indice_antisens = antisens.find('ATG')
                TGA_indice_antisens = antisens.find('TGA')
                TAG_indice_antisens = antisens.find('TAG')
                TAA_indice_antisens = antisens.find('TAA')
                orf_TGA_antisens = antisens[ATG_indice_antisens:TGA_indice_antisens]
                orf_TAG_antisens = antisens[ATG_indice_antisens:TAG_indice_antisens]
                orf_TAA_antisens = antisens[ATG_indice_antisens:TAA_indice_antisens]
                if len(str(orf_TGA_antisens)) >= 30 and len(str(orf_TGA_antisens))%3==0 :
                    orf_list.append(str(orf_TGA_antisens))
                if len(str(orf_TAG_antisens)) >= 30 and len(str(orf_TAG_antisens))%3==0:
                    orf_list.append(str(orf_TAG_antisens))
                if len(str(orf_TAA_antisens)) >= 30 and len(str(orf_TAA_antisens))%3==0:
                    orf_list.append(str(orf_TAA_antisens)) 
            
            nombre_ORF.append(len(orf_list))
            orf.append(orf_list)


#Creation dataframe
df = pandas.DataFrame(columns = ['gene', 'mRNA_name','mRNA_sequence','longueur_mRNA','nombre_ORF','ORF_sequences'])
df['gene'] = pandas.Series(gene)
df['mRNA_name'] = pandas.Series(mRNA_name)
df['mRNA_sequence'] = pandas.Series(mRNA_sequence)
df['longueur_mRNA'] = pandas.Series(longueur_mRNA)
df['nombre_ORF'] = pandas.Series(nombre_ORF)
df['ORF_sequences'] = pandas.Series(orf)

#Export vers csv
pandas.DataFrame.to_csv(df, 'ORF.csv')