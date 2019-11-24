"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""


## Import des modules ##
import os
import pandas
import re
from Bio import SeqIO


## Fonctions ##
def listdir_nohidden(path):
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f



## MAIN ##
fasta_list = listdir_nohidden('./mRNA_fasta')
path = './mRNA_fasta/'
dico={}
annovar_line = []
transcript_name = []
mutation_name = []
mRNA_sens = []
mRNA_antisens = []
gene = []
longueur_mRNA = []
nombre_ORF_sens = []
nombre_ORF_antisens = []
ORF_sens = []
ORF_antisens = []

for file_name in fasta_list:
    gene_name = (file_name.split('_')[0])
    dico[gene_name+'_mRNA'] = []
    dico[gene_name+'_mRNA'].append(SeqIO.to_dict(SeqIO.parse(path+file_name, 'fasta')))
    
    for id in dico[gene_name+'_mRNA']:
        
        for key in id.keys():
            
            if key[5] == 'N':
                annovar_line.append(key[0:5])
            else:
                annovar_line.append(key[0:6])
                
            mRNA_name = re.sub('line(\d)+','',key)
            transcript_name.append(mRNA_name[0:9])
            mutation_name.append(mRNA_name[9:])
        
        for value in id.values():
            orf_sens_list = []
            orf_antisens_list = []
            mRNA_sens.append(str(value.seq))
            mRNA_antisens.append(str(value.seq.complement()))
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
                    orf_sens_list.append(str(orf_TGA_sens))
                if len(str(orf_TAG_sens)) >= 30 and len(str(orf_TAG_sens))%3==0:
                    orf_sens_list.append(str(orf_TAG_sens))
                if len(str(orf_TAA_sens)) >= 30 and len(str(orf_TAA_sens))%3==0:
                    orf_sens_list.append(str(orf_TAA_sens))
            
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
                    orf_antisens_list.append(str(orf_TGA_antisens))
                if len(str(orf_TAG_antisens)) >= 30 and len(str(orf_TAG_antisens))%3==0:
                    orf_antisens_list.append(str(orf_TAG_antisens))
                if len(str(orf_TAA_antisens)) >= 30 and len(str(orf_TAA_antisens))%3==0:
                    orf_antisens_list.append(str(orf_TAA_antisens)) 
            
            nombre_ORF_sens.append(len(orf_sens_list))
            nombre_ORF_antisens.append(len(orf_antisens_list))
            ORF_sens.append(orf_sens_list)
            ORF_antisens.append(orf_antisens_list)


#Creation dataframe
df = pandas.DataFrame(columns = ['gene','annovar_line','transcript','mutation','mRNA_sens', 'mRNA_antisens','longueur_mRNA','nombre_ORF_sens','ORF_sens','nombre_ORF_antisens','ORF_antisens'])
df['gene'] = pandas.Series(gene)
df['annovar_line'] = pandas.Series(annovar_line)
df['transcript'] = pandas.Series(transcript_name)
df['mutation'] = pandas.Series(mutation_name)
df['mRNA_sens'] = pandas.Series(mRNA_sens)
df['mRNA_antisens'] = pandas.Series(mRNA_antisens)
df['longueur_mRNA'] = pandas.Series(longueur_mRNA)
df['nombre_ORF_sens'] = pandas.Series(nombre_ORF_sens)
df['nombre_ORF_antisens'] = pandas.Series(nombre_ORF_antisens)
df['ORF_sens'] = pandas.Series(ORF_sens)
df['ORF_antisens'] = pandas.Series(ORF_antisens)

print('\nLes mRNA ont ete traites.')

#Export vers csv
pandas.DataFrame.to_csv(df, 'ORF.csv')
print('Le fichier ORF.csv a été généré.\n')
print('\nJob done.\n')
