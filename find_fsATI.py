"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""


## Import des modules ##
import pandas


## Fonctions ##

#vérifier qu’il existe bien en 3’ une partie de séquence identique à WT
def find_same (mRNA_wildtype, orf_mut):
    if orf_mut == None:
        same = -2
    same = mRNA_wildtype.find(orf_mut[i:])
    return(same)


#déterminer où se trouve la mutation qui a décalé le site d’initialisation de la traduction en 5’
def find_fs (mRNA_wildtype, orf_mut, same):
    if same >= 0:
        fs = mRNA_wildtype.find(orf_mut[i:same],0,same)
    else:
        fs = mRNA_wildtype.find(orf_mut[i:])
    return(fs)


#regarder la présence de la séquence de Kozak à proximité de l'ATI
def find_kosak (orf):
    if orf == None:
        kosak_indice = -2
    else: 
        kosak_indice_A = orf.find('GCCACCATGG')   
        kosak_indice_G = orf.find('GCCGCCATGG')
        if kosak_indice_A > 0:
            kosak_indice = kosak_indice_A
        if kosak_indice_G > 0:
            kosak_indice = kosak_indice_G
        else: 
            kosak_indice = -1
    return(kosak_indice)



### MAIN ##
df = pandas.read_csv('ORF.csv', header=0, index_col=0)

same_sens_list = []
fs_sens_list = []
kosak_sens_list = []

same_antisens_list = []
fs_antisens_list = []
kosak_antisens_list = []

orf_list = []

#Lecture du dataframe ligne par ligne
for i in range(len(df)):

    if (df.iloc[i]['mutation']) == 'WILDTYPE':
        
        wt_sens = df.iloc[i]['mRNA_sens'] 
        wt_antisens = df.iloc[i]['mRNA_antisens']
        
        orf_list_sens = df.iloc[i+1]['ORF_sens']
        orf_list_sens = orf_list_sens.replace('[','')
        orf_list_sens = orf_list_sens.replace(']','')
        orf_list_sens = orf_list_sens.replace("'",'')
        orf_list_sens = orf_list_sens.split(',')

        orf_list_antisens = df.iloc[i+1]['ORF_antisens']
        orf_list_antisens = df.iloc[i+1]['ORF_sens']
        orf_list_antisens = orf_list_antisens.replace('[','')
        orf_list_antisens = orf_list_antisens.replace(']','')
        orf_list_antisens = orf_list_antisens.replace("'","")
        orf_list_antisens = orf_list_antisens.split(',')

        for orf_sens in orf_list_sens: 
            
            #same
            same_sens = find_same(wt_sens,orf_sens)
            same_sens_list.append(same_sens)
            #fs
            fs_sens = find_fs(wt_sens,orf_sens, same_sens)
            fs_sens_list.append(fs_sens)
            #kosac
            kosak_sens = find_kosak(orf_sens)
            kosak_sens_list.append(kosak_sens)

            orf_list.append(orf_list_sens)
        
            
df_sens_fsATI = pandas.DataFrame(columns = ['ORF', 'Pos_mut_fs', 'Pos_mut_ATI', 'Pos_Kosak'])
df_sens_fsATI['Pos_mut_ATI'] = pandas.Series(same_sens_list)
df_sens_fsATI['Pos_mut_fs'] = pandas.Series(fs_sens_list)
df_sens_fsATI['Pos_Kosak'] = pandas.Series(fs_sens_list)
df_sens_fsATI['ORF'] = pandas.Series(orf_list_sens)
print(df_sens_fsATI)


