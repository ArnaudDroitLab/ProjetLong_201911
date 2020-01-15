"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""


## Import des modules ##
import pandas


## Fonctions ##

#Vérifier qu’il existe bien en 3’ une partie de séquence identique à WT
#Colonne Pos_mut_fs
def find_ATI (mRNA_wildtype, orf_mut):
    ATI = -1
    for i in range(len(orf_mut)):
        ATI = mRNA_wildtype.find(orf_mut[i:])
    return ATI


#Déterminer où se trouve la mutation qui a décalé le site d’initialisation de la traduction en 5’
#Colonne Pos_mut_ATI
def find_fs (mRNA_wildtype, orf_mut, ATI):
    fs = -1
    if ATI > 0:
        for i in range(ATI):
            fs = mRNA_wildtype.find(orf_mut[i:ATI])
    return fs

#Regarder la présence de la séquence de Kozak à proximité de l'ATI
#Colonne Pos_Kosak
def find_kosak (orf): 
    if orf.find('A**ATGG')>=0 or orf.find('G**ATGG')>=0:
        kosak_strength = 3
    if orf.find('A**ATG')>=0 or orf.find('G**ATG')>=0 or orf.find('ATGG'):
        kosak_strength = 2
    if orf.find('ATG')>=0:
        kosak_strength = 1
    else: 
        kosak_strength = -1
    return kosak_strength

### PHYLOP ###
def find_phylop (ATI_pos):
    
    #Lecture du fichier PhyloP
    with open("../PhyloP/phylop_vert.sga", "r") as filin:
        #Extraction des infos
        for phylop_ligne in filin:
            phylop_list = phylop_ligne.split(sep='\t')
            transcript_name = phylop_list[0]
            position_genomique = phylop_list[2]
            score_phylop = int(phylop_list[4].replace('\n',''))+1        

    #Récupération des coordonnées génomiques pour chaque mRNA
    #A la main le temps de trouver une solution pour automatiser
    posgen = []
    phylop_score = -17

    for i in range(208):
        posgen.append(24663290+i)

    #Association position génomique du mRNA à celle du fichier PhyloP
    for i in range(len(position_genomique)):
        if position_genomique[i] == posgen[ATI_pos]:
            phylop_score = score_phylop[i]
        
    return (phylop_score)


############
### MAIN ###




#Lecture de orf.csv
df = pandas.read_csv('ORF.csv', header=0, index_col=0)

#Initiation des listes
ATI_list = []
fs_list = []
orf_list = []
gene_list = []
transcript_list = []
mutation_list = []
sens_list = []
kosak_list = []
phylop_serie = []

#Lecture du dataframe ligne par ligne
for i in range(len(df)):

    if (df.iloc[i]['mutation'] == 'WILDTYPE'):
        
        #Sequences WT de reference
        wt_sens = df.iloc[i]['mRNA_sens']
        wt_antisens = df.iloc[i]['mRNA_antisens']
        
        #Sequences mutee a comparer
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

        ## SENS ##
        for orf_sens in orf_list_sens: 
    
            orf_list.append(orf_sens)
            gene_list.append(df.iloc[i+1]['gene'])
            transcript_list.append(df.iloc[i+1]['transcript'])
            mutation_list.append(df.iloc[i+1]['mutation'])
            sens_list.append('sens')

            #ATI
            ATI_sens = find_ATI(wt_sens,orf_sens)
            ATI_list.append(ATI_sens)
            #fs
            fs_sens = find_fs(wt_sens,orf_sens, ATI_sens)
            fs_list.append(fs_sens)
            #Kosac
            kosak_sens = find_kosak(orf_sens[ATI_sens-15:ATI_sens+15])
            kosak_list.append(kosak_sens)
            #PhyloP
            phylop_sens = find_phylop(ATI_sens)
            phylop_serie.append(phylop_sens)
        
        ## ANTISENS ##
        for orf_antisens in orf_list_antisens: 
            
            orf_list.append(orf_antisens)
            gene_list.append(df.iloc[i+1]['gene'])
            transcript_list.append(df.iloc[i+1]['transcript'])
            mutation_list.append(df.iloc[i+1]['mutation'])
            sens_list.append('antisens')

            #ATI
            ATI_antisens = find_ATI(wt_antisens,orf_antisens)
            ATI_list.append(ATI_antisens)
            #fs
            fs_antisens = find_fs(wt_antisens,orf_antisens, ATI_antisens)
            fs_list.append(fs_antisens)
            #Kosac
            kosak_antisens = find_kosak(orf_antisens[ATI_antisens-15:ATI_antisens+15])
            kosak_list.append(kosak_antisens)
            #PhyloP
            phylop_antisens = find_phylop(ATI_antisens)
            phylop_serie.append(phylop_antisens)


#Creation dataframe resultats          
df_indices = pandas.DataFrame(columns = ['gene', 'transcript', 'mutation', 'sens', 'ORF', 'Pos_mut_fs', 'Pos_mut_ATI', 'Kosak_strength','PhyloP_score'])
df_indices['Pos_mut_ATI'] = pandas.Series(ATI_list)
df_indices['Pos_mut_fs'] = pandas.Series(fs_list)
df_indices['ORF'] = pandas.Series(orf_list)
df_indices['gene'] = pandas.Series(gene_list)
df_indices['transcript'] = pandas.Series(transcript_list)
df_indices['mutation'] = pandas.Series(mutation_list)
df_indices['sens'] = pandas.Series(sens_list)
df_indices['Kosak_strength'] = pandas.Series(kosak_list)
df_indices['PhyloP_score'] = pandas.Series(phylop_serie)

#Selectionne les lignes en fonction des indices positifs pour Pos_mut_fs et Pos_mut_ATI
df_ATI = df_indices.loc [ df_indices['Pos_mut_ATI']>0 ]
df_fsATI = df_ATI.loc [ df_ATI['Pos_mut_fs']>0 ]


print('\nfsATI trouvées.')
print('\nValidation des fsATI:')
print('***Forces de Kosak calculées.')
print('***Scores PhyloP recuperés.')

#Export vers csv
pandas.DataFrame.to_csv(df_fsATI, 'fsATI.csv')
print('\nFichier fsATI.csv généré.\n')
print('Job done.\n')