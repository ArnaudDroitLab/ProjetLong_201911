"""
@author : Julie Bogoin
Master 2 BIB Paris Diderot - 2019-2020
Projet Long
"""


## Import des modules ##
import pandas

## Fonction ##

# Selectionner les séquences avec un score Kosak=3 et PholoP>=2
def select (fsATI_file):
    #Creation d'un dataframe vide
    df_validated = pandas.DataFrame()

    #Lecture de fsATI.csv
    df_fsATI = pandas.read_csv(fsATIfile, header=0, index_col=0)

    #Tri
    for i in range(len(df_fsATI)):

        if df_fsATI.iloc[i]['Kosak_strength'] == 3 and df_fsATI.iloc[i]['PhyloP score']>=2:
            df_validated.iloc[i] = df_fsATI.iloc[i]
    
    return df_validated



### MAIN ##
df_validated = select('fsATI.csv')
print('\nSelections des fsATI réalisée.\n')

#Export vers csv
pandas.DataFrame.to_csv(df_validated, 'fsATI_validated.csv')
print('\nFichier fsATI_validated.csv généré.\n')
print('Job done.\n')





