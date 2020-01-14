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
    df_fsATI = pandas.read_csv(fsATI_file, header=0, index_col=0)

    # #Tri
    df_validated = df_fsATI.loc [ df_fsATI ['Kosak_strength']==3 , df_fsATI['PhyloP score']>=2 ]
    
    return df_validated



### MAIN ##
df_validated = select('fsATI.csv')
print('\nSelections des fsATI réalisée.')

#Export vers csv
pandas.DataFrame.to_csv(df_validated, 'fsATI_validated.csv')
print('Fichier fsATI_validated.csv généré.\n')
print('Job done.\n')





