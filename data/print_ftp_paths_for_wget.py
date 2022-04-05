import pandas as pd
from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda

#read file
filePath = 'data/gnps_processing/ProteoSAFe-DEREPLICATOR-9a43aa94-view_significant_unique/DEREPLICATOR-9a43aa94-view_significant_unique-main.tsv'
df = pd.read_csv(filePath, sep='\t')

# remove non-cyclopeptides
remove = []
smiles = list(df['SMILES'])
for i in range(len(smiles)):
    if not cda.is_unbranched_cyclopeptide_from_smiles(str(smiles[i])): remove.append(i)
df = df[~df.index.isin(remove)]

# set local file paths
ftpPaths = df['SpecFile']
fileNames = [x.split('/')[-1] for x in ftpPaths]
fileDir = 'data/gnps_files/'
df['LocalFile'] = [fileDir + file for file in fileNames]

# print lines that can be fed into `wget -i <txt list of files>`
for path in ftpPaths: print('ftp://massive.ucsd.edu/' + path)

# save filtered table
df.to_csv('data/gnps_processing/dereplicator_filtered_cyclos_with_local_files.csv', index=False)
