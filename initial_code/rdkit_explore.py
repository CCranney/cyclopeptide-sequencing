from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import numpy as np
from itertools import combinations
from pyteomics import mzxml, mgf, mass
import matplotlib.pyplot as plt
import pandas as pd

def TheoreticalSpectrumChain( chain ):
    ts = []
    for i in range( len(chain) ):
        newPeptide = list(chain[i:]) + list(chain[0:i])
        for j in range( 1,len(newPeptide) ):
            ts.append(sum(newPeptide[0:j]))
    return ts


def plot_spec(SPECTRA, COLOR):
    plt.vlines(SPECTRA.columns, np.repeat(0, len(SPECTRA.columns)), SPECTRA, colors=COLOR)

def plot_line(v, m):
    plt.hlines(y=0.7,xmin=0,xmax=m,colors='red')

########################################
# Convert SMILES to masses
########################################
#'''
tyroB1 = 'CC(C)C[C@H]1C(=O)N[C@@H](C(=O)N2CCC[C@H]2C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCCCN)C(C)C)CC3=CC=C(C=C3)O)CCC(=O)N)CC(=O)N)CC4=CC=CC=C4)CC5=CNC6=CC=CC=C65)CC7=CC=CC=C7'
suruB =  'NCCCC[C@H]1NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](C)NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@H](NC(=O)[C@H](NC1=O)[C@H](CC)C)Cc1ccccc1)CC(C)C)[C@H](CC)C)[C@H](CC)C)C(C)C'
#m = Chem.MolFromSmiles(tyroB1)
m = Chem.MolFromSmiles(suruB)

peptide_bond_representation = Chem.MolFromSmarts('NCC(=O)N')
am = np.array(Chem.GetAdjacencyMatrix(m))
print(Chem.MolFromSmiles('CC(=O)N[C@@H](CCCN=C(N)NC(=O)NC)C(=O)N(C)[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(=O)O)C(=O)O').GetSubstructMatches(peptide_bond_representation))

masses = []
for peptide_bond in m.GetSubstructMatches(peptide_bond_representation):

    alpha = peptide_bond[1]
    nitrogens = set([peptide_bond[0],peptide_bond[-1]])

    aa_atom_idx = set([alpha])
    set2 = set()
    while aa_atom_idx != set2: # while loop repeats until all atoms in the amino acid are accounted for
        set2 = aa_atom_idx.copy()
        temp_am = am[:, list(aa_atom_idx)]
        aa_atom_idx = set(np.where(temp_am==1)[0]) | aa_atom_idx
        aa_atom_idx -= nitrogens # nitrogens from the peptide bond excluded to isolate atoms of this amino acid
    aa_atom_idx.add(peptide_bond[0]) #adding a single nitrogen back in

    bonds = [] # Chem.PathToSubmol requires a list of bonds, not atoms
    for i,j in combinations(aa_atom_idx, 2):
        b = m.GetBondBetweenAtoms(int(i),int(j))
        if b: bonds.append(b.GetIdx())

    tempMol = Chem.PathToSubmol(m, bonds)
    mass = ExactMolWt(tempMol)
    masses.append(mass)
print([round(mass, 2) for mass in masses])
print(sum(masses))

theory = [(x,'x') for x in set(TheoreticalSpectrumChain( masses ))]
print(len(theory))
#'''

########################################
# Read in mzXML
########################################
'''
querySpectraFile = 'data/GNPS_outputs/surugamideB.mzXML'
with mzxml.read(querySpectraFile, use_index=True) as spectra:
    spec = spectra.get_by_id('454')
    print(spec['precursorMz'])
    precMz = float(spec['precursorMz'][0]['precursorMz'])
    exp = [(spec['m/z array'][i], spec['intensity array'][i]) for i in range(len(spec['m/z array'])) if spec['m/z array'][i] < precMz-50]
    exp.sort(key=lambda x:x[1],reverse=True)
    exp = exp[:len(exp)//2]
    exp.sort(key=lambda x:x[0])
    print(spec.keys())
all = sorted(theory + exp)

noise_i = []
matches_i = []
for i in range(1,len(all)):
    if all[i][1] == 'x': continue
    if all[i-1][1] == 'x' and Approx(all[i][0]-1,all[i-1][0],tolerance=0.01): matches_i.append(i)
    else: noise_i.append(i)

def make_spectra_format(peaks):

    print([x[1] for x in peaks])
    print([x[0] for x in peaks])
    print(len(peaks))
    df = pd.DataFrame(data=[[x[1] for x in peaks]], columns=[x[0] for x in peaks])
    return df

plot_spec(make_spectra_format([all[i] for i in sorted(noise_i)]), 'grey')
plot_spec(make_spectra_format([all[i] for i in sorted(matches_i)]), 'red')
plt.show()
'''









'''
print(atoms[0])

for i in range(len(atoms)):
    print(str(i) + ': ' + atoms[i].GetSymbol() + ': ' + str(atoms[i].GetIdx()))

bonds = m.GetBonds()
print()
for i in range(len(bonds)):
    a1 = bonds[i].GetBeginAtom().GetIdx()
    a2 = bonds[i].GetEndAtom().GetIdx()
    bt = bonds[i].GetBondType()
    print(f'i: {i}, a1: {a1}, a2: {a2}, bt: {bt}')

sm = Chem.PathToSubmol(m2, [8,88,])
sm.UpdatePropertyCache()
sm.Debug()
'''
