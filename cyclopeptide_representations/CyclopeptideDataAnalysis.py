from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import numpy as np
from itertools import combinations
from collections import Counter
import pandas as pd
import json
from pyteomics import mzxml

peptide_bond_representation = Chem.MolFromSmarts('NCC(=O)N')

class CyclopeptideDataAnalysis:
    """
    CyclopeptideDataAnalysis: Various functions useful in analyzing cyclopeptides.
    """
    def __init__(self, metadata):
        self.upload_spectrum_dictionary(metadata)

    @staticmethod
    def approx( x, y, tolerance=0.0100000001 ):
        return abs(x-y) <= tolerance

    @staticmethod
    def convert_smiles_to_amino_acid_mass_sequence(smiles):
        m = Chem.MolFromSmiles(smiles)

        am = np.array(Chem.GetAdjacencyMatrix(m))

        masses = []
        for peptide_bond in m.GetSubstructMatches(peptide_bond_representation):

            alpha = peptide_bond[1]
            nitrogens = set([peptide_bond[0],peptide_bond[-1]])

            aa_atom_idx = set([alpha])
            set2 = set()
            while aa_atom_idx != set2:
                set2 = aa_atom_idx.copy()
                temp_am = am[:, list(aa_atom_idx)]
                aa_atom_idx = set(np.where(temp_am==1)[0]) | aa_atom_idx
                aa_atom_idx -= nitrogens
            aa_atom_idx.add(peptide_bond[0])

            bonds = []
            for i,j in combinations(aa_atom_idx, 2):
                b = m.GetBondBetweenAtoms(int(i),int(j))
                if b: bonds.append(b.GetIdx())

            tempMol = Chem.PathToSubmol(m, bonds)
            mass = ExactMolWt(tempMol)
            masses.append(mass)

        # TODO: Investigate why this error is occurring. Recorded as an issue on GitHub.
        if CyclopeptideDataAnalysis.approx(sum(masses), ExactMolWt(m)): return masses
        else: return 'sequence mass does not equal total mass'

    @staticmethod
    def is_unbranched_cyclopeptide_from_smiles(smiles):
        if Chem.MolFromSmiles(smiles): peptide_bonds = Chem.MolFromSmiles(smiles).GetSubstructMatches(peptide_bond_representation)
        else: return False
        if len(peptide_bonds) < 3: return False
        nitrogen_indices = [peptide_bond[0] for peptide_bond in peptide_bonds] + [peptide_bond[-1] for peptide_bond in peptide_bonds]
        for _,i in Counter(nitrogen_indices).items():
            if i != 2: return False
        return True

    @staticmethod
    def upload_spectrum_dictionary(metadata):
        if not all(key in metadata for key in ['LocalFile','Scan','Charge']): raise Exception('metadata dictionary requires the LocalFile, Scan and Charge keys')
        print(metadata['LocalFile'])
        with mzxml.read(metadata['LocalFile'], use_index=True) as spectra: spectrum = spectra.get_by_id(str(metadata['Scan']))
        spectrum['m/z array'] = [float(x) for x in spectrum['m/z array']]
        spectrum['intensity array'] = [float(x) for x in spectrum['intensity array']]
        data = metadata
        data.update(spectrum)
        for key, value in data.items():
            if isinstance(value, np.generic):
                data[key] = value.item()
        return data

    @staticmethod
    def make_theoretical_spectrum(sequence):
        ts = []
        for i in range( len(sequence) ):
            newPeptide = list(sequence[i:]) + list(sequence[0:i])
            for j in range( 1,len(newPeptide) ):
                ts.append(sum(newPeptide[0:j]))
        return sorted(set([round(x,2) for x in ts]))
