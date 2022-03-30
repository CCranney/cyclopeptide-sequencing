from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import numpy as np
from itertools import combinations
from collections import Counter

peptide_bond_representation = Chem.MolFromSmarts('NCC(=O)N')

class CyclopeptideReference:
    """
    CyclopeptideReference: A data representation of a known cyclopeptide identified on GNPS.
    Includes data values such as name, mass spec dataset, charge, and amino acid mass sequence.
    """
    def __init__(self, ftp, scan, name, charge, smiles): pass

    @staticmethod
    def approx( x, y, tolerance=0.0100000001 ):
        return abs(x-y) <= tolerance

    @staticmethod
    def convert_smiles_to_amino_acid_mass_sequence(smiles):
        print(smiles)
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
        return masses

    @staticmethod
    def is_unbranched_cyclopeptide_from_smiles(smiles):
        peptide_bonds = Chem.MolFromSmiles(smiles).GetSubstructMatches(peptide_bond_representation)
        if len(peptide_bonds) < 3: return False
        nitrogen_indices = [peptide_bond[0] for peptide_bond in peptide_bonds] + [peptide_bond[-1] for peptide_bond in peptide_bonds]
        for _,i in Counter(nitrogen_indices).items():
            if i != 2: return False
        return True
