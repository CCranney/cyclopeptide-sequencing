import pytest

from cyclopeptide_representations.CyclopeptideReference import CyclopeptideReference

def test_cyclopeptide_reference_approx():
    assert CyclopeptideReference.approx(128.5, 128.4, tolerance = 0.1)
    assert not CyclopeptideReference.approx(128.5, 128.399999, tolerance = 0.1)
    assert not CyclopeptideReference.approx(128.5, 128.600001, tolerance = 0.1)

@pytest.fixture
def surugamide_b_smiles(): return 'NCCCC[C@H]1NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](C)NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@H](NC(=O)[C@H](NC1=O)[C@H](CC)C)Cc1ccccc1)CC(C)C)[C@H](CC)C)[C@H](CC)C)C(C)C'

@pytest.fixture
def surugamide_b_expected_masses(): return [128.09, 99.07, 113.08, 71.04, 113.08, 113.08, 147.07, 113.08]

def test_cyclopeptide_reference_convert_smiles_to_amino_acid_mass_sequence(surugamide_b_smiles, surugamide_b_expected_masses):
    masses = CyclopeptideReference.convert_smiles_to_amino_acid_mass_sequence(surugamide_b_smiles)
    assert len(surugamide_b_expected_masses) == len(masses)
    for i in range(len(masses)): assert CyclopeptideReference.approx(surugamide_b_expected_masses[i], masses[i], tolerance=0.005)

def test_cyclopeptide_reference_is_unbranched_cyclopeptide_from_smiles(surugamide_b_smiles):
    no_peptide_bonds_sapanisertib_smiles = 'CC(C)N1C2=NC=NC(=C2C(=N1)C3=CC4=C(C=C3)OC(=N4)N)N'
    one_peptide_bond_aspartame = 'COC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(=O)O)N'
    two_peptide_bonds_pubChemCID_24755479 = 'CC(=O)N[C@@H](CCCN=C(N)NC(=O)NC)C(=O)N(C)[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(=O)O)C(=O)O'
    linear_peptide_carmaphycin_b_smiles = 'CCCCCC(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCS(=O)(=O)C)C(=O)N[C@@H](CC(C)C)C(=O)[C@]1(CO1)C'
    branched_cyclopeptide_ascidiacyclamide_smiles = 'CC[C@H](C)[C@H]1C2=N[C@@H]([C@H](O2)C)C(=O)N[C@@H](C3=NC(=CS3)C(=O)N[C@H](C4=N[C@@H]([C@H](O4)C)C(=O)N[C@@H](C5=NC(=CS5)C(=O)N1)C(C)C)[C@@H](C)CC)C(C)C'
    assert not CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(no_peptide_bonds_sapanisertib_smiles)
    assert not CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(one_peptide_bond_aspartame)
    assert not CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(two_peptide_bonds_pubChemCID_24755479)
    assert not CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(linear_peptide_carmaphycin_b_smiles)
    assert not CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(branched_cyclopeptide_ascidiacyclamide_smiles)
    assert CyclopeptideReference.is_unbranched_cyclopeptide_from_smiles(surugamide_b_smiles)

@pytest.fixture
def surugamide_b_cyclopeptide_reference(surugamide_b_smiles):
    ftp = 'ftp://massive.ucsd.edu/MSV000078937/raw/Samples/ISP2-s-BuOH-P-03_B-5_P1-B-5_01_14292.mzXML'
    scan = '454'
    name = 'surugamide_b'
    charge = 1
    return CyclopeptideReference(ftp, scan, name, charge, surugamide_b_smiles)

def test_cyclopeptide_reference_init(surugamide_b_cyclopeptide_reference):
    # self.these saved variables
    # list of masses equals expected masses
    # spectrum from file saved as dict
    # work on "processing" functions, like removing low-quality or almost precMz peaks, for later
    pass
