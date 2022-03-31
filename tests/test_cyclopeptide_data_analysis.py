import pytest
from pyteomics import mzxml
from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
import numpy as np


def test_cyclopeptide_data_analysis_approx():
    assert cda.approx(128.5, 128.4, tolerance = 0.1)
    assert not cda.approx(128.5, 128.399999, tolerance = 0.1)
    assert not cda.approx(128.5, 128.600001, tolerance = 0.1)

@pytest.fixture
def surugamide_b_metadata():
    metadata = {
        'Name':'surugamideB',
        'LocalFile':'tests/data_files/surugamideB.mzXML',
        'Scan':'454',
        'Charge':'1',
        'SMILES':'NCCCC[C@H]1NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](C)NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@H](NC(=O)[C@H](NC1=O)[C@H](CC)C)Cc1ccccc1)CC(C)C)[C@H](CC)C)[C@H](CC)C)C(C)C'}
    return metadata

@pytest.fixture
def surugamide_b_expected_masses(): return [128.09, 99.07, 113.08, 71.04, 113.08, 113.08, 147.07, 113.08]

def test_cyclopeptide_data_analysis_convert_smiles_to_amino_acid_mass_sequence(surugamide_b_metadata, surugamide_b_expected_masses):
    masses = cda.convert_smiles_to_amino_acid_mass_sequence(surugamide_b_metadata['SMILES'])
    assert len(surugamide_b_expected_masses) == len(masses)
    for i in range(len(masses)): assert cda.approx(surugamide_b_expected_masses[i], masses[i], tolerance=0.005)

def test_cyclopeptide_data_analysis_is_unbranched_cyclopeptide_from_smiles(surugamide_b_metadata):
    no_peptide_bonds_sapanisertib_smiles = 'CC(C)N1C2=NC=NC(=C2C(=N1)C3=CC4=C(C=C3)OC(=N4)N)N'
    one_peptide_bond_aspartame = 'COC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(=O)O)N'
    two_peptide_bonds_pubChemCID_24755479 = 'CC(=O)N[C@@H](CCCN=C(N)NC(=O)NC)C(=O)N(C)[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CC(=O)O)C(=O)O'
    linear_peptide_carmaphycin_b_smiles = 'CCCCCC(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCS(=O)(=O)C)C(=O)N[C@@H](CC(C)C)C(=O)[C@]1(CO1)C'
    branched_cyclopeptide_ascidiacyclamide_smiles = 'CC[C@H](C)[C@H]1C2=N[C@@H]([C@H](O2)C)C(=O)N[C@@H](C3=NC(=CS3)C(=O)N[C@H](C4=N[C@@H]([C@H](O4)C)C(=O)N[C@@H](C5=NC(=CS5)C(=O)N1)C(C)C)[C@@H](C)CC)C(C)C'
    assert not cda.is_unbranched_cyclopeptide_from_smiles(no_peptide_bonds_sapanisertib_smiles)
    assert not cda.is_unbranched_cyclopeptide_from_smiles(one_peptide_bond_aspartame)
    assert not cda.is_unbranched_cyclopeptide_from_smiles(two_peptide_bonds_pubChemCID_24755479)
    assert not cda.is_unbranched_cyclopeptide_from_smiles(linear_peptide_carmaphycin_b_smiles)
    assert not cda.is_unbranched_cyclopeptide_from_smiles(branched_cyclopeptide_ascidiacyclamide_smiles)
    assert cda.is_unbranched_cyclopeptide_from_smiles(surugamide_b_metadata['SMILES'])

def test_cyclopeptide_reference_initialization(surugamide_b_metadata):
    # account for .mgf
    with mzxml.read(surugamide_b_metadata['LocalFile'], use_index=True) as spectra: expected_dict = spectra.get_by_id(surugamide_b_metadata['Scan'])
    expected_dict.update(surugamide_b_metadata)
    expected_dict['m/z array'] = [float(x) for x in expected_dict['m/z array']]
    expected_dict['intensity array'] = [float(x) for x in expected_dict['intensity array']]
    for key, value in expected_dict.items():
        if isinstance(value, np.generic):
            expected_dict[key] = value.item()

    observed_dict = cda.upload_spectrum_dictionary(surugamide_b_metadata)
    assert expected_dict == observed_dict

    metadata_without_LocalFile = surugamide_b_metadata.copy()
    del metadata_without_LocalFile['LocalFile']
    with pytest.raises(Exception, match = 'metadata dictionary requires the LocalFile, Scan and Charge keys'):
        cda.upload_spectrum_dictionary(metadata_without_LocalFile)

    metadata_without_Scan = surugamide_b_metadata.copy()
    del metadata_without_Scan['Scan']
    with pytest.raises(Exception, match = 'metadata dictionary requires the LocalFile, Scan and Charge keys'):
        cda.upload_spectrum_dictionary(metadata_without_Scan)

    metadata_without_Charge = surugamide_b_metadata.copy()
    del metadata_without_Charge['Charge']
    with pytest.raises(Exception, match = 'metadata dictionary requires the LocalFile, Scan and Charge keys'):
        cda.upload_spectrum_dictionary(metadata_without_Charge)
