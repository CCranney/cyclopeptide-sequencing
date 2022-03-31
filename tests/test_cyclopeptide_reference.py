import pytest
from cyclopeptide_representations.CyclopeptideReference import CyclopeptideReference
from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
from tests.test_cyclopeptide_data_analysis import surugamide_b_metadata, surugamide_b_expected_masses

def test_cyclopeptide_reference_initialization(surugamide_b_metadata, surugamide_b_expected_masses):
    expected_dict = cda.upload_spectrum_dictionary(surugamide_b_metadata)
    expected_sequence = cda.convert_smiles_to_amino_acid_mass_sequence(surugamide_b_metadata['SMILES'])
    cr = CyclopeptideReference.with_validation(expected_dict)

    assert cr.metadata == expected_dict

    assert len(surugamide_b_expected_masses) == len(cr.sequence)
    for i in range(len(cr.sequence)): assert cda.approx(surugamide_b_expected_masses[i], cr.sequence[i], tolerance=0.005)

    # from an error I'm investigating (see GitHub issues)
    dummy_dictionary_with_bad_smiles = {'SMILES':'CC(C)C1C(=O)NC(C(C)C)C(=O)N(C)C(Cc2ccc(cc2)O)C(=O)NC(Cc2ccc(cc2)O)C(=O)N2CCCC2C(=O)N1'}
    cr = CyclopeptideReference.with_validation(dummy_dictionary_with_bad_smiles)
    assert not cr
