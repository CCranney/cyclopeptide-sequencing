import pytest
from cyclopeptide_representations.CyclopeptideReference import CyclopeptideReference
from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
from tests.test_cyclopeptide_data_analysis import surugamide_b_metadata, surugamide_b_expected_masses

@pytest.fixture
def surugamide_b_cyclopeptide_reference(surugamide_b_metadata):
    expected_dict = cda.upload_spectrum_dictionary(surugamide_b_metadata)
    cr = CyclopeptideReference.with_validation(expected_dict)
    return cr


def test_cyclopeptide_reference_initialization(surugamide_b_metadata, surugamide_b_cyclopeptide_reference, surugamide_b_expected_masses):
    expected_dict = cda.upload_spectrum_dictionary(surugamide_b_metadata)
    expected_sequence = cda.convert_smiles_to_amino_acid_mass_sequence(surugamide_b_metadata['SMILES'])

    assert surugamide_b_cyclopeptide_reference.metadata == expected_dict

    assert len(surugamide_b_expected_masses) == len(surugamide_b_cyclopeptide_reference.sequence)
    for i in range(len(surugamide_b_cyclopeptide_reference.sequence)): assert cda.approx(surugamide_b_expected_masses[i], surugamide_b_cyclopeptide_reference.sequence[i], tolerance=0.005)

    # from an error I'm investigating (see GitHub issues)
    dummy_dictionary_with_bad_smiles = {'SMILES':'CC(C)C1C(=O)NC(C(C)C)C(=O)N(C)C(Cc2ccc(cc2)O)C(=O)NC(Cc2ccc(cc2)O)C(=O)N2CCCC2C(=O)N1'}
    cr = CyclopeptideReference.with_validation(dummy_dictionary_with_bad_smiles)
    assert not cr
