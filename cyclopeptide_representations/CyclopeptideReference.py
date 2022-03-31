from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
import pandas as pd
import json

class CyclopeptideReference:
    """
    CyclopeptideReference: A data representation of a known cyclopeptide identified on GNPS.
    Includes data values such as name, mass spec dataset, charge, and amino acid mass sequence.
    """
    def __init__(self, **kwargs):
        if 'Name' in kwargs:
            data = kwargs
            self.name = data['Name']
            self.mz = data['m/z array']
            self.intensities = data['m/z array']
            self.charge = data['Charge']
            self.precMz = data['precursorMz'][0]['precursorMz']
            self.sequence = cda.convert_smiles_to_amino_acid_mass_sequence(data['SMILES'])
            self.metadata = data
        else: #for json serialization
            self.__dict__.update(kwargs)

    @classmethod
    def with_validation(cls, data):

        sequence = cda.convert_smiles_to_amino_acid_mass_sequence(data['SMILES'])
        if sequence == 'sequence mass does not equal total mass': print(data['SMILES']); return None
        return cls(**data)

    def __repr__(self):
        return (
            'CyclopeptideReference:'
            f'\n\tName: {self.name}'
            f'\n\tMass: {sum(self.sequence)}'
            f'\n\tSequence: {self.sequence}'
        )

    def to_json(self): return self.__dict__

    @staticmethod
    def from_json(reference_json): return CyclopeptideReference(**reference_json)

    def draw_peaks(self): pass
