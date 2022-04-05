from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

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
            self.intensities = data['intensity array']
            self.charge = data['Charge']
            self.precMz = data['precursorMz'][0]['precursorMz']
            self.sequence = cda.convert_smiles_to_amino_acid_mass_sequence(data['SMILES'])
            self.metadata = data
        else: #for json serialization
            self.__dict__.update(kwargs)

    @classmethod
    def with_validation(cls, data):

        sequence = cda.convert_smiles_to_amino_acid_mass_sequence(data['SMILES'])
        if sequence == 'sequence mass does not equal total mass': return None#print(data['SMILES']); return None
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

    @staticmethod
    def plot_spec(ax, SPECTRA, COLOR):
        ax.vlines(SPECTRA.columns, np.repeat(0, len(SPECTRA.columns)), SPECTRA, colors=COLOR)

    # test this
    def convert_mz_to_mass(self): return [(x-self.charge)*self.charge for x in self.mz if (x*self.charge)-self.charge < sum(self.sequence)-50]

    def draw_peaks(self, ax=None):
        ax = ax or plt.gca()

        theoretical_peaks = cda.make_theoretical_spectrum(self.sequence)
        expected_masses = self.convert_mz_to_mass()
        i, j = 0,0
        matched_peak_indices = []
        tolerance = 0.01
        while i < len(theoretical_peaks) and j < len(expected_masses):
            if cda.approx(theoretical_peaks[i], expected_masses[j], tolerance=tolerance): matched_peak_indices.append(j); i+=1
            elif theoretical_peaks[i] > expected_masses[j]: j+=1
            else: i+=1
            if (i >= len(theoretical_peaks) or j >= len(expected_masses)) and len(matched_peak_indices) < 10:

                if tolerance < 0.1: newTolerance = tolerance + 0.01
                else: newTolerance = tolerance + 0.1

                print(f'Found {len(matched_peak_indices)} peaks at {tolerance} tolerance; increasing tolerance to {newTolerance}')
                if newTolerance >=1.0: print('ABORT: tolerance too high'); break
                tolerance = newTolerance

                i, j = 0,0
                matched_peak_indices = []
            elif i >= len(theoretical_peaks) or j >= len(expected_masses):
                print(f'Found {len(matched_peak_indices)} peaks at {tolerance} tolerance')
            if len(expected_masses) <= 10: print (f'ABORT: only {len(expected_masses)} peaks in the spectrum!'); break
        noise_peak_indices = sorted(set([i for i in range(len(expected_masses))]) - set(matched_peak_indices))
        matches_df = pd.DataFrame(data=[[self.intensities[i] for i in matched_peak_indices]], columns=[self.mz[i] for i in matched_peak_indices])
        noise_df = pd.DataFrame(data=[[self.intensities[i] for i in noise_peak_indices]], columns=[self.mz[i] for i in noise_peak_indices])

        ax.set_ylim(0,5000)
        CyclopeptideReference.plot_spec(ax, matches_df, 'red')
        CyclopeptideReference.plot_spec(ax, noise_df, 'grey')

    def draw_doughnut_plot(self, ax=None):
        ax = ax or plt.gca()

        rounded_sequences = [round(x,1) for x in self.sequence]
        cm = plt.get_cmap("tab20")
        ss_sequences = sorted(set(rounded_sequences))
        available_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        if len(ss_sequences) > len(available_colors): available_colors += [cm(i/len(ss_sequences)-len(available_colors)) for i in range(len(ss_sequences)-len(available_colors))]
        color_dict = {ss_sequences[i]:available_colors[i] for i in range(len(ss_sequences))}
        colors = [color_dict[x] for x in rounded_sequences]

        ax.pie(rounded_sequences, labels=rounded_sequences, startangle=90, wedgeprops = { 'linewidth' : 7, 'edgecolor' : 'white' }, colors=colors)

        doughnut_hole=plt.Circle( (0,0), 0.7, color='white')
        ax.add_artist(doughnut_hole)

    def draw_doughnut_and_peaks(self):
        fig, (ax1, ax2) = plt.subplots(2)
        self.draw_doughnut_plot(ax1)
        self.draw_peaks(ax2)
        plt.show()
