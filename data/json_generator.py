import pandas as pd
from cyclopeptide_representations.CyclopeptideReference import CyclopeptideReference
from cyclopeptide_representations.CyclopeptideDataAnalysis import CyclopeptideDataAnalysis as cda
import json

def generate_cyclopeptide_reference_json_from_gnps_data(gnps_input_file, output_file):
    gnps = pd.read_csv(gnps_input_file)
    cyclopeptide_references = []
    for idx, row in gnps.iterrows():
        if '.mzXML' in row['SpecFile']:
            d = cda.upload_spectrum_dictionary(row.to_dict())
            cr = CyclopeptideReference.with_validation(d)
            if cr: cyclopeptide_references.append(cr.to_json())

    with open(output_file, 'w') as f:
        json.dump(cyclopeptide_references, f)

generate_cyclopeptide_reference_json_from_gnps_data('data/gnps_processing/dereplicator_filtered_cyclos_with_local_files.csv', 'data/gnps_processing/cyclopeptide_references.json')
