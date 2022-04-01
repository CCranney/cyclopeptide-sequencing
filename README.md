# cyclopeptide-sequencing

This project is my ongoing response to the challenge issued at the end of chapter 4, Antibiotic Sequencing, of [Bioinformatics Algorithms](https://www.bioinformaticsalgorithms.org/) by Phillip Compeau.<sup>1</sup> I wrote up several introductions in English and Chinese several years ago, they can be found in the [documentation](https://github.com/CCranney/cyclopeptide-sequencing/tree/master/documentation) folder.


## Data Collection

I was unable to locate the data described at the end of the chapter for Tyrocidine B1. However, I found [a paper](https://pubmed.ncbi.nlm.nih.gov/31864964/) that identified a number of cyclopeptides from public data.<sup>2</sup> In table 5 of the supplement, they list these cyclopeptides as well as the GNPS ID that they were found in. Between that and [what I believe is the GNPS workflow they used](https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22DEREPLICATOR%22%7D) ([VarQuest](http://cab.spbu.ru/software/varquest/)<sup>3</sup>), paper link [here](https://www.nature.com/articles/s41467-018-06082-8), I was able to obtain data for a significant portion of the cyclopeptides described.

The output of my GNPS workflow is under `data/gnps_processing/DEREPLICATOR-1ae5b82b-view_significant_unique-main.tsv`. Included in this package are algorithms for identifying cyclopeptides and determining the amino acid mass sequence purely from their SMILES format as provided in the workflow output.

Mass spectrometry and metadata for several of these files is saved in JSON format in `data/cyclopeptide_references.json`. Data from this JSON file can be visually studied in the jupyter notebook `data/cyclopeptide_visuals.ipynb`. Note that not all cyclopeptides described in the Behsaz paper found their way into this analysis. In addition to just not finding all of them, I used a custom function for identifying cyclic peptides, and am debugging specific situations that I excluded from this data compilation (see Issues). I am also only currently using MZXML files, nothing for MGF or MZML files yet.

Regardless, for anyone who wanted data for developing cyclopeptide sequencing algorithms like myself, I'm hoping this repository provides sufficient data and materials to get started.

### References

1.Compeau, P. & Pevzner, P. Bioinformatics Algorithms: An Active Learning Approach. (Active Learning Publishers, 2018).

2.Behsaz, B. et al. De Novo Peptide Sequencing Reveals Many Cyclopeptides in the Human Gut and Other Environments. Cell Syst 10, 99-108.e5 (2020).

3.Mohimani, H. et al. Dereplication of microbial metabolites through database search of mass spectra. Nat Commun 9, 4035 (2018).

## Instructions

I'm building the base functionality to make future algorithm development smoother for myself and anyone with similar interests, but it's pretty bare-bones for the time being. To be updated.

**Run the tests**

```
python -m pytest tests
```

**Convert GNPS mass spec files to JSON format**

Note: you will need to manually enter the path to you GNPS output in `json_generator.py`. You will also need to have already downloaded the data files and have their respective locations saved as a new column in the GNPS output with the header "LocalFile". The GNPS output file needs to be in csv format.

I mostly put things in JSON format to slim down/simplify the data enough to upload to GitHub without too much fuss. Those data files have more scans than just the desired cyclopeptide scan, and I am effectively removing large swaths of unused data this way.

```
python -m data.json_generator
```
