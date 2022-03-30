# cyclopeptide-sequencing

This project is my ongoing response to the challenge issued at the end of chapter 4, Antibiotic Sequencing, of Bioinformatics Algorithms by Phillip Compeau^1. I wrote up several introductions in English and Chinese several years ago, they can be found in the `documentation` folder.


## Data Collection

I was unable to locate the data described at the end of the chapter for Tyrocidine B1. However, I found a paper that identified a number of cyclopeptides from public data^2. In table 5 of the supplement, they list these cyclopeptides as well as the GNPS ID that they were found in. Between that and what I believe is the GNPS workflow they used (VarQuest^3), I was able to obtain data for a significant portion of the cyclopeptides described.

The output of my GNPS workflow is under `data/DEREPLICATOR-1ae5b82b-view_significant_unique-main.tsv`. Included in this package are algorithms for identifying cyclopeptides and determining the amino acid mass sequence purely from their SMILES format as provided in the workflow output.

I will be working on providing a single datafile that contains the spectra and other important features of each cyclopeptide for future reference.

https://www.bioinformaticsalgorithms.org/
https://pubmed.ncbi.nlm.nih.gov/31864964/
http://cab.spbu.ru/software/varquest/
https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22DEREPLICATOR%22%7D
https://www.nature.com/articles/s41467-018-06082-8

## Instructions

I'm building the base functionality to make future algorithm development smoother for myself and anyone with similar interests, but it's pretty bare-bones for the time being. To be updated.

**Run the tests**

```
python -m pytest tests
```
