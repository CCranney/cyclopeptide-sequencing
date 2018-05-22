# cyclopeptide-sequencing

Please see file "Sequencing Antibiotics Project Summary" for a summary of the challenge of sequencing antibiotics and "Sequencing Antibiotics Algorithm Theory" for a basic explanation of the algorithm being developed. Chinese translations are 测序抗生素项目总结 and 测序抗生素算法理论, respectively. My native tongue is English, so these files will likely contain a number of linguistic errors.

As a basic outline, my algorithm is broken into the following steps:
1. Predict the full cyclopeptide weight and determine candidate amino acids.
2. Line up segments that overlap (using the "flip and add" technique)
3. Build crude candidate cyclopeptide sequences
4. Use the crude candidate cyclopeptide sequences to determine likely distances between candidate amino acids
5. Using frequently occuring distances, build a refined candidate cyclopeptide sequence.

The refined candidate cyclopeptide sequence is not the full cyclopeptide, but most if not all of it's theoretical spectrum would match the theoretical spectrum of the full cyclopeptide. Thus, one could create a significantly large number of new (and accurate) data points that coincide with the full cyclopeptide, improving the quality of data tremendously.
