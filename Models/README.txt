This folder contains a total of 9 models taken from the lightweight SCFG paper by EDDY (http://selab.janelia.org/software/conus/). There are actually three grammar models G4, G5, and G6 each of which are independently trained on three separate datasets rfam5, mixed80 and benchmark (refer to paper for more info.).
Each of the nine Models has three files associated with it:
1. Grammar Rules and their probabilities
2. Base-pair probability matrix
3. Nucleotide probabilities vector

G4=RUN
G5=IVO
G6=BJK
All above grammars have to be constrained on minimum loop of 2 nucleotides.
Here's an example for the naming:
Consider model of grammar 4 and training set of rfam5
corresponding files to be used are:
1. grammarRUNrfam5.txt
2. basepairRUNrfam5.txt
3. nucleotideRUNrfam5.txt

An example of using the structural entropy program using the above would be:
./Triple -s H1.txt -g grammarRUNrfam5.txt -b basepairRUNrfam5.txt -n nucleotideRUNrfam5.txt -a
Note. Certain models had to be converted to usable from by the program (Refer
to each folder for details).

folder vienna contains models from (http://www.ncbi.nlm.nih.gov/pubmed/2219430)
g6 and basic grammar along with various parameter sets.
Slight conversion approximation done.
