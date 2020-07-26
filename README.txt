There are a total of Four programs in this folder
Triple_1.0 : Original Huynen entropy (not-normalized version).
Triple_1.1 : This is the base-pairing entropy (not-normalized version).
Triple_1.2 This is the Structural Entropy (not-normalized version).
CYK outputs the optimum structure for each sequence
Triple_1.*left are the above programs using an implementation of leftmost
derivation
All above programs constraint the structures to have minimum loop of 2
nucleotides.
All above programs are parameterized so we can input the following:
1. sequences
2. grammar
3. base-pairing matrix
4. nucleotide vector

Example for any of the three entropy programs:
./Triple -s H1.txt -g grammartest_new.txt -b oribpmatrix.txt -n nvector.txt -a
