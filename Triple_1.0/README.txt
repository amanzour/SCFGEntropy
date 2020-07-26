Triple (Version 1.0)

GENERAL INFORMATION
-------------------

Triple 1.1 is an utility capable of calculating entropies of given RNA sequences (consisting of letters: A, C, G, T, U).

Triple 1.1 is implemented in C++.

C++ compiler: g++ (version 4.1.2 or above)

The software has been tested on the Linux platforms Ubuntu 9.04 or above and Redhat Enterprise 5.5. Any OS of Linux with g++ 4.1.2 or above may work for Triple 1.0.

Authors:Triple (Version 1.2)

GENERAL INFORMATION
-------------------

Triple 1.2 is an utility capable of calculating entropies of given RNA sequences (consisting of letters: A, C, G, T, U).

Triple 1.2 is implemented in C++.

C++ compiler: g++ (version 4.1.2 or above)

The software has been tested on the Linux platforms Ubuntu 9.04 or above and Redhat Enterprise 5.5. Any OS of Linux with g++ 4.1.2 or above may work for Triple 1.0.

Authors:
	Yingfeng Wang and Pooya Shareghi
	Department of Computer Science, University of Georgia, Athens, GA 30602
	Amir Manzour and Timothy I. Shaw
	Institute of Bioinformatics, University of Georgia, Athens, GA 30602

Copyright (C) 2010 RNA-Informatics@UGA. All Rights Reserved.

WHAT IS CONTAINED IN THIS PACKAGE
---------------------------------

This package is for the Linux platform. The package contains Triple 1.2 program.


INSTALLATION
---------------------------------
1. Compile Triple 1.2

	Type 'make' to build four executable program Triple.

2. Re-Compile Triple 1.2

	Type 'make clean' to remove the executable program Triple, and then type 'make'.


USAGE
-----
1.Parameters for Triple 1.2 program

	-s seqfile_name (required)
	-g grammarfile_name (required)
	-w window_size (default window size is the length of the whole sequence)
	-b base_pairmatrix_file_name (required)
	-n nucleotide_vector_file_name (required)
        -p output_file_for_storing_all_pij
	-a indicate test all sequences in the file
		(without it, program just tests the first sequence)
	-t threshold set the threshold for outputting pij matrix
		(-p should be specified if -t is specified)
	-h get help

2. Base pair matrix

	Base pair matrix is a 5x5 matrix, with real number entries.
	Both orders of rows and columns are _ACGU.
	In Triple 1.2, bulge '_' is not used. Here is an example.

 	0.00 0.00 0.00 0.00 0.00
 	0.00 0.00 0.00 0.00 0.17
 	0.00 0.00 0.00 0.25 0.00
 	0.00 0.00 0.25 0.00 0.08
 	0.00 0.17 0.00 0.08 0.00

3. Nucleotide Vector

        Nucleotide vector is a 1x4 vector containing nucleotide frequencies
used by the model. The order is ACGU. Here is an example
        0.25 0.25 0.25 0.25
4. SCF Grammar

	Non-terminals should be represented with delimiter pair '<' and '.\>';
	<Start> is the starting non-terminal.
	Colon ':' denotes the rewriting symbol "->".

	Triple 1.2 accepts 6 types of generic production rules as follows:

	<X>->a
	<X>->a<Y>
	<X>-><Y>a
	<X>->a<Y>b
	<X>->a<Y>b<Z>
	<X>-><Y>a<Z>b

	where <X>, <Y> and <Z> are non-terminals while a and b are terminals.

	The probabilities of corresponding rules are given the same line
	and separated by at least one space.

	The grammar file is ended with symbol '#'.

	Note that there is no space in the rule part.

3. Examples
-------------------

Here is a completed example of grammar,

	<Start>:a     0.499
	<Start>:a<Start>    0.499
	<Start>:a<H>b<Start>  0.001
	<Start>:a<H>b   0.001
	<H>:a<H>b   0.103
	<H>:a<Y>b   0.897
	<Y>:a<Start>b   1
	#

Examples of commands

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt -n nvector.txt -p pijtest.txt

The above command will calculate entropy of the first sequence in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. All pij values will be saved in pijtest.txt.

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt -n nvector.txt

The above command will calculate entropy of the first sequence in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. No  pij values will be outputted.

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt -n nvector.txt -a

The above command will calculate entropy of all sequences in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. No pij values will be outputted.

-s H1.txt -g grammartest_new.txt -b oribpmatrix.txt -n nvector.txt -a

	Yingfeng Wang and Pooya Shareghi
	Department of Computer Science, University of Georgia, Athens, GA 30602
	Amir Manzour and Timothy I. Shaw 
	Institute of Bioinformatics, University of Georgia, Athens, GA 30602

Copyright (C) 2010 RNA-Informatics@UGA. All Rights Reserved.

WHAT IS CONTAINED IN THIS PACKAGE
---------------------------------

This package is for the Linux platform. The package contains Triple 1.0 program.


INSTALLATION
---------------------------------
1. Compile Triple 1.0

	Type 'make' to build four executable program Triple.

2. Re-Compile Triple 1.0

	Type 'make clean' to remove the executable program Triple, and then type 'make'.


USAGE
-----
1.Parameters for Triple 1.0 program 

	-s seqfile_name (required)
	-g grammarfile_name (required)
	-w window_size (default window size is the length of the whole sequence)
	-b base_pairmatrix_file_name (required)
	-p output_file_for_storing_all_pij
	-a indicate test all sequences in the file 
		(without it, program just tests the first sequence)
	-t threshold set the threshold for outputting pij matrix 
		(-p should be specified if -t is specified)
	-h get help

2. Base pair matrix

	Base pair matrix is a 5x5 matrix, with real number entries. 
	Both orders of rows and columns are _ACGU. 
	In Triple 1.0, bulge '_' is not used. Here is an example. 

 	0.00 0.00 0.00 0.00 0.00
 	0.00 0.00 0.00 0.00 0.17
 	0.00 0.00 0.00 0.25 0.00
 	0.00 0.00 0.25 0.00 0.08 
 	0.00 0.17 0.00 0.08 0.00

2. SCF Grammar

	Non-terminals should be represented with delimiter pair '<' and '.\>'; 
	<Start> is the starting non-terminal. 
	Colon ':' denotes the rewriting symbol "->".

	Triple 1.0 accepts 6 types of generic production rules as follows:

	<X>->a
	<X>->a<Y>
	<X>-><Y>a
	<X>->a<Y>b 
	<X>->a<Y>b<Z>
	<X>-><Y>a<Z>b

	where <X>, <Y> and <Z> are non-terminals while a and b are terminals.

	The probabilities of corresponding rules are given the same line 
	and separated by at least one space. 

	The grammar file is ended with symbol '#'. 

	Note that there is no space in the rule part. 

3. Examples
-------------------

Here is a completed example of grammar,

	<Start>:a     0.499
	<Start>:a<Start>    0.499
	<Start>:a<H>b<Start>  0.001
	<Start>:a<H>b   0.001
	<H>:a<H>b   0.103
	<H>:a<Y>b   0.897
	<Y>:a<Start>b   1
	#

Examples of commands

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt -p pijtest.txt

The above command will calculate entropy of the first sequence in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. All pij values will be saved in pijtest.txt.

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt 

The above command will calculate entropy of the first sequence in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. No  pij values will be outputted.

./Triple  -s seq1_out.txt -g grammartest.txt -b oribpmatrix.txt -a

The above command will calculate entropy of all sequences in seq1_out.txt with grammar in  grammartest.txt and base pair matrix in oribpmatrix.txt. No pij values will be outputted.

-s seq1_H1.txt -g grammartest.txt -b oribpmatrix.txt -n nvector.txt