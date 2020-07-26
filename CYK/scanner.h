/**
 * Title: scanner.h 
 *
 * Description: Responsible for scanning the fasta file 
 *
 * Copyright (c) 2009 
 * Organization: UGA RNA-Informatics GROUP 
 *
 * @author: Tim Shaw gatech(a)uga.edu
 * @Created: 2/02/2009
 * @Updated: 2/02/2009
 */

#ifndef SCANNER_CLASS
#define SCANNER_CLASS
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <string>
#include "fasta.h"

class Scanner {
 public:
  Scanner();
  Scanner(const char* file);
  std::vector<fasta> getFasta();
  
 private:
  std::vector<fasta> _fastaVector;
  std::ifstream infile;

};

#endif


