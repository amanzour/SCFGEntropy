
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


