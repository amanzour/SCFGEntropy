

#include "fasta.h"
std::string fasta::getTitle() {return title;}
std::string fasta::getSequence() {return sequence;}

void fasta::setTitle(std::string str) {title = str;}
void fasta::setSequence(std::string str) {sequence = str;}
fasta::fasta(std::string str1, std::string str2) {
    title = str1;
    sequence = str2;
}


