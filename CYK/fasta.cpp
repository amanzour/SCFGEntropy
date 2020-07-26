/**
 * Title: fasta.cpp
 *
 * Description: Fasta class object
 *
 * Copyright (c) 2009 
 * Organization: UGA RNA-Informatics GROUP 
 *
 * @author: Tim Shaw gatech(a)uga.edu
 * @Created: 2/02/2009
 * @Updated: 2/02/2009
 */

#include "fasta.h"
std::string fasta::getTitle() {return title;}
std::string fasta::getSequence() {return sequence;}

void fasta::setTitle(std::string str) {title = str;}
void fasta::setSequence(std::string str) {sequence = str;}
fasta::fasta(std::string str1, std::string str2) {
    title = str1;
    sequence = str2;
}


