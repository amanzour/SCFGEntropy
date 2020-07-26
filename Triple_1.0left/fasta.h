
#ifndef FASTA_CLASS
#define FASTA_CLASS

#include <iostream>
#include <vector>
#include <string>

class fasta {
    public:
	fasta(std::string str1, std::string str2);
	std::string getTitle();
    std::string getSequence();
	void setTitle(std::string str);
	void setSequence(std::string str);
    private:
        std::string title;
        std::string sequence;
};

#endif



