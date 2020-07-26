
#include "scanner.h"

std::vector<fasta> Scanner::getFasta() {
    return _fastaVector;
}

Scanner::Scanner() {}
Scanner::Scanner(const char* file) {
    std::string str, title, sequence, empty;
    
    infile.open(file);

    getline(infile, str);
    while (!infile.eof()) {        		
	if (str.find(">", 0) == 0) {
	    int n = str.find("\r");
	    if (n != std::string::npos) {
		title = str.substr(1, str.length() - 1);
	    } else {
	        title = str.substr(1, str.length());
	    }
	}
	bool eofSeq = false; // finds the end of the sequence
        getline(infile, str);
	while (!infile.eof() && !eofSeq) {
	    
    	    if (str.find(">", 0) == 0) {		
	        eofSeq = true;
	    } else {	
	        int n = str.find("\r");
	        if (n != std::string::npos) {
	            sequence += str.substr(0, str.length() - 1);	    
	        } else {
	            sequence += str;	    
	        }
	        getline(infile, str);
	    }
	}
	if (str.find(">", 0) != 0) {		
	        int n = str.find("\r");
	        if (n != std::string::npos) {
	            sequence += str.substr(0, str.length() - 1);	    
	        } else {
	            sequence += str;	    
	        }
	}
	fasta f(title, sequence);
	_fastaVector.push_back(f);
	sequence = "";
    }
    infile.close();	
}

/**/



