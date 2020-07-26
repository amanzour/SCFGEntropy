#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include "insideoutside.h"
#include "readgrammar.h"
#include "scidouble.h"

using namespace std;
class entropy
{
public:
	entropy(int nonnum, vector<ruleinfo> thegrammar, vector< vector< int > > malpha,
		vector< vector< int > > mbeta, int maxWinWidth, int pairnum, double threshold_p, 
		int signal_pij, vector< vector<double> > bpmatrix, vector<double> nvector,string pij_filename);
	~entropy(void);
	insideoutside inout;
	double calentropy(int winbeg, int winwidth, const string &oristr);// calculate entropy 
	scidouble segprobability(int beg, int end, const string &window, scidouble p_whole); 
			// calculate the probability of the segment of sequence, Prob(s[i]^s[j]|s)
	//amir
        double maxEntropy(int winwidth);//calculate the maximum entropy of a length
        vector < vector<scidouble> > amirsegprobability(int beg, int end, const string &window, scidouble p_whole);
        vector < vector<scidouble> > amirsegprobability1(int beg, const string &window, scidouble p_whole);
        // calculate the logprobability of the segment of sequence, Prob(s[i]^s[j]|s)
        //amir done
        void plottingPij(vector<scidouble> & p, scidouble psum, int winwidth);
	void outputPij(vector<scidouble> & p, scidouble psum, int winwidth);
private:
	double thethreshold;// for plotting Pij
        scidouble SCFGcount(int N);//recursive function to calculate all SS
	int thepairnum; // the number of base pairs between alpha and beta
	//double ppair(char ntleft, char ntright);//AU, GC, GU
	//double alphapair(int beg, int end, const string &window, int winbeg);//alpha^(i,j)
	//string winstr;// winstr: string covered by the window;
	//string theoristr;//theoristr: reference of original string; 
	double callogodd (scidouble p_item, scidouble p_sum); //calculate p*log(p)
	ofstream f_p; //for outputing pij
	string pij_file; //pij file name
	int pijsignal; //0 means don't output pij
	vector< vector<double> > basepair_matrix;
        vector<double> n_vector;//nucleotide freq.
	//char underthreshold =' ';
	//char overthreshold  ='*';
};
