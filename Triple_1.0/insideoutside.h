
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include "readgrammar.h"
#include "scidouble.h"

#define MINDIST  2
using namespace std;

class insideoutside
{
public:
	insideoutside(void);
	~insideoutside(void);
	scidouble alpha(int v, int i, int j); //return result of inside algorithm
	scidouble beta(int v, int i, int j); // return result of outside algorithm
	void alphatable(const string &winstr); //calculate the alpha table with inside algorithm
	void betatable(const string &winstr, int winbeg); //calculate the alpha table with inside algorithm
	void setmaxwidth(int MaxWidth); // set the the maximum width of window
	double ppair(char ntleft, char ntright);//AU, GC, GU
	int mindist;   // the minimum distance between two nucleotides which are paired
	void betaspace(int MaxWidth); // create space for betatable
	void Initilizealpha(void); //set flagalpha to be 0;
	//double ruleprob[30]; //15 is enough for bulge version, but I set more for future
	scidouble basepair(int i, int j, const string &str);
	void setgrammar(int, vector<ruleinfo>, vector< vector< int > >, vector< vector< int > >);// set grammar
	vector<ruleinfo> grammar; // list of rules
	void setbpmatrix(vector< vector<double> > curbpmatrix);
        void setnvector(vector< double > curvector);

private:
	vector< vector< int > > map_alpha; // mapping_alpha;
	vector< vector< int > > map_beta; // mapping_beta;	
	//string nonterm; // list of all non-terminals.
	int nontnum; // # of all non-terminals

	int flagalpha; //0 indicates calculate whole alpha table; 1 indicates calculate the additional column.
	int windowbegin;  // position of the window it is the same as the parameter winbeg used in functions.
	int maxwidth;  // the maximum width of window
	int winwidth; //the width of window (the length of the subsequence)
	vector< vector< vector<scidouble> > > alphamatrix; // alpha matrix
	vector< vector< vector<scidouble> > > betamatrix; // beta matrix
	vector<scidouble> value;   //for initialize vector with two dimensions 
	void clearalpha(void); //clear alpha matrix
	void clearbeta(void);  //clear beta matrix
	scidouble qprob(int posi, const string &qstr); //q(s[i])
	void generatealpha(const string &oristr, int winlength); // calculate alpha matrix
	void generatebeta(const string &winstr); // calculate beta matrix
	int maxvalue(int valA, int valB);// return the maximum value;
	int minvalue(int valA, int valB);// return the minimum value;
	void slidealpha(const string &winstr, int winlength);// get the new table after sliding one nt
	int nt2num(char nt); //A=0, C=1, G=2, U=3;
	void cal_alpha(int v, int i, int j, const string &oristr);
	void cal_beta (int v, int i, int j, const string &winstr);
	void cal_beta_leftmost(int v, int j, const string &winstr); // calculating beta(v, -1, j) it means i=-1
	double prior_matrix[5][5];// base pair matrix
        double prior_vector[4];// nucleotide vector
};


