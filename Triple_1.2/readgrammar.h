//#pragma once
#ifndef __READGRAMMAR_H__
#define __READGRAMMAR_H__


#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>

using namespace std;

struct ruleinfo
{
   int    ruletype; // rule type 1: X:a 2: X:aH 3: X:aHb 4: X:aHbY 5: X:Ha  6: X:HaYb  
   int    leftnon;  // non-terminal in the left side
   int    non1;     // first non-terminal in the right side
   int    non2;     // second non-terminal in the right side
   double prob; // prob of the rule    ��
};

class readgrammar
{
public:
	readgrammar(string filename);
	~readgrammar(void);
	void readtext(void);
	int getnonterminalnum(void);
	vector< ruleinfo > getgrammar(void);
	vector< vector< int > > get_alpha(void);
	vector< vector< int > > get_beta(void);
private:
	vector< int > line; // for initializing map_alpha and map_beta
	 
	int nonterminalid(string therule, int non_beg, int non_end); // get the id of the given non-terminal
	void update_map_alpha(void); //update mapping_alpha;
	void update_map_beta(void); //update mapping_beta;
	double getvalue(string numstr); // input the string, return the value of the string like "123"->123
	void parserule(string therule); //parserule
	void handlerule(string rulestr); //deal with rules
	void errorreport(int errorid); //report error and exit
	vector< vector< int > > map_alpha; // mapping_alpha;
	vector< vector< int > > map_beta; // mapping_beta;
	vector<ruleinfo> grammar; // list of rules
	ruleinfo rule; // a rule	
	string gfilename, nonterm; //gfilename is the name of file of grammar; nonterm is for stroring all non-terminals.
	ifstream f_g; //file for grammar
	void identify_type(string therule, int colon_posi); // identify rule type
	void record_non(string therule, int colon_posi); // record all non-terminals
	int accout_char(string target_str, string search_str); //accout the # of times for search string appearing in the target string
};
#endif

