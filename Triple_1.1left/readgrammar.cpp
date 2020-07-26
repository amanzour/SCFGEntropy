#include "readgrammar.h"

readgrammar::readgrammar(string filename)
{
	gfilename = filename;
	nonterm   = "<Start>";
	map_alpha.push_back(line);
	map_beta.push_back(line);
	grammar.resize(0);
}

readgrammar::~readgrammar(void)
{
}

void readgrammar::readtext(void)
{
	string str_line="";
	f_g.open(gfilename.c_str(), ios::in);
	if (!f_g.is_open())
	{
		cerr<<"grammar file can't be open"<<endl;
		exit(1);
	}
	while (!f_g.eof())
	{
		getline(f_g, str_line);
		if (str_line == "")
			continue;
		if (str_line[0] == '#')
			break;
		handlerule(str_line);
	}
	f_g.close();
/*	int i, j;
	cout<<"alpha size: "<<map_alpha.size()<<endl;
	for (i=0; i<map_alpha.size(); i++)
	{
		cout<<map_alpha.at(i).size()<<":";
		for (j=0; j<map_alpha.at(i).size(); j++)
			cout<<map_alpha.at(i).at(j)<<", ";
		cout<<endl;
	}
	cout<<endl;
	cout<<"beta size: "<<map_beta.size()<<endl;
	for (i=0; i<map_beta.size(); i++)
	{
		cout<<map_beta.at(i).size()<<":";
		for (j=0; j<map_beta.at(i).size(); j++)
			cout<<map_beta.at(i).at(j)<<", ";
		cout<<endl;
	}
	cout<<endl;
	cout<<"grammar size: "<<grammar.size()<<endl;
	for (i=0; i<grammar.size(); i++)
	{
		cout<<grammar.at(i).leftnon<<", ";
		cout<<grammar.at(i).non1<<", ";
		cout<<grammar.at(i).non2<<", ";
		cout<<grammar.at(i).ruletype<<", ";
		cout<<grammar.at(i).prob<<endl;
	}*/
//	cout<<map_alpha.at(1).at(1)<<endl;
//	cout<<map_alpha.at(1).at(2)<<endl;
//	cout<<map_alpha.at(0).at(3)<<endl;
	//cout<<map_beta.at(0).size()<<endl;
}

double readgrammar::getvalue(string numstr)
{
	double theValue;
	istringstream istr;
	istr.str(numstr);
	istr>>theValue;
	return theValue;
}

int readgrammar::nonterminalid(string therule, int non_beg,
		int non_end)
{
	int oriposi, curposi;
	
	string currule;
	currule = therule;

	oriposi = nonterm.find(currule.substr(non_beg,non_end-non_beg+1));
	if (oriposi == (int)string::npos)
	{
		nonterm += currule.substr(non_beg,non_end-non_beg+1);
		curposi = accout_char(nonterm, "<")-1;

		// for mapping_alpha
		map_alpha.push_back(line);
		// for mapping_beta
		map_beta.push_back(line);
	}else
		curposi = accout_char(nonterm.substr(0, oriposi+1), "<")-1;
	return curposi;
}

void readgrammar::parserule(string therule)
{
	//int firstposi; //verify whether the first symbol is terminal
	string currule;
	int colon_posi=-100, brac_posi=-100;
	//string nonlist="ABCDEFGHIJKLMNOPQRSTUVWXYZ"; // list of all possible non-terminal
	//initialization
	currule = therule;

	colon_posi = currule.find(":");

	if (colon_posi == (int)string::npos)
		errorreport(2);
	
	if (currule.at(colon_posi+1) == '<' )
	{
		brac_posi = currule.find("<", colon_posi+2);
		if (brac_posi == (int)string::npos)
			rule.ruletype = 5; //<X>:<H>a
		else
			rule.ruletype = 6; //<X>:<H>a<Y>b
	}
	else
		identify_type(currule, colon_posi);

	record_non(currule, colon_posi);

	update_map_alpha();// update map_alpha
	
	if (rule.ruletype >=2) 
		update_map_beta(); // update map_beta

}

void readgrammar::update_map_alpha(void)
{
	map_alpha.at(rule.leftnon).push_back(grammar.size());// the current rule has not been stored in grammar
}

void readgrammar::update_map_beta(void)
{
	map_beta.at(rule.non1).push_back(grammar.size()); // the current rule has not been stored in grammar

	if ((rule.ruletype == 4)&&(rule.non1 != rule.non2)) // avoid X:aHbH
		map_beta.at(rule.non2).push_back(grammar.size());
	
	if ((rule.ruletype == 6)&&(rule.non1 != rule.non2)) // avoid X:HaHb
		map_beta.at(rule.non2).push_back(grammar.size());
}

void readgrammar::handlerule(string rulestr)
{
	int posi_space;
	posi_space = rulestr.find(" ");
	if ((posi_space == (int) string::npos)||(posi_space == (int)(rulestr.length()-1)))
		errorreport(2); //rule format error
	// probability
	rule.prob = getvalue(rulestr.substr(posi_space, rulestr.length() - posi_space));

	// parse rule
	parserule(rulestr.substr(0, posi_space));

	//store in grammar
	grammar.push_back(rule);

}

int readgrammar::getnonterminalnum(void)
{
	int num;
	//we return # of chars are as same as that of non-termals.
	num = accout_char(nonterm, "<");
	//cout<<"# of nonterminal: "<<num<<endl;
	return num;
}

vector< ruleinfo > readgrammar::getgrammar(void)
{
	return grammar;
}

vector< vector< int > > readgrammar::get_alpha(void)
{
	return map_alpha;
}

vector< vector< int > > readgrammar::get_beta(void)
{
	return map_beta;
}

void readgrammar::errorreport(int errorid)
{
	if (errorid ==2)
		cerr<<"rule should be like <X>:a<H>b<Y> 0.25"<<endl;
	if (errorid>=1000)
		cerr<<"unknown error!!!"<<endl;
	exit(1);
}

void readgrammar::identify_type(string therule, int colon_posi)
//identify rule type 1, 2, 3, 4
{
	int brac_posi;
	string str_rule;
	str_rule = therule;
	if ((int)str_rule.length() == colon_posi+2)
		rule.ruletype = 1; //<X>:a
	else
		if (str_rule.at(str_rule.length()-1) != '>')
			rule.ruletype = 3; // <X>:a<H>b
		else
		{
			brac_posi = str_rule.find("<", colon_posi+3);
			if (brac_posi == (int)string::npos)
				rule.ruletype = 2; // <X>:a<H>
			else
				rule.ruletype = 4; // <X>:a<H>b<Y>
		}
}

void readgrammar::record_non(string therule, int colon_posi)
// record all non-terminals
{
	int beg_posi, end_posi;
	string currule;
	rule.leftnon = -1;
	rule.non1    = -1;
	rule.non2    = -1;
    currule = therule;

	rule.leftnon = nonterminalid(currule, 0, colon_posi-1);
	if ((rule.ruletype >= 2)&&(rule.ruletype <= 4))
	{
		beg_posi = currule.find("<", colon_posi+1);
		end_posi = currule.find(">", colon_posi+1);
		rule.non1 = nonterminalid(currule, beg_posi, end_posi);
	}
	if (rule.ruletype == 4)
	{
		beg_posi = currule.find("<", beg_posi+1);
		end_posi = currule.find(">", end_posi+1);
		rule.non2 = nonterminalid(currule, beg_posi, end_posi);
	}
	if ((rule.ruletype == 5) || (rule.ruletype == 6))
	{
		beg_posi = currule.find("<", colon_posi+1);
		end_posi = currule.find(">", colon_posi+1);
		rule.non1 = nonterminalid(currule, beg_posi, end_posi);
	}
	if (rule.ruletype == 6)
	{
		beg_posi = currule.find("<", beg_posi+1);
		end_posi = currule.find(">", end_posi+1);
		rule.non2 = nonterminalid(currule, beg_posi, end_posi);
	}
}

int readgrammar::accout_char(string target_str, string search_str)
{
	int num = 0, posi;

	posi = target_str.find(search_str);
	while (posi != (int) string::npos)
	{
		num++;
		posi = target_str.find(search_str, posi+1);
	}

	return num;
}
