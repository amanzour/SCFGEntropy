#include "entropy.h"

char underthreshold =' ';
char overthreshold  ='*';


entropy::entropy(int nonnum, vector<ruleinfo> thegrammar, vector< vector< int > > malpha, vector< vector< int > > mbeta, int maxWinWidth,
		  int pairnum, double threshold_p, int signal_pij, vector< vector<double> > bpmatrix, vector< double > nvector, string pij_filename)
{
	inout.setmaxwidth(maxWinWidth);
	inout.setgrammar(nonnum,thegrammar,malpha,mbeta);
	inout.betaspace(maxWinWidth);
	inout.Initilizealpha();
	thepairnum = pairnum;
	thethreshold = threshold_p;
	pij_file = pij_filename;
	pijsignal = signal_pij;
	basepair_matrix = bpmatrix;
        n_vector=nvector;
	inout.setbpmatrix(basepair_matrix);
        inout.setnvector(n_vector);
	if (pijsignal != 0)
	{
	    f_p.open(pij_file.c_str(), ios::out);
	    if (!f_p.is_open())
	    {
		cerr<<"pij file can't be open"<<endl;
		exit(1);
	    }
	}
	themalpha = malpha;
	wholegrammar = thegrammar;
}

entropy::~entropy(void)
{
}

double entropy::calentropy(int winbeg, int winwidth, const string &oristr)
//begin position of the window in whole genome, length, original sequence
{
	int i, j;//, pcount=1;
	double entr=0;
	scidouble p(0), p_whole, sum_p, sum_temp;
	vector<scidouble> p_array;
	string winstr = oristr.substr(winbeg, winwidth);
	inout.alphatable(winstr);
	inout.betatable(winstr, winbeg);// generate beta table for each window

	p_whole = inout.alpha(0, 0, winwidth-1);//Prob[s]//should be executed after betatable, otherwise winwidth has no value
	if (pijsignal != 0)
	{
	    f_p<<"alpha[1, N]: "<<p_whole.Sci2Double()<<endl;
	    f_p<<"i j pij"<<endl;
	}
	//cout<<(-1)*log(p_whole)<<endl;
	//cout<< inout.alpha(0, 6, 65)<<endl;
	//cout<< inout.alpha(0, 5, 67)<<endl;
	p_array.clear();
	sum_p.AssignDouble(0);
	for (i=0; i<winwidth; i++)
	{
		//for (j=(i+1+inout.mindist); j<winwidth; j++)
		for (j=i; j<winwidth; j++)
		{
			p = segprobability(i, j, winstr, p_whole); //P[i,j]=Prob(s[i]^s[j]|s)
			sum_p = sum_p + p;
			p_array.push_back(p);
			if (pijsignal != 0)
			    f_p<<i+1<<" "<<j+1<<" "<<p.Sci2Double()<<endl;
		}
		//cout<<endl;
	}
	sum_temp = sum_p;
	//cout << "sum_p: " << sum_p.coef << " " << sum_p.base << endl;
	//sum_p = 1; //no normalization;
	if (thethreshold > -50)
	{
		plottingPij(p_array, sum_p, winwidth);
		outputPij(p_array, sum_p, winwidth);
	}
	for (i=0; i<(int)p_array.size(); i++)
		if (p_array.at(i).coef >= 1)
		//if ((p_array.at(i)/sum_p)>0.01)
		    //entr -= (p_array.at(i)/sum_p)*log((p_array.at(i)/sum_p));
		    entr -= callogodd (p_array.at(i), sum_p);
	entr = (double)(entr);

	if (pijsignal != 0)
	   f_p.close();
	return entr;
	//return sum_temp;
}

void entropy::traceSS(string & curSS, int nonterm, int beg, int end)
{
    ruleinfo currule;
    int left_nt, right_nt, ruleid;
    left_nt  = inout.cykchoice(nonterm, beg, end, 0);
    right_nt = inout.cykchoice(nonterm, beg, end, 1);
    ruleid   = inout.cykchoice(nonterm, beg, end, 2);
    currule  = wholegrammar.at(themalpha.at(nonterm).at(ruleid));
    if ((currule.ruletype == 1)&&(beg==end)) //X:a
    {
        curSS.at(left_nt)='.';
    }
    if ((currule.ruletype == 2)&&(beg<end)) //X:aH
    {
        curSS.at(left_nt)='.';
        traceSS(curSS, currule.non1, beg+1, end);
    }
    if ((currule.ruletype == 3)&&(right_nt>=(left_nt+inout.mindist+1))) //X:aHb
    {
        curSS.at(left_nt)='(';
        curSS.at(right_nt)=')';
        traceSS(curSS, currule.non1, beg+1, end-1);
    }
    if ((currule.ruletype == 4)&&(right_nt>=(left_nt+inout.mindist+1))) //X:aHbY
    {
        curSS.at(left_nt)='(';
        curSS.at(right_nt)=')';
        traceSS(curSS, currule.non1, left_nt+1, right_nt-1);
        traceSS(curSS, currule.non2, right_nt+1, end);
    }
    if ((currule.ruletype == 5)&&(beg<end)) //X:Ha
    {
        curSS.at(right_nt)='.';
        traceSS(curSS, currule.non1, beg, end-1);
    }
    if ((currule.ruletype == 6)&&(right_nt>=(left_nt+inout.mindist+1))) //X:HaYb
    {
        curSS.at(left_nt)='(';
        curSS.at(right_nt)=')';
        traceSS(curSS, currule.non1, beg, left_nt-1);
        traceSS(curSS, currule.non2, left_nt+1, right_nt-1);
    }

}

string entropy::getSS(int winbeg, int winwidth, const string &oristr)
//begin position of the window in whole genome, length, original sequence
{
	//int i, j;
	string winstr = oristr.substr(winbeg, winwidth);
	string strSS=winstr; //Secondary Structure
	inout.alphatable(winstr);
	//inout.betatable(winstr, winbeg);// generate beta table for each window
	//p_whole = inout.alpha(0, 0, winwidth-1);//Prob[s]//should be executed after betatable, otherwise winwidth has no value
    traceSS(strSS, 0, 0, (int)(strSS.length()-1));



	return strSS;
}


double entropy::callogodd (scidouble p_item, scidouble p_sum)
{
	double logvalue;
	scidouble curp;
	curp = p_item / p_sum;
	logvalue = log(curp.coef) + curp.base * log(double(10));
	logvalue *= curp.Sci2Double();
	return logvalue;
}

scidouble entropy::segprobability(int beg, int end, const string &window, scidouble p_whole)//Prob(s[i]^s[j]|s)
//begin position, length, sequence within window, begin position of the window in whole genome, probability of sequence within window
{
	scidouble probability(0), sita;
	int length, k, v;
	ruleinfo currule;

	//if ((beg==44)&&(end==63))
	//	beg=beg;

	length=window.length();

	sita = inout.basepair(beg,end,window)/p_whole;

	if ((end-beg)<(inout.mindist+1))
		probability.AssignDouble(0);
	else
	{
		for (v=0; v<(int)inout.grammar.size(); v++)
		{
			currule = inout.grammar.at(v);
			if ((currule.leftnon!=0) && (beg==0) && (end==(length-1)))
				continue;
			if ((currule.ruletype!=3) && (beg==0) && (end==(length-1)))
				continue;
			if (currule.ruletype == 3) // v:av'b
				probability = probability + scidouble(currule.prob) * inout.beta(currule.leftnon,beg-1,end+1)
								* inout.alpha(currule.non1,beg+1,end-1);
			if (currule.ruletype == 4) // v:av_1bv_2
				for (k=end+1; k<length; k++)
					probability = probability +
									scidouble(currule.prob)
									* inout.beta(currule.leftnon,beg-1,k+1)
									* inout.alpha(currule.non1,beg+1,end-1)
									* inout.alpha(currule.non2,end+1,k);
			if (currule.ruletype == 6) // v:v_1av_2b
				for (k=0; k<beg; k++)
					probability = probability +
									scidouble(currule.prob)
									* inout.beta(currule.leftnon,k-1,end+1)
									* inout.alpha(currule.non1, k, beg-1)
									* inout.alpha(currule.non2, beg+1, end-1);
		}
		probability = sita*probability;
	}
	//if(probability>1)
	//	probability=probability;
	return probability;
}

/*
double entropy::alphapair(int beg, int end, const string &window, int winbeg)//alpha^(i,j)
//begin position, length, sequence within window, begin position of the window in whole genome
{
	int i;
	double prop=0, p_in;
	if ((end-beg)>=(inout.mindist+1))
	{
		//for (i=0; i<thepairnum; i++)
		p_in = inout.alpha(beg+1, end-1, winbeg); //alpha(i+1,j-1)
		prop = p_in;
		for (i=0; i<thepairnum; i++)
			prop *= inout.basepair(beg-i, end+i, window);
	}
	return prop;
}
*/
void entropy::plottingPij(vector<scidouble> & p, scidouble psum, int winwidth)
{
	int i, j, account=0;
	double sum_below=0, sum_above=0;

	f_p.width(3);
	f_p<<" "<<" ";
	for (j=0; j<winwidth; j++)
	{
		if (j%5==0)
		{
			f_p.width(3);
			f_p<<(j+1)<<" ";
			continue;
		}
			f_p.width(3);
			f_p<<"-"<<" ";
	}
	f_p<<endl;
	for (i=0; i<winwidth; i++)
	{
		f_p.width(3);
		f_p<<(i+1)<<" ";
		for (j=0; j<winwidth; j++)
		{
			if (j<i)
			{
				f_p.width(3);
				f_p<<underthreshold<<" ";
				continue;
			}
			if ((p.at(account)).Sci2Double()>=thethreshold)
			{
				f_p.width(3);
				f_p<<overthreshold<<" ";
				sum_above += (p.at(account)).Sci2Double();
			}
			else
			{
				f_p.width(3);
				f_p<<underthreshold<<" ";
				sum_below += (p.at(account)).Sci2Double();
			}
			account++;
		}
		f_p<<endl;
	}
	f_p<<"below threshold: "<<sum_below<<endl;
	f_p<<"above threshold: "<<sum_above<<endl;
}

void entropy::outputPij(vector<scidouble> & p, scidouble psum, int winwidth)
{
	int i, j, account=0;
	//cout.width(3);
        f_p<<" "<<" ";
	for (j=0; j<winwidth; j++)
	{
		if (j%5==0)
		{
			//cout.width(3);
			f_p<<(j+1)<<" ";
			continue;
		}
			//cout.width(3);
			f_p<<"-"<<" ";
	}
	f_p<<endl;
	for (i=0; i<winwidth; i++)
	{
		//cout.width(3);
		f_p<<(i+1)<<" ";
		for (j=0; j<winwidth; j++)
		{
			if (j<i)
			{
				//cout.width(3);
				f_p<<underthreshold<<" ";
				continue;
			}
			if ((p.at(account)).Sci2Double()>=thethreshold)
			{
				//cout.width(3);
				f_p<<(p.at(account)).Sci2Double()<<" ";
			}
			else
			{
				//cout.width(3);
				f_p<<underthreshold<<" ";
			}
			account++;
		}
		f_p<<endl;
	}
}


