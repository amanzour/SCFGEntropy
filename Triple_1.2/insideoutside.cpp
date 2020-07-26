#include "insideoutside.h"

/*
double prior_matrix[5][5]={{0,0.0004,0.0003,0.0004,0.0003},
							{0.0005,0.0052,0.0081,0.0055,0.1593},
							{0.0004,0.0089,0.004,0.2474,0.0053},
							{0.0004,0.0062,0.2503,0.0057,0.0516},
							{0.0004,0.1696,0.0055,0.0536,0.0107}};
*/
/*
double prior_matrix[5][5]=
{
 {0, 0.0000,0.0000,0.0000,0.0000},
 {0.0000,0.002, 0.002, 0.002,0.16},
 {0.0000,0.002, 0.002, 0.22, 0.002},
 {0.0000,0.002, 0.22, 0.002, 0.08},
 {0.0000,0.16,0.002,0.08,0.002}
};
*/


/*
double prior_matrix[5][5]=
{
 {0, 0.0000,0.0000,0.0000,0.0000},
 {0.0000,0.02, 0.02, 0.02,0.16},
 {0.0000,0.02, 0.02, 0.22, 0.02},
 {0.0000,0.02, 0.22, 0.02, 0.08},
 {0.0000,0.16,0.02,0.08,0.02}
};
*/


/*
double prior_matrix[5][5]=
{
 {0, 0.0000,0.0000,0.0000,0.0000},
 {0.0000,0.00001, 0.00001, 0.00001,0.17},
 {0.0000,0.00001, 0.00001, 0.25,  0.00001},
 {0.0000,0.00001, 0.25,   0.00001,0.08},
 {0.0000,0.17,0.00001,0.08,0.00001}
};
*/

/*
double prior_matrix[5][5]=
{
 {0, 0.0000,0.0000,0.0000,0.0000},
 {0.0000,0.00000, 0.00000, 0.00000,0.17},
 {0.0000,0.00000, 0.00000, 0.25,  0.00000},
 {0.0000,0.00000, 0.25,   0.00000,0.08},
 {0.0000,0.17,0.00000,0.08,0.00000}
};  // latest 10/09/2009
*/
/*
double prior_matrix[5][5]=
{
 {0, 0.0000,0.0000,0.0000,0.0000},
 {0.0000,0.00001, 0.00001, 0.00001,0.17},
 {0.0000,0.00001, 0.00001, 0.25,  0.0001},
 {0.0000,0.00001, 0.25,   0.00001,0.08},
 {0.0000,0.17,0.00001,0.08,0.00001}
};*/


insideoutside::insideoutside(void)
{
	//double ruleprob[] = {0.25,0.25,0.25,0.25};
	mindist=MINDIST;
/*
	ruleprob[0]=0.48; //X-->a
	ruleprob[1]=0.48; //X-->aX
	ruleprob[2]=0.02; //X-->aHaX
	ruleprob[3]=0.02; //X-->aHb

	ruleprob[4]=0.9;  //H-->aHb
	ruleprob[5]=0.1;  //H-->aYb
	
	ruleprob[6]=1.0;  //Y-->aZb

	ruleprob[7]=0.95;  //Z-->aX
	ruleprob[8]=0.05;  //Z-->aHbX
*/
/*
	ruleprob[0]=0.4; //X-->a
	ruleprob[1]=0.4; //X-->aX
	ruleprob[2]=0.06; //X-->aHaX
	ruleprob[3]=0.06; //X-->aHb

	ruleprob[4]=0.06;  //H-->aHb
	ruleprob[5]=0.06;  //H-->aYb
	
	ruleprob[6]=0.06;  //Y-->aZb

	ruleprob[7]=0.4;  //Z-->aX
	ruleprob[8]=0.06;  //Z-->aHbX
*/

/*
	ruleprob[0]=0.25; //X-->a
	ruleprob[1]=0.25; //X-->aX
	ruleprob[2]=0.25; //X-->aHaX
	ruleprob[3]=0.25; //X-->aHb

	ruleprob[4]=0.5;  //H-->aHb
	ruleprob[5]=0.5;  //H-->aYb
	
	ruleprob[6]=1;  //Y-->aZb

	ruleprob[7]=0.5;  //Z-->aX
	ruleprob[8]=0.5;  //Z-->aHbX
*/	
	
}

insideoutside::~insideoutside(void)
{
}

scidouble insideoutside::alpha(int v, int i, int j)
{
	scidouble reval;
	if ((i >= 0)&&(i <= j)&&(j < winwidth))
		//reval = alphamatrix.at(winbeg+i).at(j-i);
		reval = alphamatrix.at(v).at(i).at(j);
	else
		if (i>j)
			reval.AssignDouble(1);
		else
			reval.AssignDouble(0);
	return reval;
}
scidouble insideoutside::beta(int v, int i, int j)
{
	  if ((v==0) && (i < 0) && (j > winwidth - 1)) 
		  return scidouble(1);
	  if ((v>=1) && (i < 0) && (j > winwidth - 1)) 
		  return scidouble(0);
	  if (i>j)
		  return scidouble(1);
	  
	  if (((i==-1)&&(j==-1))||((i==winwidth)&&(j==winwidth)))
		  return scidouble(1);//!!!

	  if ((i==-1)&&(j<=winwidth-1))
			return betamatrix.at(v).at(winwidth).at(j);//!!!

      if (i < -1) 
		  return scidouble(0);
	  	  
	  /*
      if ((v>=1)&&(j > winwidth - 1))
		    return 0;
	  if (j>winwidth)
		  return 0;
	  */
	  if (j>winwidth-1)
		  return betamatrix.at(v).at(i).at(winwidth);

	  return betamatrix.at(v).at(i).at(j);
}

void insideoutside::alphatable(const string &winstr)
{
	int i, winlength;
	vector< vector<scidouble> > alphat;
	//clearalpha();
	winlength = winstr.length();
	if (flagalpha == 0)
	{
		value.resize(maxwidth);
		for (i=0; i<winlength; i++)
			alphat.push_back(value);

		for (i=0; i< nontnum; i++)
			alphamatrix.push_back(alphat);
		
		generatealpha(winstr, winlength);
		flagalpha = 1;
	} else
		generatealpha(winstr, winlength);
		//slidealpha(winstr, winlength);
}

void insideoutside::cal_alpha(int v, int i, int j, const string &oristr)
{
	scidouble curalpha(0);
	int r, k;
	ruleinfo currule;
	for (r=0; r<(int)map_alpha.at(v).size(); r++)
	{
		currule = grammar.at(map_alpha.at(v).at(r));
		
		if ((currule.ruletype == 1)&&(i==j)) //X:a
			curalpha = curalpha + scidouble(currule.prob) * qprob(i, oristr);

		if ((currule.ruletype == 2)&&(i<j)) //X:aH
			curalpha = curalpha + scidouble(currule.prob) * qprob(i, oristr) * alpha(currule.non1, i+1, j);

		if ((currule.ruletype == 3)&&(j>=(i+mindist+1))) //X:aHb
			curalpha = curalpha + scidouble(currule.prob) * basepair(i,j,oristr) * alpha(currule.non1, i+1, j-1);

		if ((currule.ruletype == 4)&&(j>=(i+mindist+1))) //X:aHbY
			for (k=i+1; k<=j-1; k++)
				curalpha = curalpha + scidouble(currule.prob) 
							* basepair(i,k,oristr) 
							* alpha(currule.non1, i+1, k-1)
							* alpha(currule.non2, k+1, j);

		if ((currule.ruletype == 5)&&(i<j)) //X:Ha
			curalpha = curalpha + scidouble(currule.prob) * qprob(j, oristr) * alpha(currule.non1, i, j-1);
		
		if ((currule.ruletype == 6)&&(j>=(i+mindist+1))) //X:HaYb
			for (k=i+1; k<=j-1; k++)
				curalpha = curalpha + scidouble(currule.prob)
							* basepair(k, j, oristr)
							* alpha(currule.non1, i, k-1)
							* alpha(currule.non2, k+1, j-1);

	}
	alphamatrix.at(v).at(i).at(j) = curalpha;
}

void insideoutside::generatealpha(const string &oristr, int winlength)
{
	int i, j, v, incre;
	for (incre = 0; incre < winlength; incre++)
		for (i = 0; i<winlength-incre; i++)
		{
			j= i+incre;
			for (v = 0; v < nontnum; v++)
				cal_alpha (v, i, j, oristr);
		}
}

void insideoutside::cal_beta(int v, int i, int j, const string &winstr)
{
	scidouble curbeta(0);
	int r, k;
	ruleinfo currule;

	for (r=0; r<(int)map_beta.at(v).size(); r++)
	{
		currule = grammar.at(map_beta.at(v).at(r));

		if (j == winwidth)
		{
			if (currule.ruletype == 2) //v':av 
				curbeta = curbeta + scidouble(currule.prob) * qprob(i, winstr) * beta(currule.leftnon,i-1,winwidth);
			if ((currule.ruletype == 4) && (currule.non2 == v)) //v':av''bv
				for (k=0; k<=(i-mindist-1); k++)
					curbeta = curbeta + scidouble(currule.prob)*basepair(k,i,winstr)*beta(currule.leftnon,k-1,winwidth)
								*alpha(currule.non1,k+1,i-1);
		}
	
		if (j <winwidth)
		{
			if (currule.ruletype == 2) // v':av
				curbeta = curbeta + scidouble(currule.prob) * qprob(i, winstr) * beta(currule.leftnon, i-1, j);

			if ((currule.ruletype == 3) && (j>=(i+mindist+1))) // v':avb
				curbeta = curbeta + scidouble(currule.prob) * basepair(i,j,winstr) * beta(currule.leftnon, i-1, j+1);

			if ((currule.ruletype == 4)) //v':av_1bv_2
			{
				if (currule.non1 == v) // v_1 == v
					for (k=j+1; k<winwidth; k++)
						curbeta = curbeta + scidouble(currule.prob) * basepair(i,j,winstr)
									* beta(currule.leftnon,i-1,k+1)
									* alpha(currule.non2,j+1,k);

				if (currule.non2 == v) // v_2 == v
					for (k=0; k<=(i-mindist-1); k++)
						curbeta = curbeta + scidouble(currule.prob)*basepair(k,i,winstr)
									*beta(currule.leftnon,k-1,j)
									*alpha(currule.non1,k+1,i-1); 
			}
			if (currule.ruletype == 5) // v':va
				curbeta = curbeta + scidouble(currule.prob) * qprob(j, winstr) * beta(currule.leftnon,i, j+1);

			if (currule.ruletype == 6) //v':v_1av_2b
			{
				if (currule.non1 == v) // v_1 == v
					for (k=j+mindist+1; k<winwidth; k++)
						curbeta = curbeta + scidouble(currule.prob)*basepair(j,k,winstr)
									* beta(currule.leftnon,i, k+1)
									* alpha(currule.non2, j+1, k-1);
				if (currule.non2 == v) // v_2 == v
					for (k=0; k<=(i-1); k++)
						curbeta = curbeta + scidouble(currule.prob)*basepair(i,j,winstr)
								* beta(currule.leftnon,k-1,j+1)
								* alpha(currule.non1, k, i-1);
			}
		}
	}

	betamatrix.at(v).at(i).at(j) = curbeta;
}

void insideoutside::cal_beta_leftmost(int v, int j, const string &winstr)
{
	scidouble curbeta(0);
	int r, k;
	ruleinfo currule;
	for (r=0; r<(int)map_beta.at(v).size(); r++)
	{
		currule = grammar.at(map_beta.at(v).at(r));
		if (currule.ruletype == 5) //v':va
			curbeta = curbeta + scidouble(currule.prob) * qprob(j,winstr) * beta(currule.leftnon,-1,j+1);
		if ((currule.ruletype == 6) && (currule.non1 == v)) //v':vav"b
			for (k=j+mindist+1; k<winwidth; k++)
				curbeta = curbeta + scidouble(currule.prob)*basepair(j,k,winstr)
							*beta(currule.leftnon,-1,k+1) * alpha(currule.non2,j+1,k-1);
	}
	betamatrix.at(v).at(winwidth).at(j) = curbeta;
}

void insideoutside::generatebeta(const string &winstr)
{
	int i, j, v;
	for (j=winwidth-1; j>=0; j--)
		for (v = 0; v < nontnum; v++)
			cal_beta_leftmost(v, j, winstr); //calculate beta(v, -1, j);

	for (i=0; i<winwidth; i++)
		//for (j=i; j<=winwidth; j++)
		for (j=winwidth; j>=i; j--)
			for (v = 0; v < nontnum; v++)
				cal_beta (v, i, j, winstr);
}


void insideoutside::betatable(const string &winstr, int winbeg)// winbeg starts with 0 instead of 1
{
//	int i;
	windowbegin = winbeg;
	winwidth = (int)winstr.length();
	clearbeta();
//	value.resize((winwidth+1));
//	for (i=0; i<winwidth; i++)
//		betamatrix.push_back(value);
	generatebeta(winstr);
}

void insideoutside::betaspace(int MaxWidth)
{
	int i;
	vector< vector<scidouble> > betat;
	winwidth = MaxWidth;
	value.resize((winwidth+1));
	for (i=0; i<=winwidth; i++)
		betat.push_back(value);

	for (i=0; i<nontnum; i++)
		betamatrix.push_back(betat);
}

void insideoutside::setmaxwidth(int MaxWidth)
{
	maxwidth = MaxWidth;
}

void insideoutside::clearalpha()
{
	//alphamatrix.clear();
	//value.clear();
}

void insideoutside::clearbeta()
{
	//value.clear();
	//betamatrix.clear();
}

scidouble insideoutside::qprob(int posi, const string &qstr)
// A, C, G, U
{
	double probval;
	char Nt;

	if ((posi<0)||(posi>=(int)qstr.length()))
		return scidouble(0);
	Nt=qstr[posi]; 
    switch (Nt)
	{
		case 'A': 
			probval = prior_vector[0];//0.253088577; //0.25;
			break;
		case 'C':
			probval = prior_vector[1]; //0.230557837; //0.25;
			break;
		case 'G':
			probval = prior_vector[2]; //0.258991666; //0.25;
			break;
		case 'U':
			probval = prior_vector[3]; //0.25736192; //0.25;
			break;
		default:
			probval = 0;
			break;
	}
	return scidouble(probval);
}

int insideoutside::maxvalue(int valA, int valB)
{
	int reval;
	if (valA>valB)
		reval=valA;
	else
		reval=valB;
	return reval;
}

int insideoutside::minvalue(int valA, int valB)
{
	int reval;
	if (valA<valB)
		reval=valA;
	else
		reval=valB;
	return reval;
}

double insideoutside::ppair(char ntleft, char ntright) 
// left character of the base pair, right character of the base pair,
{
	int index_left, index_right;
	double prop=0;

	index_left  = nt2num(ntleft);
	index_right = nt2num(ntright);
	if ((index_left>=0)||(index_right>=0))
		prop=prior_matrix[index_left][index_right];
	else
		prop=0;

	return prop;
}

void insideoutside::Initilizealpha(void)
{
	flagalpha= 0;
}

scidouble insideoutside::basepair(int i, int j, const string & str)
{
	double reval;
	if ((i<0)||(j>=(int)str.length()))
		reval = 0;
	else
		if (i==j)
			reval =0;//qprob(str[i]);
		else
			if (j-i<(mindist+1))
				reval = 0;
			else
				reval=ppair(str[i], str[j]);
	return scidouble(reval);
}

int insideoutside::nt2num(char nt)
{
	int value;
	switch (nt)
	{
		case 'A': value =1; break;
		case 'C': value =2; break;
		case 'G': value =3; break;
		case 'U': value =4; break;
		default : value =0;
	}
	return value;
}

void insideoutside::setgrammar(int nonnum, vector<ruleinfo> thegrammar,
							   vector< vector< int > > malpha, vector< vector< int > > mbeta)
{
	nontnum = nonnum;
	grammar = thegrammar;
	map_alpha = malpha;
	map_beta  = mbeta;
}

void insideoutside::setbpmatrix(vector< vector<double> > curbpmatrix)
{
    int i, j;
    for (i=0; i<5; i++)
    {
	for (j=0; j<5; j++)
	{
	    prior_matrix[i][j] = curbpmatrix.at(i).at(j);
	    //cout<<prior_matrix[i][j]<<" ";
	}
	//cout<<endl;
    }
}
void insideoutside::setnvector(vector< double > curvector)
{
    int i;
    for (i=0; i<4; i++)
         prior_vector[i] = curvector.at(i);
	    //cout<<prior_matrix[i][j]<<" ";
	//cout<<endl;
}

