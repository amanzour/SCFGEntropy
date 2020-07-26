#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <unistd.h>
//#include <getopt.h>
#include "entropy.h"
#include "scanner.h"
#include "fasta.h"
#include "readgrammar.h"
//#include "insideoutside.h"

//#define WINWIDTH 20

using namespace std;


void transfer(string &str) // change 'T' into 'G'
// string including the whole genome
{
	int i, length;
	length = str.length();
	for (i=0; i<length; i++)
		if (str.at(i) == 'T')
			str.at(i) = 'U';
}

vector< vector<double> > readbpmatrix(string & str_bpmatrix_file)
//str_bpmatrix_file is the name of base pair matrix file
{
      vector< vector<double> > matrix( 5, vector<double>(0) );
      //istringstream istr[5][5];
      ifstream f_bm;
      int i, j;
      double cur_val;
      f_bm.open(str_bpmatrix_file.c_str(), ios::in);
      if (!f_bm.is_open())
      {
	  cerr<<"base pair matrix file can't be open"<<endl;
	  exit(1);
      }
      for (i=0; i<5; i++)
	  for (j=0; j<5; j++)
	  {
	      f_bm >> cur_val;
	      matrix.at(i).push_back(cur_val);
	  }
      f_bm.close();
      return matrix;
}
vector< double > readnvector(string & str_nvector_file)
//str_bpmatrix_file is the name of base pair matrix file
{
      vector< double > vector( vector< double >(0) );
      ifstream f_bm;
      int i;
      double cur_val;
      f_bm.open(str_nvector_file.c_str(), ios::in);
      if (!f_bm.is_open())
      {
	  cerr<<"nucleotide file can't be open"<<endl;
	  exit(1);
      }
      for (i=0; i<4; i++)
	  {
	      f_bm >> cur_val;
	      vector.push_back(cur_val);
	  }
      f_bm.close();
      return vector;
}

int main (int argc, char * argv[])
{
	int i, winwidth, pairnum=1;//pairnum: the number of base pairs between alpha and beta
	double threshold=-100;
	istringstream istr,istr4;
	string str_struct, seq_grammar, seqstr="AAAUUUUUUUUAAAAAAAAAUUUUUUUUAAAAAAAAAUUUUUUUU";
	int optchr, seq_index, signal_winsize=0, signal_testallseq=0, signal_pij=0;
	string seq_filename, pij_filename, matrix_filename,vector_filename;
	vector< vector<double> > bpmatrix;
        vector< double > nvector;
//	int temp_nu_i, temp_nu_j; // for the loop of output bp matrix
	//winwidth = WINWIDTH;


	while ((optchr = getopt(argc, argv, "s:g:w:b:n:p:t:ah")) != -1)
	{
	    switch (optchr)
	    {
		case 'h':
		    cout<<"-s seqfile_name (required)"<<endl;
		    cout<<"-g grammarfile_name (required)"<<endl;
                    cout<<"-w window_size (default window size is the length of the whole sequence)"<<endl;
                    cout<<"-b base_pairmatrix_file_name (required)"<<endl;
                    cout<<"-n nucleotide_file_name (required)"<<endl;
		    cout<<"-p output_file_for_storing_all_pij"<<endl;
		    cout<<"-a indicate test all sequences in the file (without it, program just tests the first sequence)"<<endl;
		    cout<<"-t threshold  set the threshold for outputting pij matrix (-p should be specified  if -t is specified)"<<endl;
		    cout<<"-h get help"<<endl;
		    exit(1);
		    break;
		case 's':
		    seq_filename = string (optarg);
		    break;
		case 'g':
		    seq_grammar = string(optarg);
                    break;
		case 'w':
		    istr.str(string(optarg));
		    istr>>winwidth;
                    signal_winsize=1;
		    break;
		case 'b':
		    matrix_filename = string(optarg);
                    break;
                case 'n':
		    vector_filename = string(optarg);
                    break;
		case 'p':
		    signal_pij = 1;
		    pij_filename = string(optarg);
		    break;
		case 'a':
		    signal_testallseq=1;
		    break;
		case 't':
		    istr4.str(string(optarg));
		    istr4>>threshold;
		    break;
		case '?':
		    cerr<<"wrong parameter, please execute ./ncRNA -h"<<endl;
		    break;
	    }
	}

        Scanner scanner(seq_filename.c_str());
	vector<fasta> fasta_vector = scanner.getFasta();

	//cout<<"original length: "<<fasta_vector.at(0).getSequence().length()<<endl;

	readgrammar readGram(seq_grammar);

	readGram.readtext();

	bpmatrix = readbpmatrix(matrix_filename);
        nvector = readnvector(vector_filename);
/* // test of reading base pair matrix
	for (temp_nu_i=0; temp_nu_i<5; temp_nu_i++)
	{
	    for (temp_nu_j=0; temp_nu_j<5; temp_nu_j++)
		cout<<bpmatrix.at(temp_nu_i).at(temp_nu_j)<<" ";
	    cout<<endl;
	}
*/
	if (signal_testallseq == 0)
	{
	    seqstr=fasta_vector.at(0).getSequence();
	    //cout<<"seqstr length: "<<seqstr.length()<<endl;
	    transfer(seqstr);
	    //cout<<"length: "<<seqstr.length()<<endl;
	    if (signal_winsize == 0)
		winwidth = seqstr.length();

	    entropy entro(readGram.getnonterminalnum(), readGram.getgrammar(), readGram.get_alpha(),
					readGram.get_beta(),  winwidth, pairnum, threshold,
					signal_pij, bpmatrix,nvector, pij_filename);
	//entro.inout.alphatable(seqstr); // generate alpha table for the genome
        cout<<">"<<fasta_vector.at(0).getTitle()<<endl;
        cout<<seqstr<<endl;
	    for (i=0; i<=(int)(seqstr.length()-winwidth); i++)
	//for (i=0; i<=10; i++)
	    {
            str_struct = entro.getSS(i, winwidth, seqstr);
            cout<<str_struct<<endl;
	    }
	} else
	    for (seq_index=0; seq_index<(int)fasta_vector.size(); seq_index++)
	    {
	    	seqstr=fasta_vector.at(seq_index).getSequence();
	    	//cout<<"seqstr length: "<<seqstr.length()<<endl;
	    	transfer(seqstr);
	    	//cout<<"length: "<<seqstr.length()<<endl;
	    	if (signal_winsize == 0)
	    		winwidth = seqstr.length();
	    	//cout<<"seq "<<seq_index<<endl;
	    	entropy entro(readGram.getnonterminalnum(), readGram.getgrammar(), readGram.get_alpha(),
	    			readGram.get_beta(), winwidth, pairnum, threshold,
	    			signal_pij, bpmatrix,nvector, pij_filename);
	    	//entro.inout.alphatable(seqstr); // generate alpha table for the genome

                cout<<">"<<fasta_vector.at(seq_index).getTitle()<<endl;
	    	cout<<seqstr<<endl;
	    	for (i=0; i<=(int)(seqstr.length()-winwidth); i++)
	    		//for (i=0; i<=10; i++)
	    	{

	    		//temp = entro.calentropy(i, winwidth, seqstr);
	    		//temp = temp/(2*log((double)winwidth));
	    		//cout<<temp<<endl; //for outputing entropy
	    		str_struct = entro.getSS(i, winwidth, seqstr);
	    		cout<<str_struct<<endl;
	    	}

	    }
	return 0;
}



