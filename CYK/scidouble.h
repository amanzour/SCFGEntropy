//#pragma once

#ifndef __SCIDOUBLE_H__
#define __SCIDOUBLE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

class scidouble
{
public:
	scidouble(void);
	scidouble(double value);
	~scidouble(void);
	//1230: coef = 1.23, base = 3
	double coef;
	int base;
	void   AssignDouble (const double & value); // assign coef and base
	void   Double2Sci (double invalue, double & outcoef, int & outbase);// transfer double into scientific notation
	double Sci2Double (const scidouble & scivalue); // transfer scientific notation into double
	double Sci2Double (); //get the double value of the object itself
	scidouble operator + (const scidouble & scivalue);
	scidouble operator * (const scidouble & scivalue);
	scidouble operator / (const scidouble & scivalue);
	bool operator > (const scidouble & scivalue);
};
#endif


