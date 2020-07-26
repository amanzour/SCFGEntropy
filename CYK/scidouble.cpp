#include "scidouble.h"

scidouble::scidouble(void)
{
}

scidouble::~scidouble(void)
{
}

scidouble::scidouble(double value)
{
	AssignDouble(value);
}

void scidouble::AssignDouble (const double & value)
{
	Double2Sci (value, coef, base);
}

void scidouble::Double2Sci (double invalue, double & outcoef, int & outbase)
{
	outcoef = invalue;
	outbase = 0;

	while ((fabs(outcoef) < 1)&&(outcoef != 0))
	{
		outcoef *= 10;
		outbase--;
	}

	while (fabs(outcoef) >= 10)
	{
		outcoef /= 10;
		outbase++;
	}
}

double scidouble::Sci2Double(const scidouble & scivalue)
{
	double dvalue;
	dvalue = scivalue.coef*pow((double)10,scivalue.base);
	return dvalue;
}

double scidouble::Sci2Double()
{
	double dvalue;
	dvalue = coef*pow((double)10,base);
	return dvalue;
}

scidouble scidouble::operator * (const scidouble & scivalue)
{
	double realcoef;
	int    realbase;
	scidouble cursci;
	//cursci.coef = this->coef * scivalue.coef;
	//cursci.base = this->base + scivalue.base;
	Double2Sci (this->coef * scivalue.coef, realcoef, realbase);
	cursci.coef = realcoef;
	if (realcoef >= 1)
		cursci.base = this->base + scivalue.base + realbase;
	else
		cursci.base = 0;
	return cursci;
}

scidouble scidouble::operator / (const scidouble & scivalue)
{
	double realcoef;
	int    realbase;
	scidouble cursci;
	if (scivalue.coef >= 1)
	{
		Double2Sci ( this->coef / scivalue.coef, realcoef, realbase);
		cursci.coef = realcoef;
		if (realcoef >= 1)
			cursci.base = this->base - scivalue.base + realbase;
		else
			cursci.base = 0;
	}
	else
	{
		cursci.coef = 0;
		cursci.base = 0;
	}
	return cursci;
}

scidouble scidouble::operator + (const scidouble & scivalue)
{
	double realcoef;
	int    realbase;
	scidouble cursci, smallsci;
	if ((this->coef < 1 ) || (scivalue.coef < 1  ))
	{
		if (this->coef < 1 )
			cursci = scivalue;
		else
			cursci = *this;
	}
	else
	{
		if (this->base > scivalue.base)
		{
			smallsci.coef = scivalue.coef;
			smallsci.base = scivalue.base - this->base;
			Double2Sci (this->coef + Sci2Double(smallsci), realcoef, realbase);
			cursci.coef = realcoef;
			if (realcoef >= 1)
				cursci.base = this->base + realbase;
			else
				cursci.base = 0;
		}
		else
		{
			smallsci.coef = this->coef;
			smallsci.base = this->base - scivalue.base;
			Double2Sci (scivalue.coef + Sci2Double(smallsci), realcoef, realbase);
			cursci.coef = realcoef;
			if (realcoef >= 1)
				cursci.base = scivalue.base + realbase;
			else
				cursci.base = 0;
		}
	}
	return cursci;
}

bool scidouble::operator > (const scidouble & scivalue)
{
    bool reval;
    if ((fabs(this->base)<10) && (fabs(scivalue.base)<10))
        reval = (Sci2Double() > Sci2Double(scivalue)) ;
        else
            if ((this->coef)*(scivalue.coef) <= 0)
            {
                if (this->coef > scivalue.coef)
                    reval = true;
                else
                    reval = false;
            }
            else
                if (this->base == scivalue.base)
                    if (this->coef > scivalue.coef)
                        reval = true;
                    else
                        reval = false;
                else
                    if ((this->base > scivalue.base)&&(this->coef > 0))
                        reval = true;
                    else
                        if ((this->base > scivalue.base)&&(this->coef < 0))
                            reval = false;
                        else
                            if ((this->base < scivalue.base)&&(this->coef > 0))
                                reval = false;
                            else
                                if ((this->base < scivalue.base)&&(this->coef < 0))
                                    reval = true;
                                else
                                    cout<<"wrong case"<<endl;
    return reval;
}
