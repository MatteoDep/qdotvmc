#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>
#include "matrix.h"
#include "vector.h"

// Matrix functions

matrix::matrix () {r=NULL; row=0; col=0;}

matrix::matrix (int nr, int nc) {allocate(nr,nc);}

void matrix::allocate (int nr, int nc)
{
	row=nr; 
	col=nc;
	r=new double*[row];
	for (int i=0; i<row; i++) r[i]=new double[col];
}

matrix::~matrix () {deallocate();}

void matrix::deallocate ()
{
	for (int i=0; i<row; i++) delete[] r[i];
	delete[] r;
}

matrix::matrix (const matrix& a) {allocate(a.getrow(),a.getcol()); *this = a;}

std::ostream &operator<< (std::ostream& out, const matrix& a) 
{
    if (a.getrow()>0 && a.getcol()>0)
    {
	for (int i=0; i<a.getrow(); i++)
	{
	    for (int j=0; j<a.getcol()-1; j++) out<<a(i,j)<<", ";
	    out<<a(i,a.getcol()-1)<<"\n";
	}
    }
    else out<<"Err: empty matrix!";
    return out;
}

std::istream &operator>> (std::istream& in, matrix& a) 
{
	for (int i=0; i<a.getrow(); i++)
	{
		for (int j=0; j<a.getcol(); j++)
		{
			std::cout<<"element "<<i+1<<","<<j+1<<": ";
			in>>a(i,j);
		}
	}
	return in;
}

int matrix::getrow () const {return row;}

int matrix::getcol () const {return col;}

double const matrix::operator() (int i, int j) const {return r[i][j];}

double& matrix::operator() (int i, int j) {return r[i][j];}

void matrix::zeros () 
{
	for (int i=0; i<row; i++)
	{
		for (int j=0; j<col; j++) r[i][j]=0;
	}
}

bool matrix::redim (int nr, int nc) 
{
	if (row==nr && col==nc) return false;
	else 
	{
		if (r!=NULL)
		{
			deallocate();
		}
		allocate(nr,nc);
		return true;
	}
}


