#ifndef __MATRIX_H
#define __MATRIX_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>
#include "vector.h"

// Matrix class -----------------------------------------------------------
class matrix
{
    private:
    int row;			    // number of rows
    int col;			    // number of columns
    double** r;			    // related array
    void allocate (int nr,int nc);  // allocate memory in heap
    void deallocate ();		    // free heap memory
    public:
    //basic methods
    matrix ();					    // constructor - create a NULL matrix - ex: matrix m; 
    matrix (int nr, int nc);			    // constructor - create an (nr x nc)-dimensional matrix - ex: matrix m(nr,nc);
    matrix (const matrix& a);			    // constructor - create a matrix equal to a - ex: matrix m(a);
    ~matrix ();					    // deconstructor
    int getrow () const;			    // get the number of rows of the matrix
    int getcol () const;			    // get the number of columns of the matrix
    bool redim (int n,int m);			    // change matrix dimensions
    double const operator() (int i, int j) const;   // read the i-j component of the matrix
    double& operator() (int i, int j);		    // write the i-j component of the matrix
    // overloading operators
    matrix& operator= (const matrix& a);
    matrix& operator+= (const matrix& a);
    friend matrix outer(const vector& a, const vector& b);
    friend matrix operator+ (const matrix& a, const matrix& b);
    friend matrix operator- (const matrix& a, const matrix& b);
    friend matrix operator* (const matrix& a, const double k);
    friend matrix operator/ (const matrix& a, const double k);
    friend matrix operator+ (const double k, const matrix& a);
    friend matrix operator- (const double k, const matrix& a);
    friend matrix operator* (const double k, const matrix& a);
    friend matrix operator/ (const double k, const matrix& a);
    friend std::ostream &operator<< (std::ostream &out, const matrix& a);
    friend std::istream &operator>> (std::istream &in, matrix& a);
    // other functions
    void zeros ();
    matrix inner (const matrix& a) const;   // inner product - v·a
    vector inner (const vector& a) const;   // inner product - v·a
    matrix minor (int i, int j) const;	    // minor i-j
    double det () const;		    // determinant
    double det2 () const;		    // optimized determinant for 2x2 matrix
    double det3 () const;		    // optimized determinant for 3x3 matrix
    matrix inverse () const;		    //inverse
};

// Matrix inline functions

inline matrix& matrix::operator= (const matrix& a) 
{
    redim(a.getrow(),a.getcol()); 
    for (int i=0; i<row; i++)
    {
	for (int j=0; j<col; j++) r[i][j]=a(i,j);
    }
    return *this;
}

inline matrix& matrix::operator+= (const matrix& a) 
{
    redim(a.getrow(),a.getcol()); 
    for (int i=0; i<row; i++)
    {
	for (int j=0; j<col; j++) r[i][j]+=a(i,j);
    }
    return *this;
}

inline matrix outer(const vector& a, const vector& b)
{
    matrix tmp(a.getlen(),b.getlen());
    for (int i=0; i<a.getlen(); i++)
    {
	for (int j=0; j<b.getlen(); j++) tmp(i,j)=a(i)*b(j);
    }
    return tmp;
}

inline matrix operator+ (const matrix& a, const matrix& b) 
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=a(i,j)+b(i,j);
    }
    return tmp;
}

inline matrix operator+ (const matrix& a, const double k) 
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=a(i,j)+k;
    }
    return tmp;
}

inline matrix operator- (const matrix& a, const matrix& b) 
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=a(i,j)-b(i,j);
    }
    return tmp;
}

inline matrix operator- (const matrix& a, const double k) 
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=a(i,j)-k;
    }
    return tmp;
}

inline matrix operator* (const double k, const matrix& a) 
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=k*a(i,j);
    }
    return tmp;
}

inline matrix operator* (const matrix& a, const matrix& b)
{
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0;i<a.getrow(); i++)
    {
	for (int j=0;j<a.getcol(); j++) tmp(i,j)=a(i,j)*b(i,j);
    }
    return tmp;
}

inline matrix operator/ (const matrix& a, const double k)
{ 
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0; i<a.getrow(); i++)
    {
	for (int j=0; j<a.getcol(); j++) tmp(i,j)=a(i,j)/k;
    }
    return tmp;
}

inline matrix operator/ (const matrix& a, const matrix& b)
{
    matrix tmp(a.getrow(),a.getcol());
    for (int i=0;i<a.getrow(); i++)
    {
	for (int j=0;j<a.getcol(); j++) tmp(i,j)=a(i,j)/b(i,j);
    }
    return tmp;
}

inline matrix matrix::inner (const matrix& a) const 
{
    matrix tmp(row,a.getcol());
    tmp.zeros();
    if (a.getrow()!=col)
    {
	std::cout<<"Nothing happened: dimensions mismatch\n";
	return *this;
    }
    else 
    {
	for (int i=0; i<row; i++)
	    {
		for (int j=0; j<a.getcol(); j++) 
		{
		    for (int k=0; k<col; k++) tmp(i,j)+=r[i][k]*a(k,j);
		}
	    }
    }
    return tmp;
}

inline vector matrix::inner (const vector& a) const 
{
    vector tmp(row);
    tmp.zeros();
    if (a.getlen()!=col)
    {
	std::cout<<"Nothing happened: dimensions mismatch\nreturning zero";
	return tmp;
    }
    else 
    {
	for (int i=0; i<row; i++)
	    {
		for (int k=0; k<col; k++) tmp(i)+=r[i][k]*a(k);
	    }
    }
    return tmp;
}

inline matrix matrix::minor (int i, int j) const
{
    matrix tmp(row-1,col-1);
    for (int k=0; k<i; k++)
    {
	for (int p=0; p<j; p++) tmp(k,p)=r[k][p];
	for (int p=j+1; p<col; p++) tmp(k,p-1)=r[k][p];
    }
    for (int k=i+1; k<row; k++)
    {
	for (int p=0; p<j; p++) tmp(k-1,p)=r[k][p];
	for (int p=j+1; p<col; p++) tmp(k-1,p-1)=r[k][p];
    }
    return tmp;
}

inline double matrix::det () const
{
    double det=0;
    if (row!=col) std::cout<<"Error!: Cannot find determinant of a non-square matrix\n";
    else
    {
	switch (row)
	{
	case 1 :
	    det=r[0][0];
	    break;
	case 2 :
	    det=det2();
	    break;
	case 3 :
	    det=det3();
	    break;
	default :
	    for (int i=0; i<row; i++)
	    {
		if (i%2==0) det+=r[i][0]*(this->minor(i,0)).det();
		else det-=r[i][0]*(this->minor(i,0)).det();
	    }
	}
    }
    return det;
}

inline double matrix::det2 () const
{
    double det;
    det=r[0][0]*r[1][1]-r[0][1]*r[1][0];
    return det;
}

inline double matrix::det3 () const
{
    double det;
    det=r[0][0]*(r[1][1]*r[2][2]-r[2][1]*r[1][2])
	- r[1][0]*(r[0][1]*r[2][2]-r[2][1]*r[0][2])
	+ r[2][0]*(r[0][1]*r[1][2]-r[1][1]*r[0][2]);
    return det;
}

inline matrix matrix::inverse () const
{
    double det=this->det();
    matrix tmp(row,row);
    if (row!=col) std::cout<<"Error!: Cannot find determinant of a non-square matrix\n";
    else
    {
	if (det==0) std::cout<<"Error!: matrix is singular\n";
	else
	{
	    if (row==1) tmp(0,0)=1/r[0][0];
	    else
	    {
		for (int i=0; i<row; i++)
		{
		    for (int j=0; j<col; j++)
		    {
			if ((i+j)%2==0) tmp(j,i)=(this->minor(i,j)).det()/det;
			else tmp(j,i)=-(this->minor(i,j)).det()/det;
		    }
		}
	    }
	}
    }
    return tmp;
}

#endif
