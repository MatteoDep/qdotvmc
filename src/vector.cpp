#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>
#include "vector.h"

// Vector functions

vector::vector () {r=NULL; length=0;}

vector::vector (int n) {allocate(n);}

void vector::allocate (int n) {length=n; r=new double[n];}

vector::~vector () {deallocate();}

void vector::deallocate () {delete[] r;}

vector::vector (const vector& a) {allocate(a.getlen()); *this = a;}

std::ostream &operator<< (std::ostream& out, const vector& a) 
{
	if (a.getlen()>0)
	{
		for (int i=0; i<a.getlen()-1; i++)
		{
			out<<a(i)<<", ";
		}
		out<<a(a.getlen()-1);
	}
	else out<<"Err: empty vector!";
	return out;
}

std::istream &operator>> (std::istream& in, vector& a) 
{
	for (int i=0; i<a.getlen(); i++)
	{
		std::cout<<"element "<<i+1<<": ";
		in>>a(i);
	}
	return in;
}

int vector::getlen () const {return length;}

double const vector::operator() (int i) const {return r[i];}

double& vector::operator() (int i) {return r[i];}

void vector::zeros () 
{
	for (int i=0; i<length; i++) r[i]=0;
}

bool vector::redim (int n) 
{
	if (length==n) return false;
	else 
	{
		if (r!=NULL)
		{
			deallocate();
		}
		allocate(n);
		return true;
	}
}
