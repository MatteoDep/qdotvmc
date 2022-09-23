#ifndef __VECTOR_H
#define __VECTOR_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>

// Vector class -----------------------------------------------------------
class vector
{
    private:
    int length;			    // length of the vector
    double* r;			    // related array
    void allocate (int n);	    // allocate memory in heap
    void deallocate ();		    // free heap memory
    public:
    //basic functions
    vector ();				    // constructor - create a NULL vector - ex: vector v; 
    vector (int n);			    // constructor - create an n-dimensional vector - ex: vector v(n);
    vector (const vector& a);		    // constructor - create a vector equal to a - ex: vector v(a);
    ~vector ();				    // deconstructor
    int getlen () const;		    // get the length of the vector
    bool redim (int n);			    // change vector dimension
    double const operator() (int i) const;  // read the i-th component of the vector
    double& operator() (int i);		    // write the i-th component of the vector
    // overloading operators
    vector& operator= (const vector& a);
    vector& operator+= (const vector& a);
    vector& operator-= (const vector& a);
    friend vector sqrt (const vector& a);
    friend vector exp (const vector& a);
    friend vector round (const vector& a);
    friend vector operator+ (const vector& a, const vector& b);
    friend vector operator+ (const double k, const vector& a);
    friend vector operator+ (const vector& a, const double k);
    friend vector operator- (const vector& a, const vector& b);
    friend vector operator- (const double k, const vector& a);
    friend vector operator- (const vector& a, const double k);
    friend vector operator* (const vector& a, const vector& b);
    friend vector operator* (const double k, const vector& a);
    friend vector operator* (const vector& a, const double k);
    friend vector operator/ (const vector& a, const vector& b);
    friend vector operator/ (const double k, const vector& a);
    friend vector operator/ (const vector& a, const double k);
    friend std::ostream &operator<< (std::ostream &out, const vector& a);
    friend std::istream &operator>> (std::istream &in, vector& a);
    // other functions
    void zeros ();
    double norm () const;		    // euclidean norm - sqrt(v·v)
    double norm2 () const;		    // euclidean norm^2 - (v·v)
    double rel (const vector& a) const;	    // relative distance - ||v-a|| 
    double inner (const vector& a) const;   // inner product - v·a
};

// Vector inline functions

inline vector round (const vector& a)
{
	vector tmp(a.getlen());
	for (int i=0; i<a.getlen(); i++) tmp(i)=rint(a(i));
	return tmp;
}

inline vector sqrt (const vector& a)
{
	vector tmp(a.getlen());
	for (int i=0; i<a.getlen(); i++) tmp(i)=sqrt(a(i));
	return tmp;
}

inline vector exp (const vector& a)
{
	vector tmp(a.getlen());
	for (int i=0; i<a.getlen(); i++) tmp(i)=exp(a(i));
	return tmp;
}

inline vector& vector::operator= (const vector& a) 
{
    redim(a.getlen()); 
    for (int i=0; i<length; i++) r[i]=a(i);
    return *this;
}

inline vector& vector::operator+= (const vector& a) 
{
    for (int i=0; i<length; i++) r[i]+=a(i);
    return *this;
}

inline vector& vector::operator-= (const vector& a) 
{
    for (int i=0; i<length; i++) r[i]-=a(i);
    return *this;
}

inline vector operator+ (const vector& a, const vector& b) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)+b(i);
    return tmp;
}

inline vector operator+ (const vector& a, const double k) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)+k;
    return tmp;
}

inline vector operator+ (const double k, const vector& a) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)+k;
    return tmp;
}

inline vector operator- (const vector& a, const vector& b) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)-b(i);
    return tmp;
}

inline vector operator- (const vector& a, const double k) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)-k;
    return tmp;
}

inline vector operator- (const double k, const vector& a) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=k-a(i);
    return tmp;
}

inline vector operator* (const vector& a, const vector& b)
{
    vector tmp(a.getlen());
    for (int i=0;i<a.getlen(); i++) tmp(i)=a(i)*b(i);
    return tmp;
}

inline vector operator* (const double k, const vector& a) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=k*a(i);
    return tmp;
}

inline vector operator* (const vector& a, const double k) 
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=k*a(i);
    return tmp;
}

inline vector operator/ (const vector& a, const vector& b)
{
    vector tmp(a.getlen());
    for (int i=0;i<a.getlen(); i++) tmp(i)=a(i)/b(i);
    return tmp;
}

inline vector operator/ (const vector& a, const double k)
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=a(i)/k;
    return tmp;
}

inline vector operator/ (const double k, const vector& a)
{ 
    vector tmp(a.getlen());
    for (int i=0; i<a.getlen(); i++) tmp(i)=k/a(i);
    return tmp;
}

inline double vector::rel (const vector& a) const {vector tmp=*this-a; return tmp.norm();}

inline double vector::norm () const {vector tmp(*this); return sqrt(tmp.inner(tmp));}

inline double vector::norm2 () const {vector tmp(*this); return tmp.inner(tmp);}

inline double vector::inner (const vector& a) const 
{
    double sum=0;
    for (int i=0; i<length; i++) sum+=r[i]*a(i);
    return sum;
}

#endif
