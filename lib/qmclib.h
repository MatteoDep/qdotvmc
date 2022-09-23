#ifndef __QMCLIB_H
#define __QMCLIB_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <complex>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>
#include"vector.h"
#include"matrix.h"

class Parameters;
class System;
class Slater;
class Jastrow;
class Walker;

// Parameters class -----------------------------------------------------------
class Parameters
{
    private:
    int n_conf;	    //number of set of parameters
    int count, count1D;
    bool next;
    vector* p;
    vector p_old;
    vector s;
    vector delta_p;
    vector grad_new;    // gradient in the parameters space
    vector grad_old;
    vector delta_grad;
    matrix H;
    double threshold;	// step for moving in the parameters space
    double grad_step;	// step for gradient computation
    double alpha;
    double der, der2;
    double E_0, der_0;
    void allocate ();
    void deallocate ();
    public:
    Parameters (int n);	// constructor
    void init_dfp (double guess_a, double guess_b, double m_step, double threshold_);    // init for dfp
    void init (double a_min, double a_max, double bmin, double b_max);
    ~Parameters ();
    void init_jastrow (Jastrow* jst, System& sys);  // initialiaze 2 jastrow objects with the correct parameters
    bool move (vector E);   // compute gradient and move parameters proportionally
    void iterate ();
    void updateH ();
    double calcder (vector E);
    double calcder2 (vector E);
    double geta () const;
    double getb () const;
    friend std::ostream &operator<< (std::ostream &out, const Parameters& par);
};

// System class -----------------------------------------------------------
class System
{
    private:
    int N;		    // number of particle of the system
    matrix rel;	    // Stores particles relative distances
    vector r2;	    // Stores particles distances from origin
    void allocate ();  // allocate memory in heap
    void deallocate ();	    // free heap memory
    public:
    vector* p;		    // Stores the particle positions
    // functions
    System (int n);				// constructor - create a system of N particles
    void init (double a);	    // initialize positions and distances
    void update (Walker& wlk);
    double psi (int* qn, int cn);    // calculate single particle wave function using quantum numbers qn and coordinates number cn
    vector gradpsi (int* qn, int cn);    // calculate gradient wave function using quantum numbers qn and coordinates number cn
    double lappsi (int* qn, int cn);    // calculate laplacian of the wave function using quantum numbers qn and coordinates number cn
    double getr2(int i);    // get r(i)
    double getrel(int i, int j);    // get rel(i,j)
    ~System ();					// deconstructor
};

// Slater class -----------------------------------------------------------
class Slater
{
    private:
    int N;		    // number of particle of the system
    int n_up;		    // Number of particle with spin up (the first n_up particles have spin up)
    double ratio;	// Stores determinant ratio
    int** qn;	    // Quantum numbers
    matrix Dup;	    // Stores slater determinant for spin up particles
    matrix Ddown;	    // Stores slater determinant for spin down particles
    vector** gradDup;	    // Stores slater determinant for spin up particles
    vector** gradDdown;	    // Stores slater determinant for spin down particles
    matrix lapDup;	    // Stores slater determinant for spin up particles
    matrix lapDdown;	    // Stores slater determinant for spin down particles
    matrix Dup_inv;    // Stores inverse slater determinant for spin up particles
    matrix Ddown_inv;  // Stores inverse slater determinant for spin down particles
    void allocate ();  // allocate memory in heap
    void deallocate ();	    // free heap memory
    public:
    // functions
    Slater (int n);			// constructor - create a system of N particles
    void init (System& sys);	// initialize slater determinants
    int getn_up ();		// read n_up
    double getratioup (Walker& wlk);		// get ratio det_new/det_curr
    double getratiodown (Walker& wlk);	// get ratio det_new/det_curr
    vector getgradup (int i);	// get ratio (nabla det)/det
    vector getgraddown (int i);	// get ratio (nabla det)/det
    double getlapup (int i);	// get ratio (nabla^2 det)/det
    double getlapdown (int i);	// get ratio (nabla^2 det)/det
    void updateup (Walker& wlk, System& sys);    // update matrices
    void updatedown (Walker& wlk, System& sys);    // update matrices
    ~Slater ();					// deconstructor
};

// Jastrow class -----------------------------------------------------------
class Jastrow
{
    private:
    int N;		    // number of particle of the system
    double a,b;		    // Variational parameters
    vector grad;	// Store the grad ratio
    matrix G;	    // Stores exponents of Jastrow factor g(r_ij)
    public:
    //basic functions
    Jastrow ();			// constructor - create a system of N particles
    void allocate (int n);  // allocate memory in heap
    double getG (int i, int j) const;	    // read G(i,j)
    double getratio (Walker& wlk, System& sys);	// get ratio J_new/J_curr
    vector getgrad (int i, System& sys);	// get ratio (nabla J)/J
    double getlap (int i, System& sys);	// get ratio (nabla^2 J)/J
    void update (Walker& wlk);    // update matrices
    void init (System& sys, double A, double B);	    // initialize matrices of exponents
    double g(double r_ij);		//calculate exponent of Jastrow factor
    ~Jastrow ();					// deconstructor
};

// Walker class -----------------------------------------------------------
class Walker
{
    private:
    int N;		    // number of particle of the system
    int curr;		// current particle
    double delta;	// metropolis box side
    double r2;	    // walker distance from origin
    void allocate ();  // allocate memory in heap
    void deallocate ();	    // free heap memory
    public:
    vector rel;
    vector p;	    // walker coordinates
    vector psi;	    // Stores walker wave functions
    //basic functions
    Walker (int n);			// constructor - create a system of N particles
    ~Walker ();					// deconstructor
    int current ();	// read curr
    double getr2 ();	// read r
    void calcpsi (int i, int* qn);	// calculate wave function with quantum numbers qn
    void setdelta (double d);	    // set delta
    void attempt(int i, System& sys);	// move attempt for the ith particle
    double univar ();	// uniform distributed stochastic variable
};

// Inline functions

inline double H(int n, double x)
{
	double H_=0;
	switch(n)
	{
		case 0: 
		    H_=1;
		    break;
		case 1:
		    H_=2*x;
		    break;
	}
	return H_;
}

#endif
