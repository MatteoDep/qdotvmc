/* ------------------------------------------------------ */
/* VARIATIONAL MONTE CARLO STUDY OF CIRCULAR QUANTUM DOTS */
/* ------------------------------------------------------ */

#include <cstdlib>
#include <chrono>
#include <iostream>
#include <sstream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string>
#include <fstream>
#include <time.h>
#include <complex>
// user defined headers 
#include "vector.h"
#include "matrix.h"
#include "qmclib.h"

// Reduced units ----------------------------------------------------------
#define m_star 6.7E-2		// [me]
#define a_star 97.93		// [Å]
#define H_star 11.86 		// [meV]

// Global variables --------------------------------------------
int N;
int n_conf;	// n_conf=5 for DFP
int Ref;

// Functions declaration ---------------------------------------
// Metropolis algorithm
int metropolis (System& sys, Slater& slt, Jastrow* jst, Walker& wlk);
// Calculate the total energy
vector calc_local_energy (System& sys, Slater& slt, Jastrow* jst);
// Calculate weights for reweighing method
vector reweight (Jastrow* jst);

// main function -----------------------------------------------
int main(int argc, char* argv[])
{
    // Metropolis parameters
    double delta=3.5;
    int acc, Ntot=1000, Nblock=1000;

    // read arguments
    bool dfp=false;
    std::ostringstream f_name;
    if (argc>2)
    {
	N=atoi(argv[1]);
	if ((N!=2)&&(N!=3)&&(N!=4)&&(N!=5)&&(N!=6))
	{
	    std::cout<<"invalid N: type just \"./qdotvmc\" for info\n";
	    return 1;
	}
	if (strcmp(argv[2],"surf")==0)
	{
	    dfp=false;
	    n_conf=9;
	    Ref=n_conf/2;
	    f_name<<argv[2]<<N;
	}
	else
	{
	    if (strcmp(argv[2],"dfp")==0)
	    {
		dfp=true;
		n_conf=5;
		Ref=0;
		Ntot*=50;
		f_name<<argv[2]<<N;
	    }
	    else
	    {
		std::cout<<"invalid mode: type just \"./qdotvmc\" for info\n";
		return 1;
	    }
	}
	char c;
	std::string s;
	for (int i=3; i<argc; i++)
	{
	    c=argv[i][0];
	    s.assign(argv[i]);
	    s.assign(s.begin()+2,s.end());
	    switch (c)
	    {
		case 'f':
		    f_name.str(s);
		    break;
		case 'N':
		    Ntot=stoi(s);
		    break;
		case 'n':
		    Nblock=stoi(s);
		    break;
	    }
	}
    }
    else
    {
	std::cout<<"Please rerun with arguments as follows:\n";
	std::cout<<"\"./qdotvmc N mode {f=filename} {N=Ntot} {n=Nblock}\"\n";
	std::cout<<"where:\n\tN=number of electrons, must be between 2 and 6\n";
	std::cout<<"\tmode=surf or dfp\n";
	std::cout<<"\tOptional: file_name, Ntot, Nblock (MC_cicles=Ntot*Nblock)\n";
	return 1;
    }

    //start chronometer
    auto start=std::chrono::high_resolution_clock::now();

    // Set Outputs ----------------------------------------------
    // Results
    std::ofstream f;
    f.open("Results/"+f_name.str(), std::ofstream::out | std::ofstream::trunc);
    // plot parameters
    std::ofstream plt_par;
    plt_par.open("Results/par_"+f_name.str(), std::ofstream::out | std::ofstream::trunc);
    // cout
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout.precision(10);
    f.setf( std::ios::fixed, std:: ios::floatfield );
    f.precision(10);

    // Objects initializzation ---------------------------------
    Parameters par(n_conf);
    // par.init_grad(a0,b0,grad_step,move_step); //requires n_conf=3
    System sys(N);
    sys.init(0.5);
    Slater slt(N);
    Walker wlk(N);
    wlk.setdelta(delta);
    slt.init(sys);
    Jastrow* jst;
    jst = new Jastrow[n_conf];   // allocate in heap
    for (int i=0; i<n_conf; i++) jst[i].allocate(N);

    // cumulators and averages
    vector Eloc(n_conf), Ecum_block(n_conf), Eave_block(n_conf), Ecum(n_conf), E2cum(n_conf), Eave(n_conf), E2ave(n_conf), Eerr(n_conf), norm_block(n_conf), W(n_conf);

    if (dfp)
    {
	// Parameters ------------------------------------------
	// Jastrow parameters
	Ref=0;
	double a0, b0;
	double threshold;
	double grad_step;
	std::cout<<"guess values for parameters:\n\ta=";
	std::cin>>a0; std::cout<<"\tb="; std::cin>>b0;
	std::cout<<"dfp parameters:\n\tgrad_step=";
	std::cin>>grad_step; std::cout<<"\tthreshold="; std::cin>>threshold;
        plt_par<<"data_dfp=\""<<f_name.str()<<"\""<<std::endl;

	// Main loop -------------------------------------------
	par.init_dfp(a0,b0,grad_step,threshold);
	bool k=true;
	while (k)
	{
	    acc=0;
	    Ecum.zeros(); E2cum.zeros();
	    par.init_jastrow(jst,sys);
	    for (int i=0; i<Ntot; i++)
	    {
		Ecum_block.zeros(); norm_block.zeros();
		for (int j=0; j<Nblock; j++)
		{
		    acc+=metropolis(sys,slt,jst,wlk);
		    Eloc=calc_local_energy(sys,slt,jst);
		    W=reweight(jst);
		    Ecum_block+=Eloc*W;
		    norm_block+=W;
		}
		Eave_block=Ecum_block/norm_block;
		Ecum+=Eave_block;
		E2cum+=Eave_block*Eave_block;
	    }
	    Eave=Ecum/(double)Ntot;
	    E2ave=E2cum/(double)Ntot;
	    Eerr=sqrt((E2ave-Eave*Eave)/(double)Ntot);
	    std::cout<<"\nAcceptance = "<<acc*100.0/(double)(N*Nblock*Ntot)<<"%,\n";
	    std::cout<<"E="<<Eave(0)<<"+-"<<Eerr(0)<<",\t(a="<<par.geta()<<", b="<<par.getb()<<")\n";
	    f<<par.geta()<<"\t"<<par.getb()<<"\t"<<Eave(0)<<"\t"<<Eerr(0)<<std::endl;
	    k=par.move(Eave);
	}
    }
    else
    {
	// Parameters ------------------------------------------
	// Jastrow parameters
	int n_rep;
	double a_min=0.8, a_max=1.0, b_min=0.5, b_max=0.6;
	std::cout<<"ranges for parameters:\n";
	std::cout<<"\ta_min="; std::cin>>a_min; std::cout<<"\ta_max="; std::cin>>a_max;
	std::cout<<"\tb_min="; std::cin>>b_min; std::cout<<"\tb_max="; std::cin>>b_max;
	std::cout<<"grid density:\n\tNpoints=3*"; std::cin>>n_rep;
	int Npoints=n_rep*sqrt(n_conf);
	vector a(Npoints), b(Npoints);
	for (int i=0; i<Npoints; i++)
	{
	    a(i)=a_min + i*(a_max-a_min)/(double)(Npoints-1);
	    b(i)=b_min + i*(b_max-b_min)/(double)(Npoints-1);
	}
	// Pass parameters to file
        plt_par<<"data_surf=\""<<f_name.str()<<"\"\n";
	plt_par<<"Npoints="<<Npoints<<"\n";
	plt_par<<"a_min="<<a_min<<"\na_max="<<a_max<<"\n";
	plt_par<<"b_min="<<b_min<<"\nb_max="<<b_max<<std::endl;
	// Min value
	double min_a, min_b, min_E=100, min_Eerr;

	// Main loop -------------------------------------------
	for (int k=0; k<Npoints; k+=sqrt(n_conf))
	{
	    for (int l=0; l<Npoints; l+=sqrt(n_conf))
	    {
		acc=0;
		Ecum.zeros(); E2cum.zeros();
		par.init(a(k),a(k+sqrt(n_conf)-1),b(l),b(l+sqrt(n_conf)-1));
		par.init_jastrow(jst,sys);
		for (int i=0; i<Ntot; i++)
		{
		    Ecum_block.zeros(); norm_block.zeros();
		    for (int j=0; j<Nblock; j++)
		    {
			acc+=metropolis(sys,slt,jst,wlk);
			Eloc=calc_local_energy(sys,slt,jst);
			W=reweight(jst);
			Ecum_block+=Eloc*W;
			norm_block+=W;
		    }
		    Eave_block=Ecum_block/norm_block;
		    Ecum+=Eave_block;
		    E2cum+=Eave_block*Eave_block;
		}
		Eave=Ecum/(double)Ntot;
		E2ave=E2cum/(double)Ntot;
		Eerr=sqrt((E2ave-Eave*Eave)/(double)Ntot);
		std::cout<<"\nAcceptance = "<<acc*100.0/(double)(N*Nblock*Ntot)<<"%,\n";
		for (int i=0; i<sqrt(n_conf); i++)
		{
		    for (int j=0; j<sqrt(n_conf); j++)
		    {
			std::cout<<"E(a="<<a(k+i)<<",b="<<b(l+j)<<") = "<<Eave(i*sqrt(n_conf)+j)<<" +- "<<Eerr(i*sqrt(n_conf)+j)<<"\n";
			f<<a(k+i)<<"\t"<<b(l+j)<<"\t"<<Eave(i*sqrt(n_conf)+j)<<"\t"<<Eerr(i*sqrt(n_conf)+j)<<"\n";
			if (Eave(i*sqrt(n_conf)+j)<min_E)
			{
			    min_E=Eave(i*sqrt(n_conf)+j);
			    min_Eerr=Eerr(i*sqrt(n_conf)+j);
			    min_a=a(k+i);
			    min_b=b(l+j);
			}
		    }
		}
	    }
	}

	std::cout<<"\nThe minimum calculated Energy is "<<min_E<<" +- "<<min_Eerr<<" for parameters (a="<<min_a<<",b="<<min_b<<")\n";
    }

    // close files
    f.close();
    plt_par.close();

    // print the execution time
    auto end=std::chrono::high_resolution_clock::now();
    int ex_time=( std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() )*1e-9;
    std::cout<<"\nTotal simulation time = "<<ex_time/3600<<"h "<<(ex_time%3600)/60<<"m "<<(ex_time%3600)%60<<"s\n"<<std::endl;
	
    // free heap
    delete[] jst; 

    return 0;
}

int metropolis (System& sys, Slater& slt, Jastrow* jst, Walker& wlk)
{
    int acc=0;
    double P_acc;
    
    for (int i=0; i<slt.getn_up(); i++)
    {
	wlk.attempt(i,sys);
	P_acc = slt.getratioup(wlk) * jst[Ref].getratio(wlk,sys);
	if (wlk.univar()<P_acc*P_acc)
	{
	    acc++;
	    sys.update(wlk);
	    slt.updateup(wlk,sys);
	    for (int l=0; l<n_conf; l++) jst[l].update(wlk);
	}
    }

    for (int i=slt.getn_up(); i<N; i++)
    {
	wlk.attempt(i,sys);
	P_acc = slt.getratiodown(wlk) * jst[Ref].getratio(wlk,sys);
	if (wlk.univar()<P_acc*P_acc)
	{
	    acc++;
	    sys.update(wlk);
	    slt.updatedown(wlk,sys);
	    for (int l=0; l<n_conf; l++) jst[l].update(wlk);
	}
    }

    return acc;
}

// Calculate the total energy
vector calc_local_energy (System& sys, Slater& slt, Jastrow* jst)
{
    vector E(n_conf);
    E.zeros();
    double Ecomm=0;
    vector grad_slt(2);

    for (int i=0; i<slt.getn_up(); i++)
    {
	for (int j=0; j<i; j++) Ecomm+=1.0/sys.getrel(i,j);
	Ecomm += 0.5*sys.getr2(i) - 0.5*slt.getlapup(i);
	grad_slt=slt.getgradup(i);
	for (int l=0; l<n_conf; l++) E(l)-=0.5*jst[l].getlap(i,sys) + grad_slt.inner(jst[l].getgrad(i,sys));
    }

    for (int i=slt.getn_up(); i<N; i++)
    {
	for (int j=0; j<i; j++) Ecomm+=1.0/sys.getrel(i,j);
	Ecomm += 0.5*sys.getr2(i) - 0.5*slt.getlapdown(i);
	grad_slt=slt.getgraddown(i);
	for (int l=0; l<n_conf; l++) E(l)-=0.5*jst[l].getlap(i,sys) + grad_slt.inner(jst[l].getgrad(i,sys));
    }
    return E+Ecomm;
}

vector reweight (Jastrow* jst)
{
    vector W(n_conf);
    W.zeros();
    for (int i=0; i<N; i++)
    {
	for (int j=0; j<i; j++)
	{
	    for (int l=0; l<Ref; l++) W(l)+=jst[l].getG(i,j)-jst[Ref].getG(i,j);
	    for (int l=Ref+1; l<n_conf; l++) W(l)+=jst[l].getG(i,j)-jst[Ref].getG(i,j);
	}
    }
    W=exp(2*W);
    return W;
}
