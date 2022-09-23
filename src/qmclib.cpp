#include <cstdlib>
#include <complex>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <fstream>
#include"vector.h"
#include"matrix.h"
#include"qmclib.h"

// Parameters class functions

Parameters::Parameters (int n) {n_conf=n; allocate();}

void Parameters::allocate () 
{
    p = new vector[n_conf];
    for (int i=0; i<n_conf; i++) p[i].redim(2);
    p = new vector[n_conf];
    for (int i=0; i<n_conf; i++) p[i].redim(2);
    grad_old.redim(2);
    grad_new.redim(2);
    delta_grad.redim(2);
    p_old.redim(2);
    s.redim(2);
    delta_p.redim(2);
    H.redim(2,2);
}

void Parameters::init (double a_min, double a_max, double b_min, double b_max)
{
    if (n_conf>1)
    {
	for (int i=0; i<sqrt(n_conf); i++)
	{
	    for (int j=0; j<sqrt(n_conf); j++)
	    {
		p[i*(int)sqrt(n_conf)+j](0)=a_min + i*(a_max-a_min)/(double)(sqrt(n_conf)-1);
		p[i*(int)sqrt(n_conf)+j](1)=b_min + j*(b_max-b_min)/(double)(sqrt(n_conf)-1);
	    }
	}
    }
    else
    {
	p[0](0)=a_min;
	p[0](1)=b_min;
    }
}

void Parameters::init_dfp (double guess_a, double guess_b, double g_step, double threshold_)
{
    if (n_conf==5)
    {
	count=0;
	next=true;
	grad_step=g_step;
	threshold=threshold_;
	p[0](0)=guess_a;	p[0](1)=guess_b;
	for (int i=1; i<n_conf; i++)
	{
	    p[i]=p[0];
	    p[i]((i-1)/2)+=pow(-1.0,i)*grad_step;
	}
	p_old=p[0];
	H.zeros();
	for (int i=0; i<2; i++) H(i,i)=1.0;
    }
    else std::cout<<"Error! For gradient computation n_conf=3 is required\n";
}

void Parameters::init_jastrow (Jastrow* jst, System& sys)
{for (int i=0; i<n_conf; i++) jst[i].init(sys,p[i](0),p[i](1));}

bool Parameters::move (vector E)
{
    bool k=true;
    if (next)
    {
	E_0=E(0);
	std::cout<<"count="<<count<<"\n";
	for (int i=1; i<n_conf; i+=2)
	    grad_new(i/2)=(E(i+1)-E(i))/(2*grad_step);
	std::cout<<"grad="<<grad_new<<"\n";
	if (grad_new.norm2()<=threshold*threshold) k=false;
	if (count>0) updateH();
	alpha=grad_step;
	s=-1.0*H.inner(grad_new);
	std::cout<<"s="<<s<<"\n";
	der_0=grad_new.inner(s);
	std::cout<<"der_0="<<der_0<<"\n";
	next=false;
	count1D=0;
	p[0]=p_old+alpha*s;
	p[1]=p[0]-2*grad_step*s;
	p[2]=p[0]-grad_step*s;
	p[3]=p[0]+grad_step*s;
	p[4]=p[0]+2*grad_step*s;
    }
    else
    {
	der=calcder(E);
	std::cout<<"der="<<der<<"\n";
	der2=calcder2(E);
	std::cout<<"der2="<<der2<<"\n";
	if (der2<0)
	{
	    if (der<=0) alpha+=1.0;
	    if (der>=0) alpha-=1.0;
	    std::cout<<"alpha="<<alpha<<"\n";
	    p[0]=p_old+alpha*s;
	    p[1]=p[0]-2*grad_step*s;
	    p[2]=p[0]-grad_step*s;
	    p[3]=p[0]+grad_step*s;
	    p[4]=p[0]+2*grad_step*s;
	}
	else
	{
	    if (der*der<=threshold*threshold)
	    {
		if (E(0)<E_0)
		{
		    next=true;
		    iterate();
		}
		else
		{
		    k=false;
		    std::cout<<"Error: Accuracy too low!\nExiting\n";
		}
	    }
	    else
	    {
		alpha-=der/der2;
		std::cout<<"alpha="<<alpha<<"\n";
		p[0]=p_old+alpha*s;
		p[1]=p[0]-2*grad_step*s;
		p[2]=p[0]-grad_step*s;
		p[3]=p[0]+grad_step*s;
		p[4]=p[0]+2*grad_step*s;
	    }
	}
	count1D++;
    }
    count++;
    return k;
}

void Parameters::iterate()
{
    p_old=p[0];
    grad_old=grad_new;
    for (int i=1; i<n_conf; i++)
    {
	p[i]=p[0];
	p[i]((i-1)/2)+=pow(-1.0,i)*grad_step;
    }
}

void Parameters::updateH ()
{
	// delta gradient
	delta_grad=grad_new-grad_old;
	std::cout<<"delta_grad="<<delta_grad<<"\n";
	// delta point
	delta_p=alpha*s;
	std::cout<<"delta_p="<<delta_p<<"\n";
	// update hessian inverse
	H+=( (outer(delta_p,delta_p))/(delta_p.inner(delta_grad)) ) - ( (outer(H.inner(delta_grad),H.inner(delta_grad)))/(delta_grad.inner(H.inner(delta_grad))) );
}

double Parameters::calcder (vector E) {return (E(1)-8*E(2)+8*E(3)-E(4))/(12*grad_step);}

double Parameters::calcder2 (vector E) {return (-E(1)+16*E(2)-30*E(0)+16*E(3)-E(4))/(12*grad_step*grad_step);}

std::ostream &operator<< (std::ostream& out, const Parameters& par) {out<<"("<<par.geta()<<", "<<par.getb()<<")"; return out;}

double Parameters::geta () const {return p[0](0);}

double Parameters::getb () const {return p[0](1);}

Parameters::~Parameters () {deallocate();}

void Parameters::deallocate () {delete[] p;}

// System class functions

System::System (int n) {N=n; allocate();}

void System::allocate ()
{
    p=new vector[N];
    for (int i=0; i<N; i++) p[i].redim(2);
    rel.redim(N,N);
    r2.redim(N);
}

void System::init(double a)
{
    // set initial positions
    switch (N)
    {
	case 2 :
	    p[0](0)=-a/2.;  p[0](1)=-a/2.;
	    p[1](0)=a/2.;   p[1](1)=a/2.;
	    break;
	case 3 :
	    p[0](0)=-a/2.;  p[0](1)=-a/2.;
	    p[1](0)=-a/2.;  p[1](1)=a/2.;
	    p[2](0)=a;	    p[2](1)=0;
	    break;
	case 4 :
	    p[0](0)=-a;	    p[0](1)=-a;
	    p[1](0)=a;	    p[1](1)=a;
	    p[2](0)=-a;	    p[2](1)=a;
	    p[3](0)=a;	    p[3](1)=-a;
	    break;
	case 5 :
	    p[0](0)=-a;	    p[0](1)=-a;
	    p[1](0)=a;	    p[1](1)=a;
	    p[2](0)=-a;	    p[2](1)=a;
	    p[3](0)=a;	    p[3](1)=-a;
	    p[4](0)=0;	    p[4](1)=0;
	    break;
	case 6 :
	    p[0](0)=-a;	    p[0](1)=-a;
	    p[1](0)=a;	    p[1](1)=a;
	    p[2](0)=-a;	    p[2](1)=a;
	    p[3](0)=a;	    p[3](1)=-a;
	    p[4](0)=0;	    p[4](1)=0;
	    p[5](0)=2*a;    p[5](1)=0;
	    break;
    }
    srand(997);
    for (int i=0; i<N; i++)
    {
	vector tmp(2);
	tmp(0)=0.05*a*rand()/(double)RAND_MAX;
	tmp(1)=0.05*a*rand()/(double)RAND_MAX;
	p[i]+=tmp;
    }
    // fill matrices
    rel.zeros();
    for (int i=0; i<N; i++)
    {
	for (int j=0; j<i; j++)
	{
	    rel(i,j)=p[i].rel(p[j]);
	}
	r2(i)=p[i].norm2();
    }
}

double System::psi(int* qn, int cn)
{
    double psi_=0;
    psi_=H(qn[0],p[cn](0))*H(qn[1],p[cn](1))*exp(-0.5*r2(cn));
    return psi_;
}

vector System::gradpsi(int* qn, int cn)
{
    vector gradpsi_(2);
    gradpsi_(0)=(2*qn[0]*H(qn[0]-1,p[cn](0))/H(qn[0],p[cn](0)) - p[cn](0));
    gradpsi_(1)=(2*qn[1]*H(qn[1]-1,p[cn](1))/H(qn[1],p[cn](1)) - p[cn](1));
    return gradpsi_;
}

double System::lappsi(int* qn, int cn)
{
    double lappsi_=0;
    double Hx=H(qn[0],p[cn](0));
    double Hy=H(qn[1],p[cn](1));
    lappsi_=(4*qn[0]*(qn[0]-1)*H(qn[0]-2,p[cn](0))/Hx + 4*qn[1]*(qn[1]-1)*H(qn[1]-2,p[cn](1))/Hy + r2(cn) - 4*qn[0]*p[cn](0)*H(qn[0]-1,p[cn](0))/Hx - 4*qn[1]*p[cn](1)*H(qn[1]-1,p[cn](1))/Hy - 2);
    return lappsi_;
}

double System::getrel(int i, int j) {return rel(i,j);}

double System::getr2(int i) {return r2(i);}

void System::update (Walker& wlk) 
{
    int curr=wlk.current();
    p[curr]=wlk.p;
    r2(curr)=wlk.getr2();
    for (int i=0; i<curr; i++) rel(curr,i)=wlk.rel(i);
    for (int i=curr+1; i<N; i++) rel(i,curr)=wlk.rel(i);
}
System::~System () {deallocate();}

void System::deallocate () {delete[] p;}

// Slater class functions

Slater::Slater (int n) 
{
    N=n;
    qn=new int*[N];
    for (int i=0; i<N; i++) qn[i]=new int[2];

    switch (N)
    {
	case 2:
	    n_up=1;
	    qn[0][0]=0;	    qn[0][1]=0;
	    qn[1][0]=0;	    qn[1][1]=0;
	    break;
	case 3:
	    n_up=2;
	    qn[0][0]=0;	    qn[0][1]=0;
	    qn[1][0]=1;	    qn[1][1]=0;
	    qn[2][0]=0;	    qn[2][1]=0;
	    break;
	case 4:				    //TODO
	{
	    int L=0, S=0, count=0;
	    std::cout<<"For N=4 there are 3 possibilities for L (total angular momentum) and S (total spin).\n";
		std::cout<<"Select of which one you want to calculate the energy:\n";
	    while (((L!=0&&S!=0)&&(L!=0&&S!=1)&&(L!=2&&S!=0))||count<1)
	    {
		if (count>0) std::cout<<"invalid configuration..\n choose from:\n";
		std::cout<<"L=0, S=1\tL=0, S=0\tL=2, S=0\n";
		std::cout<<"L=";
		std::cin>>L;
		std::cout<<"S=";
		std::cin>>S;
		std::cout<<"\r";
		count++;
	    }
	    qn[0][0]=0;	    qn[0][1]=0;
	    qn[1][0]=1;	    qn[1][1]=0;
	    if (S==1)
	    {
		n_up=3;
		qn[2][0]=0;	    qn[2][1]=1;
		qn[3][0]=0;	    qn[3][1]=0;
	    }
	    else
	    {
		n_up=2;
		qn[2][0]=0;	    qn[2][1]=0;
		if (L==0)
		{
		    qn[3][0]=0;	    qn[3][1]=1;
		}
		else
		{
		    qn[3][0]=1;	    qn[3][1]=0;
		}
	    }
	    break;
	}
	case 5:
	    n_up=3;
	    qn[0][0]=0;	    qn[0][1]=0;
	    qn[1][0]=1;	    qn[1][1]=0;
	    qn[2][0]=0;	    qn[2][1]=1;
	    qn[3][0]=0;	    qn[3][1]=0;
	    qn[4][0]=1;	    qn[4][1]=0;
	    break;
	case 6:
	    n_up=3;
	    qn[0][0]=0;	    qn[0][1]=0;
	    qn[1][0]=1;	    qn[1][1]=0;
	    qn[2][0]=0;	    qn[2][1]=1;
	    qn[3][0]=0;	    qn[3][1]=0;
	    qn[4][0]=1;	    qn[4][1]=0;
	    qn[5][0]=0;	    qn[5][1]=1;
	    break;
    }

    allocate();
}

void Slater::allocate ()
{
    Dup.redim(n_up,n_up);
    Ddown.redim(N-n_up,N-n_up);
    gradDup=new vector*[n_up];
    for (int i=0; i<n_up; i++)
    {
	gradDup[i]=new vector[n_up];
	for (int j=0; j<n_up; j++) gradDup[i][j].redim(2);
    }
    gradDdown=new vector*[N-n_up];
    for (int i=0; i<N-n_up; i++)
    {
	gradDdown[i]=new vector[N-n_up];
	for (int j=0; j<N-n_up; j++) gradDdown[i][j].redim(2);
    }
    lapDup.redim(n_up,n_up);
    lapDdown.redim(N-n_up,N-n_up);
    Dup_inv.redim(n_up,n_up);
    Ddown_inv.redim(N-n_up,N-n_up);
}

void Slater::init(System& sys)
{
    // set quantum numbers
    // fill matrices
    for (int i=0; i<n_up; i++)
    {
	for (int j=0; j<n_up; j++)
	{
	    Dup(i,j)=sys.psi(qn[j],i);
	    gradDup[i][j]=sys.gradpsi(qn[j],i)*Dup(i,j);
	    lapDup(i,j)=sys.lappsi(qn[j],i)*Dup(i,j);
	}
    }
    for (int i=0; i<N-n_up; i++)
    {
	for (int j=0; j<N-n_up; j++)
	{
	    Ddown(i,j)=sys.psi(qn[j+n_up],i+n_up);
	    gradDdown[i][j]=sys.gradpsi(qn[j+n_up],i+n_up)*Ddown(i,j);
	    lapDdown(i,j)=sys.lappsi(qn[j+n_up],i+n_up)*Ddown(i,j);
	}
    }
    Dup_inv=Dup.inverse();
    Ddown_inv=Ddown.inverse();
}
int Slater::getn_up () {return n_up;}

double Slater::getratioup (Walker& wlk)
{
    ratio=0;
    for (int i=0; i<n_up; i++)
    {
	wlk.calcpsi(i,qn[i]);
	ratio+=wlk.psi(i)*Dup_inv(i,wlk.current());
    }
    return ratio;
}

double Slater::getratiodown (Walker& wlk)
{
    ratio=0;
    for (int i=0; i<N-n_up; i++)
    {
	wlk.calcpsi(i+n_up,qn[i+n_up]);
	ratio+=wlk.psi(i+n_up)*Ddown_inv(i,wlk.current()-n_up);
    }
    return ratio;
}

vector Slater::getgradup (int i)
{
    vector grad(2);
    grad.zeros();
    for (int j=0; j<n_up; j++)
    {
	grad+=gradDup[i][j]*Dup_inv(j,i);
    }
    return grad;
}

vector Slater::getgraddown (int i)
{
    vector grad(2);
    grad.zeros();
    for (int j=0; j<N-n_up; j++)
    {
	grad+=gradDdown[i-n_up][j]*Ddown_inv(j,i-n_up);
    }
    return grad;
}

double Slater::getlapup (int i)
{
    double lap=0;
    for (int j=0; j<n_up; j++)
    {
	lap+=lapDup(i,j)*Dup_inv(j,i);
    }
    return lap;
}

double Slater::getlapdown (int i)
{
    double lap=0;
    for (int j=0; j<N-n_up; j++)
    {
	lap+=lapDdown(i-n_up,j)*Ddown_inv(j,i-n_up);
    }
    return lap;
}

void Slater::updateup (Walker& wlk, System& sys)
{
    int curr=wlk.current();
    double sum;
    for (int i=0; i<curr; i++)
    {
	sum=0;
	for (int j=0; j<n_up; j++) sum+=wlk.psi(j)*Dup_inv(j,i);
	for (int j=0; j<n_up; j++)Dup_inv(j,i)-=Dup_inv(j,curr)*sum/ratio;
    }
    for (int i=curr+1; i<n_up; i++)
    {
	sum=0;
	for (int j=0; j<n_up; j++) sum+=wlk.psi(j)*Dup_inv(j,i);
	for (int j=0; j<n_up; j++)Dup_inv(j,i)-=Dup_inv(j,curr)*sum/ratio;
    }
    for (int i=0; i<n_up; i++)
    {
	Dup(curr,i)=wlk.psi(i);	// Maybe for N=4
	Dup_inv(i,curr)=Dup_inv(i,curr)/ratio;
	gradDup[curr][i]=sys.gradpsi(qn[i],curr)*wlk.psi(i);
	lapDup(curr,i)=sys.lappsi(qn[i],curr)*wlk.psi(i);
    }
}

void Slater::updatedown (Walker& wlk, System& sys)
{
    int curr=wlk.current();
    double sum;
    for (int i=0; i<curr-n_up; i++)
    {
	sum=0;
	for (int j=0; j<N-n_up; j++) sum+=wlk.psi(j+n_up)*Ddown_inv(j,i);
	for (int j=0; j<N-n_up; j++)Ddown_inv(j,i)-=Ddown_inv(j,curr-n_up)*sum/ratio;
    }
    for (int i=curr+1-n_up; i<N-n_up; i++)
    {
	sum=0;
	for (int j=0; j<N-n_up; j++) sum+=wlk.psi(j+n_up)*Ddown_inv(j,i);
	for (int j=0; j<N-n_up; j++)Ddown_inv(j,i)-=Ddown_inv(j,curr-n_up)*sum/ratio;
    }
    for (int i=0; i<N-n_up; i++)
    {
	Ddown(curr-n_up,i)=wlk.psi(i+n_up);	// Maybe for N=4
	Ddown_inv(i,curr-n_up)=Ddown_inv(i,curr-n_up)/ratio;
	gradDdown[curr-n_up][i]=sys.gradpsi(qn[i+n_up],curr)*wlk.psi(i+n_up);
	lapDdown(curr-n_up,i)=sys.lappsi(qn[i+n_up],curr)*wlk.psi(i+n_up);
    }
}

Slater::~Slater () {deallocate();}

void Slater::deallocate () 
{
    for (int i=0; i<N; i++) delete[] qn[i];
    delete[] qn;
    for (int i=0; i<n_up; i++) delete[] gradDup[i];
    delete[] gradDup;
    for (int i=0; i<N-n_up; i++) delete[] gradDdown[i];
    delete[] gradDdown;
}

// Jastrow class functions

Jastrow::Jastrow () {}

void Jastrow::allocate (int n)
{
    N=n;
    G.redim(N,N);
    grad.redim(2);
}

void Jastrow::init(System& sys, double A, double B)
{
    a=A; b=B;
    // fill matrices
    G.zeros();
    for (int i=0; i<N; i++)
    {
	for (int j=0; j<i; j++) G(i,j)=g(sys.getrel(i,j));
    }
}

double Jastrow::g(double r_ij)
{
    double g_;
    g_=a*r_ij/(1+b*r_ij);
    return g_;
}

double Jastrow::getG (int i, int j) const {return G(i,j);}

double Jastrow::getratio (Walker& wlk, System& sys)
{
    double delta_G=0;
    int curr=wlk.current();
    for (int i=0; i<curr; i++)
    {
	wlk.rel(i)=wlk.p.rel(sys.p[i]);
	delta_G+=g(wlk.rel(i))-G(curr,i);
    }
    for (int i=curr+1; i<N; i++)
    {
	wlk.rel(i)=sys.p[i].rel(wlk.p);
	delta_G+=g(wlk.rel(i))-G(i,curr);
    }
    return exp(delta_G);
}

vector Jastrow::getgrad (int i, System& sys)
{
    return grad;
}

double Jastrow::getlap (int i, System& sys)
{
    grad.zeros();
    double lap=0;
    double r_ij, fact;
    for (int j=0; j<i; j++)
    {
	r_ij=sys.getrel(i,j);
	fact=(a/(r_ij*(1+b*r_ij)*(1+b*r_ij)));
	grad+=(sys.p[i]-sys.p[j])*fact;
	lap+=fact*(1-b*r_ij)/(1+b*r_ij);
    }
    for (int j=i+1; j<N; j++)
    {
	r_ij=sys.getrel(j,i);
	fact=(a/(r_ij*(1+b*r_ij)*(1+b*r_ij)));
	grad+=(sys.p[i]-sys.p[j])*fact;
	lap+=fact*(1-b*r_ij)/(1+b*r_ij);
    }
    return lap + grad.inner(grad);
}

void Jastrow::update (Walker& wlk)
{
    int curr=wlk.current();
    for (int i=0; i<curr; i++) G(curr,i)=g(wlk.rel(i));
    for (int i=curr+1; i<N; i++) G(i,curr)=g(wlk.rel(i));
}

Jastrow::~Jastrow () {}

// Walker class functions

Walker::Walker (int n) {N=n; allocate(); srand(997);}

void Walker::allocate ()
{
    p.redim(2);
    rel.redim(N);
    psi.redim(N);
}

void Walker::setdelta (double d) {delta=d;}

double Walker::getr2 () {return r2;}

int Walker::current () {return curr;}

void Walker::calcpsi (int i, int* qn)
{
    psi(i)=H(qn[0],p(0))*H(qn[1],p(1))*exp(-0.5*r2);
}

void Walker::attempt (int i, System& sys)
{
    curr=i;
    p(0) = sys.p[i](0) + delta*(rand()/(double)RAND_MAX - 0.5);
    p(1) = sys.p[i](1) + delta*(rand()/(double)RAND_MAX - 0.5);
    r2=p.norm2();
}

double Walker::univar () {return rand()/(double)RAND_MAX;}

Walker::~Walker () {}

