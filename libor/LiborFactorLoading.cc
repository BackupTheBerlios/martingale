/* WARANTY NOTICE AND COPYRIGHT
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Copyright (C) Michael J. Meyer

matmjm@mindspring.com
spyqqqdia@yahoo.com

*/



#include "Matrices.h"
#include "Array.h"
#include "LiborFactorLoading.h"
#include "Random.h"
#include "Utils.h"
#include <cmath>
#include <iostream>
using namespace Martingale;





/*******************************************************************************
 *
 *            Volatility Surface
 * 
 ******************************************************************************/
 
 
string 
VolSurface::
volSurfaceType(int volType)
{
     switch(volType){
			 
		 case VolSurface::CONST : return "CONST";
         case VolSurface::JR : return "JR";
		 case VolSurface::M : return "M";
	}
	return "  ";
}


VolSurface* 
VolSurface::
sample(int type)
{
   switch(type){
		   
	   case M     : return M_VolSurface::sample();
	   case JR    : return JR_VolSurface::sample(); 
	   case CONST : return CONST_VolSurface::sample(); 
   }
   std::cout << "\n\nVolSurface::sample(): unknown VolSurface type = " << type
        << ". Exiting.";
   exit(1);
}


// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, const VolSurface& vols){ return vols.printSelf(os); }



// M_VolSurface

Real 
M_VolSurface::
F(Real D, Real s) const { return D*std::exp(s/D); }  
    

Real 
M_VolSurface::
G(Real D, Real s) const { return D*std::exp(s/D)*(s-D); }   


Real 
M_VolSurface::	
H(Real D, Real s) const { return D*std::exp(s/D)*((s-D)*(s-D)+D*D); }  
         

Real 
M_VolSurface::	
g(Real x) const { return 1+a*x*std::exp(-x/d); }

//<---------- Verify integral against quasi Monte Carlo ---------->
Real 
M_VolSurface::
integral_sgsg(Real t, Real T1, Real T2) const     
{
   	Real R,f,f1,f2,D1,D2,D12;
    
    f=a*std::exp(-1/d); f1=f/T1; f2=f/T2;       
	D1=d*T1; D2=d*T2; D12=D1*D2/(D1+D2);
  
    R=t;
    R+=f*F(D1,t);
    R+=f*F(D2,t);
    R+=f*f*F(D12,t);
  
    R-=f1*G(D1,t);
    R-=f2*G(D2,t);
    R-=f*(f1+f2)*G(D12,t);
  
    R+=f1*f2*H(D12,t);      
    //now R=\int g(1-u/a)g(1-u/b)du, see book, 6.11.1
  
    return R;
    
} //end integral_sgsg 


std::ostream& 
M_VolSurface::
printSelf(std::ostream& os) const
{
  	return os << "\n\nVolSurface, type book, 6.11.1, " 
	          << "a=" << a << ", d=" << d << endl;
}
	

VolSurface* 
M_VolSurface::
sample()
{ 
	Real _a=1.5, _d=2.0; 
	return new M_VolSurface(_a,0.0,0.0,_d); 
}

   
     
// JR_VolSurface

Real 
JR_VolSurface::
sigma(Real t, Real T) const 
{ 
   Real s=(T-t);
   return d+(a+b*s)*std::exp(-c*s);
}


//<---------- Verify integral against quasi Monte Carlo ---------->
Real 
JR_VolSurface::
integral_sgsg(Real t, Real T1, Real T2) const
{
   Real f,A,B,C,
        ctmT1=c*(t-T1),
	    ctmT2=c*(t-T2),
        q=ctmT1+ctmT2,                    // c(2t-ti-tj)
        ac=a*c,
        cd=c*d;
              
   f=1.0/(c*c*c);
   A=ac*cd*(std::exp(ctmT2)+std::exp(ctmT1))+c*cd*cd*t;
   B=b*cd*(std::exp(ctmT1)*(ctmT1-1)+std::exp(ctmT2)*(ctmT2-1));
   C=std::exp(q)*(ac*(ac+b*(1-q))+b*b*(0.5*(1-q)+ctmT1*ctmT2))/2;
       
   return f*(A-B+C);
} // end integral_sgsg
   

VolSurface* 
JR_VolSurface::
sample()
{ 
	Real _a=-0.05, _b=0.5, _c=1.5, _d=0.15; 
	return new JR_VolSurface(_a,_b,_c,_d); 
}


std::ostream& 
JR_VolSurface::
printSelf(std::ostream& os) const 
{
   	return
	os << "\n\nVolSurface, type Jaeckel-Rebonato, " 
	   << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d 
	   << endl;
}


// M_VolSurface



// CONST_VolSurface

std::ostream& 
CONST_VolSurface::
printSelf(std::ostream& os) const
{
   	return os << "\n\nVolSurface, type: constant." << endl;
}



/*******************************************************************************
 *
 *                CORRELATIONS
 *
 ******************************************************************************/


Correlations::
Correlations(int _n, Real _alpha, Real _beta, Real _r_oo, int correlationType) : 
n(_n), corrType(correlationType),
alpha(_alpha), beta(_beta), r_oo(_r_oo), 
correlationMatrix(n-1,1) 
{   }

   
Correlations* 
Correlations::
sample(int n, int type)
{
    switch(type){
		   
	    case JR : return JR_Correlations::sample(n); 
	    default : return CS_Correlations::sample(n); 
    }
}


Real&
Correlations::
operator()(int i, int j) { return correlationMatrix(i,j); }


   
string 
Correlations::
correlationType(int corrType)
{
     switch(corrType){
			 
		 case Correlations::JR : return "JR";
         case Correlations::CS : return "CS";
		 default :               return "Unknown-Correlations";
	}
	return "   ";   // makes compiler happy
}

// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, const Correlations& vols){ return vols.printSelf(os); }



// JR_Correlations

JR_Correlations::
JR_Correlations(const RealArray1D& _T, Real beta) : 
Correlations(_T.getDimension()-1,0.0,beta,0.0,JR), T(_T) 	
{ 
	setCorrelations(); 
}


void 
JR_Correlations::
setCorrelations()
{
    for(int i=1;i<n;i++)
    for(int j=i;j<n;j++) 
	   correlationMatrix(i,j)=std::exp(beta*(T[i]-T[j])); 
}


Correlations* 
JR_Correlations::
sample(int n, Real delta=0.25)
{ 
	RealArray1D T(n+1); T[0]=0.0;
	for(int i=0;i<n;i++) T[i+1]=T[i]+delta;
		
	Real beta=0.1;
	return new JR_Correlations(T,beta); 
}
	

std::ostream& 
JR_Correlations::
printSelf(std::ostream& os) const
{
	return
	os << "\n\nCorrelations: Jaeckel-Rebonato, " 
	   << "beta=" << beta << endl;
}



// CS_Correlations

void 
CS_Correlations::
setCorrelations()
{
     RealArray1D b(n-1,1); 
     for(int i=1;i<n;i++)
     { Real x=((Real)i)/(n-1); b[i]=std::exp(-f(x)); }
           
     // initialize the correlation matrix rho
     for(int i=1;i<n;i++)
     for(int j=i;j<n;j++) correlationMatrix(i,j)=b[i]/b[j];  
}
	
		
Real 
CS_Correlations::
f(Real x) const 
{
    Real a=alpha/2, b=(beta-alpha)/6, c=log(r_oo)-(a+b);
    return x*(c+x*(a+b*x));  //cx+ax^2+bx^3
}


std::ostream& 
CS_Correlations::
printSelf(std::ostream& os) const
{
	return
	os << "\n\nCorrelations: Coffee-Shoenmakers, " 
	   << "alpha=" << alpha << "beta=" << beta << "r_oo=" << r_oo 
	   << endl;
}


Correlations* 
CS_Correlations::
sample(int _n)
{ 
	Real _alpha=1.8, _beta=0.1, _r_oo=0.4;
	return new CS_Correlations(_n,_alpha,_beta,_r_oo);
}





/*******************************************************************************
 *
 *                     LiborFactorLoading
 * 
 ******************************************************************************/


 
LiborFactorLoading::
LiborFactorLoading
(const RealArray1D& L0, const RealArray1D& deltas, const RealArray1D& _k,
VolSurface* vols, Correlations* corrs) :
n(L0.getDimension()), 
delta(deltas), T(n+1), l(L0), x(n), k(_k),
vol(vols), corr(corrs), rho(*corrs)
{   
    // tenor structure, initial XLibors
	T[0]=0;
	for(int i=0;i<n;i++){ 
			
		T[i+1]=T[i]+delta[i];
		x[i]=delta[i]*l[i];
	}
		
	if(n!=corrs->getDimension()){
	    std::cout << "\n\nLiborfactorLoading(): Libors-Correlations dimensions don't match."
		          << "Libors n = " << n << ", correlations n = " << corrs->getDimension()
	              << "\nTerminating.";
	    exit(1);
	}

} // end constructor
	
	
	
LiborFactorLoading* 
LiborFactorLoading::
sample(int n, int volType=VolSurface::JR, int corrType=Correlations::CS)
{
	RealArray1D deltas(n);
	RealArray1D c(n);        // k_j
	RealArray1D L0(n);
    for(int i=0;i<n;i++)
	{ deltas[i]=0.25; c[i]=0.2+0.1*Random::U01(); L0[i]=0.04; }
				
	VolSurface* vols = VolSurface::sample(volType);
	Correlations* corrs = Correlations::sample(n,corrType);
			
	return new LiborFactorLoading(L0,deltas,c,vols,corrs);
}
	
	
	
std::ostream& 
LiborFactorLoading::
printSelf(std::ostream& os) const
{
	RealVector L0(n,l.getData());      // initial Libors		
	RealVector vols(n-1,1);            // annualized volatilities
	for(int i=1;i<n;i++) vols[i]=annualVol(i);
		
	return 
	os << "\n\nLiborFactorLoading: dimension = " << n 
	   << "Initial Libors: " << L0
	   << "\nAnnualized volatilities: " << vols
	   << "\nVolatility surface: " << vol
	   << "\nCorrelations: " << corr;
}
	
	
Real 
LiborFactorLoading::	
sigma(int i, Real t) const 
{ 
	return k[i]*(vol->sigma(t,T[i])); 
}
	
	
Real 
LiborFactorLoading::
annualVol(int i) const
{
	if(i==0) return 0.0;
	Real   T_i=T[i],
	       volSqr=integral_sgi_sgj_rhoij(i,i,0.0,T_i);
	return sqrt(volSqr/T_i);
}
					

// LOG-COVARIATION MATRICES 
	
	
Real 
LiborFactorLoading::
integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const
{  
   return k[i]*k[j]*rho(i,j)*(vol->integral_sgsg(s,t,T[i],T[j])); 
}

   	
	
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrix(int p,int q, Real s, Real t) const
{
    int size=q-p;       // matrix size
	UTRRealMatrix& logCVM=*(new UTRRealMatrix(size,p));

    for(int i=p;i<q;i++)
    for(int j=i;j<q;j++)
       logCVM(i,j)=integral_sgi_sgj_rhoij(i,j,s,t);

    return logCVM;
}// end logCovariationMatrix
   
   
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrix(int t) const
{
   return logLiborCovariationMatrix(t+1,n,T[t],T[t+1]);
}
   
   
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrixRoot(int t) const
{
   return logLiborCovariationMatrix(t).utrRoot();
}


const RealMatrix& 
LiborFactorLoading::
reducedRankLogLiborCovariationMatrixRoot(int t, int r) const
{
   return logLiborCovariationMatrix(t).rankReducedRoot(r);
}

   
// TEST PROGRAMS 
		 
void 
LiborFactorLoading::
selfTest() const
{
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
    std::cout << "\nLIBOR FACTORLOADING SELFTEST:" << endl << this;
	
	cout << "\nTesting the root L of the matrix C=logLiborCovariationMatrix(t):" << endl;
	for(int t=0;t<n-1;t++){
		
		const UTRRealMatrix& C=logLiborCovariationMatrix(t);
		const UTRRealMatrix& L=logLiborCovariationMatrixRoot(t);
		cout << "\nt = " << t << ": ";
		C.testEquals(L.aat(),precision,epsilon,"C=L*L'");
	}

} // end test


void 
LiborFactorLoading::
factorizationTest(int r) const
{
	cerr << "\n\nApproximate factorization of all single time step log-Libor"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank " << r
	     << "\n\nRelative errors (trace norm): " << endl;
	
	for(int t=0;t<n-2;t++){
		
		const UTRRealMatrix& Ct=logLiborCovariationMatrix(t);
		Ct.testFactorization(r);
	}
}


// CALIBRATION

// use the first coordinates of x to set the Libor volatility scaling factors 
// k[j] (including k[0] although its useless) then the next coordinates to set 
// VolSurface::a,b,c,d and the follwing coordinates to set 
// Correlations::alpha,beta,r_oo
void 
LiborFactorLoading::
setParameters(const RealArray1D& X)
{
	Real* u=X.getData();
	for(int j=0;j<n;j++) k[j]=u[j];
	vol->setParameters(u[n],u[n+1],u[n+2],u[n+3]);
	corr->setParameters(u[n+4],u[n+5],u[n+6]);
}
	

// Global Insertion

std::ostream& operator << 
(std::ostream& os, const LiborFactorLoading& vols){ return vols.printSelf(os); }









