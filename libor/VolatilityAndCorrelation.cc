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
#include "VolatilityAndCorrelation.h"
#include "Random.h"
#include "Utils.h"
#include <string>
#include <cmath>
#include <iostream>

MTGL_BEGIN_NAMESPACE(Martingale)





/*******************************************************************************
 *
 *            Volatility Surface
 * 
 ******************************************************************************/
 
 
std::string 
VolSurface::
volSurfaceType()
{
     switch(volType){
			 
		 case CONST : return "CONST";
         case JR    : return "JR";
		 case M     : return "M";
	}
	return "Unknown volatility surface";
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


// TEST OF VOLATILITY INTEGRALS

void 
VolSurface::
testVolSurfaceIntegrals(int N, Real precision)
{
	std::cout << "\n\n\n\nTesting volatility integrals: " << endl << *this;
		
	Real analytic, montecarlo, sum, ft, u, t,
	     sg_tT1, sg_tT2;  // sigma(t,T)
	
	// loop over some integration intervals [a,b] (b=a+i*0.3)
	for(Real a=1.0;a<7.0;a+=0.5)
	for(int i=1;i<5;i++){
		
		Real b=a+i*0.3,
		     T1=b+1.23, T2=b+2.13;
	    // Monte  carlo expectation
	    sum=0;
	    for(int i=0;i<N;i++){
					   
		    u=Random::U01();
			t=a+u*(b-a);
			sg_tT1=sigma(t,T1);
			sg_tT2=sigma(t,T2);
			ft=sg_tT1*sg_tT2;
		        
			sum+=ft;
	    }
	    montecarlo=(b-a)*sum/N;
		analytic=integral_sgsg(a,b,T1,T2);
		
		Real error=100*abs(montecarlo-analytic)/analytic;
		if(error>precision){
			
		     std::cout << "\n\nTest failed."
			           << "\nAnalytic integral: " << analytic
			           << "\nMonte Carlo integral: " << montecarlo
			           << "\nRelative error (%): " << error
			           << "Aborting.";
			exit(1);
		}
	} // end main loop
	std::cout << "\n\n Test passed.";
} // end testVolSurfaceIntegrals
			    
			         
			






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
  	return os << "\nVolSurface, type book, 6.11.1, " 
	          << "a=" << a << ", d=" << d;
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
	os << "\nVolSurface, type Jaeckel-Rebonato, " 
	   << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d;
}


// M_VolSurface



// CONST_VolSurface

std::ostream& 
CONST_VolSurface::
printSelf(std::ostream& os) const
{
   	return os << "\nVolSurface, type: constant.";
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
correlationType()
{
     switch(corrType){
			 
		 case Correlations::JR : return "JR"; 
         case Correlations::CS : return "CS";   
	}
	return "Unknown-Correlations";
}


void 
Correlations::
setParameters(Real _alpha, Real _beta, Real _r_oo)
{ 
	alpha=_alpha; beta=_beta; r_oo=_r_oo; 
    setCorrelations();
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
	os << "\nCorrelations: Jaeckel-Rebonato, " 
	   << "beta=" << beta;
}



// CS_Correlations
CS_Correlations::
CS_Correlations(int _n, Real _alpha, Real _beta, Real _r_oo) : 
Correlations(_n,_alpha,_beta,_r_oo,CS) 	
{ 
	setCorrelations(); 
}


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
	os << "\nCorrelations: Coffee-Shoenmakers, " 
	   << "alpha=" << alpha << ", beta=" << beta << ", r_oo=" << r_oo;
}


Correlations* 
CS_Correlations::
sample(int _n)
{ 
	Real _alpha=1.8, _beta=0.1, _r_oo=0.4;
	return new CS_Correlations(_n,_alpha,_beta,_r_oo);
}


MTGL_END_NAMESPACE(Martingale)
