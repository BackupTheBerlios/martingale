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



//#include "Matrices.h"
#include "Matrix.h"
#include "Array.h"
#include "VolatilityAndCorrelation.h"
#include "Random.h"
#include "Utils.h"
#include <string>
#include <cmath>
#include <iostream>

MTGL_BEGIN_NAMESPACE(Martingale)


using std::exp;
using std::ostream;
using std::string;
using std::cout;
using std::endl;





/*******************************************************************************
 *
 *            Volatility Surface
 * 
 ******************************************************************************/
 
 
string 
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


string 
VolSurface::
volSurfaceType(int type)
{
     switch(type){
			 
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
   cout << "\n\nVolSurface::sample(): unknown VolSurface type = " << type
        << ". Exiting.";
   exit(1);
}



// TEST OF VOLATILITY INTEGRALS

void 
VolSurface::
testVolSurfaceIntegrals(int N, Real precision)
{
	cout << "\n\n\n\nTesting volatility integrals: " << endl << *this;
		
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
			
		     cout << "\n\nTest failed."
	              << "\nAnalytic integral: " << analytic
	              << "\nMonte Carlo integral: " << montecarlo
		          << "\nRelative error (%): " << error
		          << "Aborting.";
			 exit(1);
		}
	} // end main loop
	cout << "\n\n Test passed.";
} // end testVolSurfaceIntegrals
			    
			         
			






// M_VolSurface

Real 
M_VolSurface::
F(Real D, Real s) const { return D*exp(s/D); }  
    

Real 
M_VolSurface::
G(Real D, Real s) const { return D*exp(s/D)*(s-D); }   


Real 
M_VolSurface::	
H(Real D, Real s) const { return D*exp(s/D)*((s-D)*(s-D)+D*D); }  
         

Real 
M_VolSurface::	
g(Real x) const { return 1+a*x*exp(-x/d); }

	
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


ostream& 
M_VolSurface::
printSelf(ostream& os) const
{
  	return os << "\nVolSurface, type book, 6.11.1, " 
	          << "a=" << a << ", d=" << d;
}
	

VolSurface* 
M_VolSurface::
sample()
{ 
	Real a_=1.5, d_=2.0; 
	return new M_VolSurface(a_,0.0,0.0,d_); 
}

   
     
// JR_VolSurface

Real 
JR_VolSurface::
sigma(Real t, Real T) const 
{ 
   Real s=(T-t);
   return d+(a+b*s)*exp(-c*s);
}


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
   A=ac*cd*(exp(ctmT2)+exp(ctmT1))+c*cd*cd*t;
   B=b*cd*(exp(ctmT1)*(ctmT1-1)+exp(ctmT2)*(ctmT2-1));
   C=exp(q)*(ac*(ac+b*(1-q))+b*b*(0.5*(1-q)+ctmT1*ctmT2))/2;
       
   return f*(A-B+C);
} // end integral_sgsg
   

VolSurface* 
JR_VolSurface::
sample()
{ 
	Real a_=-0.05, b_=0.5, c_=1.5, d_=0.15; 
	return new JR_VolSurface(a_,b_,c_,d_); 
}


ostream& 
JR_VolSurface::
printSelf(ostream& os) const 
{
   	return
	os << "\nVolSurface, type Jaeckel-Rebonato, " 
	   << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d;
}


// M_VolSurface



// CONST_VolSurface

ostream& 
CONST_VolSurface::
printSelf(ostream& os) const
{
   	return os << "\nVolSurface, type: constant.";
}



/*******************************************************************************
 *
 *                CORRELATIONS
 *
 ******************************************************************************/


Correlations::
Correlations(int n_, Real alpha_, Real beta_, Real r_oo_, int correlationType) : 
n(n_), corrType(correlationType),
alpha(alpha_), beta(beta_), r_oo(r_oo_), 
correlationMatrix(n-1,1) 
{   }

   
Correlations* 
Correlations::
sample(int m, int type)
{
    switch(type){
		   
	    case JR : return JR_Correlations::sample(m); 
	    default : return CS_Correlations::sample(m); 
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


string 
Correlations::
correlationType(int type)
{
     switch(type){
			 
		 case Correlations::JR : return "JR"; 
         case Correlations::CS : return "CS";   
	}
	return "Unknown-Correlations";
}



void 
Correlations::
setParameters(Real alpha_, Real beta_, Real r_oo_)
{ 
	alpha=alpha_; beta=beta_; r_oo=r_oo_; 
    setCorrelations();
}



// JR_Correlations

JR_Correlations::
JR_Correlations(const RealArray1D& T_, Real beta_) : 
Correlations(T_.getDimension()-1,0.0,beta_,0.0,JR), T(T_) 	
{ 
	setCorrelations(); 
}


void 
JR_Correlations::
setCorrelations()
{
    for(int i=1;i<n;i++)
    for(int j=i;j<n;j++) 
	   correlationMatrix(i,j)=exp(beta*(T[i]-T[j])); 
}


Correlations* 
JR_Correlations::
sample(int m, Real delta=0.25)
{ 
	RealArray1D T(m+1); T[0]=0.0;
	for(int i=0;i<m;i++) T[i+1]=T[i]+delta;
		
	Real beta_=0.1;
	return new JR_Correlations(T,beta_); 
}
	

ostream& 
JR_Correlations::
printSelf(ostream& os) const
{
	return
	os << "\nCorrelations: Jaeckel-Rebonato, " 
	   << "beta=" << beta;
}



// CS_Correlations
CS_Correlations::
CS_Correlations(int m, Real alpha_, Real beta_, Real r_oo_) : 
Correlations(m,alpha_,beta_,r_oo_,CS) 	
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


ostream& 
CS_Correlations::
printSelf(ostream& os) const
{
	return
	os << "\nCorrelations: Coffee-Shoenmakers, " 
	   << "alpha=" << alpha << ", beta=" << beta << ", r_oo=" << r_oo;
}


Correlations* 
CS_Correlations::
sample(int m)
{ 
	Real alpha_=1.8, beta_=0.1, r_oo_=0.4;
	return new CS_Correlations(m,alpha_,beta_,r_oo_);
}


MTGL_END_NAMESPACE(Martingale)
