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




#include "LiborFactorLoading.h"
#include "LiborCalibrator.h"
#include "Utils.h"
#include "Random.h"
#include <cmath>
using namespace Martingale;
using namespace std;




/*******************************************************************************
 *
 *                     Class: LiborFactorLoading
 * 
 ******************************************************************************/


 
    LiborFactorLoading::LiborFactorLoading
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
		    cout << "\n\nLiborfactorLoading(): Libors-Correlations dimensions don't match."
		        << "\nTerminating.";
		    exit(1);
		}

	} // end constructor
	
	
	
	LiborFactorLoading* LiborFactorLoading::
	sample(int n, int volType=VolSurface::JR, int corrType=Correlation::CS)
    {
		RealArray1D deltas(n);
		RealArray1D k(n);
		RealArray1D L0(n);
        for(int i=0;i<n;i++)
		{ deltas[i]=delta; k[i]=0.2+0.1*Random::U01(); L0[i]=0.04; }
				
		VolSurface* vols = VolSurface::sample(volType);
		Correlations* corrs = Correlations::sample(n,corrType);
			
		return new LiborFactorLoading(L0,deltas,_k,vols,corrs);
	}
	
	
	
// PRINT FIELDS
	
	std::ostream& LiborFactorLoading::printSelf(std::ostream& os)
    {
		vector<Real> L0(n,l.getData());      // initial Libors		
		vector<Real> vols(n-1,1);            // annualized volatilities
		for(int i=1;i<n;i++) vols[i]=annualVol(i);
		
		return 
		os << "\n\nLiborFactorLoading: dimension = " << n 
		   << "Initial Libors: " << L0
		   << "\nAnnualized volatilities: " << vols
		   << "\nVolatility surface: " << vol
		   << "\nCorrelations: " << corr;
	}
	
	
	Real LiborFactorLoading::annualVol(int i)
    {
		Real   T_i=T[i];
		       volSqr=integral_sgi_sgj_rhoij(i,i,0.0,T_i);
		return sqrt(volSqr/T_i);
	}
					

// LOG-COVARIATION MATRICES 
	
	
   Real LiborFactorLoading::
   integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const
   {  
	   return k[i]*k[j]*rho(i,j)*(vol->integral_sgsg(s,t,T[i],T[j])); 
   }

   	
	
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int p,int q, Real s, Real t) const
   {
       int size=q-p;       // matrix size
	   UTRMatrix<Real>& logCVM=*(new UTRMatrix<Real>(size,p));

       for(int i=p;i<q;i++)
       for(int j=i;j<q;j++)
       logCVM(i,j)=integral_sgi_sgj_rhoij(i,j,s,t);

       return logCVM;
   }// end logCovariationMatrix
   

   
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int t) const
   {
	   return logLiborCovariationMatrix(t+1,n,T[t],T[t+1]);
   }
   

   
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrixRoot(int t) const
   {
	   return logLiborCovariationMatrix(t).utrRoot();
   }


   Matrix<Real>& LiborFactorLoading::
   reducedRankLogLiborCovariationMatrixRoot(int t, int r) const
   {
	   return logLiborCovariationMatrix(t).rankReducedRoot(r);
   }

   
   
   

// TEST PROGRAMS 
		 

void LiborFactorLoading::selfTest()
{
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
    cout << "\nLIBOR FACTORLOADING SELFTEST:" << endl << toString();
	printFields();
	
	cout << "\nTesting the root L of the matrix C=logLiborCovariationMatrix(t):" << endl;
	for(int t=0;t<n-1;t++){
		
		UTRMatrix<Real>& C=logLiborCovariationMatrix(t);
		UTRMatrix<Real>& L=logLiborCovariationMatrixRoot(t);
		cout << "\nt = " << t << ": ";
		C.testEquals(L.aat(),precision,epsilon,"C=L*L'");
	}

} // end test


void LiborFactorLoading::factorizationTest(int r)
{
	cerr << "\n\nApproximate factorization of all single time step log-Libor"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank " << r
	     << "\n\nRelative errors (trace norm): " << endl;
	
	for(int t=0;t<n-2;t++){
		
		UTRMatrix<Real>& Ct=logLiborCovariationMatrix(t);
		Ct.testFactorization(r);
	}
}



/*******************************************************************************
 *
 *                 VOL SURFACES
 *
 ******************************************************************************/



// OUR OWN VOL SURFACE

    // The integral \f$\int_t^T g(1-s/T1)g(1-s/T2)ds\f$. 
    // See book, 6.11.1  

    Real M_VolSurface::
	integral_sgsg(Real t, Real T1, Real T2) const     
    {
	   	Real R,f,f1,f2,D1,D2,D12;
    
        f=a*exp(-1/d); f1=f/T1; f2=f/T2;       
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
   
   
   
   
// JAECKEL-REBONATO VOLATILITY SURFACE

   Real JR_VolSurface::
   integral_sgsg(Real t, Real T1, Real T2) const
   {
	   Real f,A,B,C,
              ctmT1=c*(t-T1),
              ctmT2=c*(t-T2),
              T1mT2=fabs(T1-T2),
              q=ctmT1+ctmT2,                    // c(2t-ti-tj)
              ac=a*c,
              cd=c*d,
              bcd=b*c*d;
              
       f=exp(-Beta*T1mT2)/(c*c*c);
       A=ac*cd*(exp(ctmT2)+exp(ctmT1))+c*cd*cd*t;
       B=bcd*(exp(ctmT1)*(ctmT1-1)+exp(ctmT2)*(ctmT2-1));
       C=exp(q)*(ac*(ac+b*(1-q))+b*b*(0.5*(1-q)+ctmT1*ctmT2))/2;
       
       return f*(A-B+C);
   } // end integral_sgsg
   
   
   
   VolSurface* VolSurface::sample(int type)
   {
	   switch(type){
		   
		   case M  : return M_VolSurface::sample();
		   case JR : return JR_VolSurface::sample(); 
		   default : return CONST_VolSurface::sample(); 
	   }
   }



/*******************************************************************************
 *
 *                CORRELATIONS
 *
 ******************************************************************************/

    
    Correlations* Correlations::sample(int n, int type)
    {
	    switch(type){
		   
		    case JR : return JR_Correlations::sample(n); 
		    default : return CS_Correlations::sample(n); 
	    }
    }


   	Correlations* JR_Correlations::sample(int n, Real delta=0.25)
	{ 
		RealArray1D T(n+1); T[0]=0.0;
		for(int i=0;i<n;i++) T[i+1]=T[i]+delta;
		
		Real beta=0.1;
		return new JR_Corrleations(T,beta); 
	}
	
	
	void CS_Correlations::setCorrelations()
	{
        RealArray1D b(n-1,1); 
        for(int i=1;i<n;i++)
        { Real x=((Real)i)/(n-1); b[i]=exp(-f(x)); }
           
        // initialize the correlation matrix rho
        for(int i=1;i<n;i++)
        for(int j=i;j<n;j++) rho(i,j)=b[i]/b[j];  
	}
	
	
	Real CS_Correlations::f(Real x) const 
    {
        Real a=alpha/2, b=(beta-alpha)/6,
               c=log(r_oo)-(a+b);
              
        return x*(c+x*(a+b*x));  //cx+ax^2+bx^3
    }


