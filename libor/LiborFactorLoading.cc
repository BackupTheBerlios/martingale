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
	
	
	
// PRINT FIELDS
	
	std::ostream& LiborFactorLoading::printSelf(std::ostream& os) const
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
	
	
	Real LiborFactorLoading::annualVol(int i) const
    {
		if(i==0) return 0.0;
		Real   T_i=T[i],
		       volSqr=integral_sgi_sgj_rhoij(i,i,0.0,T_i);
		return sqrt(volSqr/T_i);
	}
					

// LOG-COVARIATION MATRICES 
	
	
   Real LiborFactorLoading::
   integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const
   {  
	   return k[i]*k[j]*rho(i,j)*(vol->integral_sgsg(s,t,T[i],T[j])); 
   }

   	
	
   const UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int p,int q, Real s, Real t) const
   {
       int size=q-p;       // matrix size
	   UTRMatrix<Real>& logCVM=*(new UTRMatrix<Real>(size,p));

       for(int i=p;i<q;i++)
       for(int j=i;j<q;j++)
       logCVM(i,j)=integral_sgi_sgj_rhoij(i,j,s,t);

       return logCVM;
   }// end logCovariationMatrix
   

   
   const UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int t) const
   {
	   return logLiborCovariationMatrix(t+1,n,T[t],T[t+1]);
   }
   

   
   const UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrixRoot(int t) const
   {
	   return logLiborCovariationMatrix(t).utrRoot();
   }


   const Matrix<Real>& LiborFactorLoading::
   reducedRankLogLiborCovariationMatrixRoot(int t, int r) const
   {
	   return logLiborCovariationMatrix(t).rankReducedRoot(r);
   }

   
   
   

// TEST PROGRAMS 
		 

void LiborFactorLoading::selfTest() const
{
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
    cout << "\nLIBOR FACTORLOADING SELFTEST:" << endl << this;
	
	cout << "\nTesting the root L of the matrix C=logLiborCovariationMatrix(t):" << endl;
	for(int t=0;t<n-1;t++){
		
		const UTRMatrix<Real>& C=logLiborCovariationMatrix(t);
		const UTRMatrix<Real>& L=logLiborCovariationMatrixRoot(t);
		cout << "\nt = " << t << ": ";
		C.testEquals(L.aat(),precision,epsilon,"C=L*L'");
	}

} // end test


void LiborFactorLoading::factorizationTest(int r) const
{
	cerr << "\n\nApproximate factorization of all single time step log-Libor"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank " << r
	     << "\n\nRelative errors (trace norm): " << endl;
	
	for(int t=0;t<n-2;t++){
		
		const UTRMatrix<Real>& Ct=logLiborCovariationMatrix(t);
		Ct.testFactorization(r);
	}
}


// CALIBRATION

// use the first coordinates of x to set the Libor volatility scaling factors 
// k[j] (including k[0] although its useless) then the next coordinates to set 
// VolSurface::a,b,c,d and the follwing coordinates to set 
// Correlations::alpha,beta,r_oo
void LiborFactorLoading::setParameters(const RealArray1D& X)
{
	Real* u=X.getData();
	for(int j=0;j<n;j++) k[j]=u[j];
	vol->setParameters(u[n],u[n+1],u[n+2],u[n+3]);
	corr->setParameters(u[n+4],u[n+5],u[n+6]);
}
	



/*******************************************************************************
 *
 *                 VOL SURFACES
 *
 ******************************************************************************/


	string VolSurface::volSurfaceType(int volType)
    {
         switch(volType){
			 
			 case VolSurface::CONST : return "CONST";
             case VolSurface::JR : return "JR";
			 case VolSurface::M : return "M";
		}
		return "  ";
	}



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
              q=ctmT1+ctmT2,                    // c(2t-ti-tj)
              ac=a*c,
              cd=c*d;
              
       f=1.0/(c*c*c);
       A=ac*cd*(exp(ctmT2)+exp(ctmT1))+c*cd*cd*t;
       B=b*cd*(exp(ctmT1)*(ctmT1-1)+exp(ctmT2)*(ctmT2-1));
       C=exp(q)*(ac*(ac+b*(1-q))+b*b*(0.5*(1-q)+ctmT1*ctmT2))/2;
       
       return f*(A-B+C);
   } // end integral_sgsg
   
   
   
   VolSurface* VolSurface::sample(int type)
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



/*******************************************************************************
 *
 *                CORRELATIONS
 *
 ******************************************************************************/
   
   
   	string Correlations::correlationType(int corrType)
    {
         switch(corrType){
			 
			 case Correlations::JR : return "JR";
             case Correlations::CS : return "CS";
			 default :               return "Unknown-Correlations";
		}
		return "   ";   // makes compiler happy
	}


    
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
		return new JR_Correlations(T,beta); 
	}
	
	
	void CS_Correlations::setCorrelations()
	{
        RealArray1D b(n-1,1); 
        for(int i=1;i<n;i++)
        { Real x=((Real)i)/(n-1); b[i]=exp(-f(x)); }
           
        // initialize the correlation matrix rho
        for(int i=1;i<n;i++)
        for(int j=i;j<n;j++) correlationMatrix(i,j)=b[i]/b[j];  
	}
	
	
	Real CS_Correlations::f(Real x) const 
    {
        Real a=alpha/2, b=(beta-alpha)/6,
               c=log(r_oo)-(a+b);
              
        return x*(c+x*(a+b*x));  //cx+ax^2+bx^3
    }


