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

#include "DriftlessLMM.h"
#include "Utils.h"
//#include "Matrices.h"
#include <cmath>
using namespace Martingale;


/*******************************************************************************
 *
 *                                 DriftlessLMM
 *
 ******************************************************************************/ 

	 
	 
	 // (X_p(T_t),...,X_{n-1}(T_t)) 
     vector<Real>& DriftlessLMM::
	 XLvect(int t, int p) 
	 { 
		 XVec.setDimension(n-p);
		 XVec.setIndexBase(p);
		 for(int j=p;j<n-1;j++) XVec[j]=U(t,j)/H(t,j+1);
	     XVec[n-1]=U(t,n-1);
		 return XVec;
	 }
		 

	
//  CONSTRUCTOR

 
    DriftlessLMM::
	DriftlessLMM(LiborFactorLoading* fl) : LiborMarketModel(fl),
    Z(n), U(n), Y(n), H(n+1), m(n), 
	V(new Real[n]),
	logLiborCovariationMatrices(n-1),
	logLiborCovariationMatrixRoots(n-1),
	SG(new MonteCarloLiborDriver(n)),
	XVec(*(new vector<Real>(n)))
	{        
        // initialize U,Y path arrays
        for(int j=0;j<n;j++){ 
			
			U(0,j)=x[j]; for(int k=j+1;k<n;k++) U(0,j)*=1+x[k];
			Y(0,j)=log(U(0,j)); 
		}
        
		// set pointers to matrices for time step simulation and initialize
		// the deterministic drift steps
		for(int t=0;t<n-1;t++){
			
			UTRMatrix<Real>& cvm_t=factorLoading->logLiborCovariationMatrix(t);
			UTRMatrix<Real>& cvmr_t=factorLoading->logLiborCovariationMatrixRoot(t);
			UTRMatrix<Real> *ptr_cvm_t=&cvm_t, *ptr_cvmr_t=&cvmr_t;
			logLiborCovariationMatrices.setMatrix(t,ptr_cvm_t);
            logLiborCovariationMatrixRoots.setMatrix(t,ptr_cvmr_t);
			// deterministic drift steps
			for(int j=t+1;j<n;j++) m(t,j)=-0.5*cvm_t(j,j);
			
		} // end for t
		
		// accrual factors at time zero
        Real f=1;
	    for(int j=n-1;j>=0;j--){ f+=U(0,j); H(0,j)=f; }
		// final accrual factor T_n->T_n
		for(int t=0;t<n;t++)H(t,n)=1;
		
	} // end constructor
	
	
		
	LiborMarketModel* DriftlessLMM::
	sample(int n, Real delta)
    {
		LiborFactorLoading* 
		fl=LiborMarketModel::sampleFactorLoading(n,delta,LiborMarketModel::CS);
		return new DriftlessLMM(fl);
	}


    

// WIENER INCREMENTS	
	

	void DriftlessLMM::printWienerIncrements(int t, int T)
    {
         // recall that Z is an upper triangular array
		 for(int s=t;s<T;s++){
			 
             for(int k=s+1;k<n;k++){ cout << Z(s,k) << " "; 
				                     if(k==n-1) cout << endl; }
		}
	} // end printWienerIncrements
	

		
		
// TIME STEP

    

     // Time step T_t->T_{t+1} for Libors X_j, j>=p
     void DriftlessLMM::timeStep(int t, int p)
     {   
		 /* The matrices needed for the time step T_t->T_{t+1}.
          */
         UTRMatrix<Real>& R=logLiborCovariationMatrixRoots.getMatrix(t); 
                    
         int q=max(t+1,p);   
         // only Libors U_j, j>=q make the step. To check the index shifts 
         // to zero based array indices below consider the case p=t+1 
         // (all Libors), then q=t+1.
         
         
         // volatility step vector V
         for(int j=q;j<n;j++)
         {
             V[j]=0;
             for(int k=j;k<n;k++) V[j]+=R(j,k)*Z(t,k);  
         }
    
		 // the drift step
		 
         // compute Y_j=log(U_j) using the cached deterministic drift
		 // steps m(t,j)
         for(int j=q;j<n;j++){ Y(t+1,j)=Y(t,j)+m(t,j)+V[j];
                               U(t+1,j)=exp(Y(t+1,j)); }
							  
		 // write the accrual factors H_j=B_j/B_n
         Real f=1;
		 for(int j=n-1;j>=q;j--){ f+=U(t+1,j); H(t+1,j)=f; }
 
      
   }  // end timeStep
     

   
     // H_pq(T_t)
     Real DriftlessLMM::H_pq(int p, int q, int t)
     {
		 Real sum=0;
		 for(int k=p;k<q;k++) sum+=delta[k]*H_it(k+1,t);
		 return sum;
	 }


// FORWARD SWAP RATES                      
             

     
     // S_{pq}(T_t)=k(T_t,[T_p,T_q])
     Real DriftlessLMM::swapRate(int p, int q, int t)
     { 
         Real num=U(t,p), denom=delta[p]*H(t,p+1);
		 for(int j=p+1;j<q;j++){ num+=U(t,j); denom+=delta[j]*H(t,j+1); }
		 return num/denom;
     } //end Swap_Rate


	 
	 
// ANNUITY NUMERAIRE
                  
     
     
     // B_{pq}(T_t)=\sum_{k=p}^{q-1}\delta_kB_{k+1}(T_t).
     Real DriftlessLMM::B_pq(int p, int q, int t)
     { 
         Real S=delta[p]*H(t,p+1);
		 for(int j=p+1;j<q;j++) S+=delta[j]*H(t,j+1);
	     return S/H(t,t);
     } //end B_pq

 
	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 

     Real DriftlessLMM::
	 capletAggregateVolatility(int i)
     { 
		 Real Ti=Tc[i];
		 UTRMatrix<Real>& R=factorLoading->logLiborCovariationMatrix(i,n,0,Ti).utrRoot();
		 vector<Real> x(n-i,i);
		 x[i]=1;
		 Real f=H(0,i+1); 
		 for(int j=i+1;j<n;j++)x[j]=-U(0,j)/f;
		 x*=R.transpose();
		 return x.norm();
     } // end capletAggregateVolatility
	 
	 

     Real DriftlessLMM::
	 swaptionAggregateVolatility(int p, int q, int t)
     { 
         Real Tt=Tc[t];
		 UTRMatrix<Real>& Q=factorLoading->logLiborCovariationMatrix(p,n,0,Tt).utrRoot();
		 vector<Real> x(n-p,p);
		 Real denom1=H(0,p)-H(0,q),
		      denom2=0;
		 for(int j=p;j<q;j++) denom2+=delta[j]*H(0,j+1);
			 
		 x[p]=U(0,p)/denom1;
		 for(int j=p+1;j<q;j++)x[j]=U(0,j)/denom1-(Tc[j]-Tc[p])*U(0,j)/denom2;
		 for(int j=q;j<n;j++)  x[j]=-(Tc[q]-Tc[p])*U(0,j)/denom2;
		 x*=Q.transpose();
		 return x.norm();
     } // end swaptionAggregateVolatility



     Real DriftlessLMM::
	 bondAggregateVolatility(Bond* B, int t)
     { 
         int p=B->get_p(), 
		     q=B->get_q();
		 Real F=B->forwardPrice();
		 vector<Real>& b=B->get_b();
		 		 
		 Real Tt=Tc[t];
		 UTRMatrix<Real>& R=factorLoading->logLiborCovariationMatrix(t,n,0,Tt).utrRoot();
		 vector<Real> x(n-t,t);
			 
		 for(int j=t;j<p;j++)  x[j]=-U(0,j)/H_i0(t);
		 for(int j=p;j<q;j++)  x[j]=b[j]*U(0,j)/F-U(0,j)/H_i0(t);
		 for(int j=q;j<n;j++)  x[j]=b[q-1]*U(0,j)/F-U(0,j)/H_i0(t);
		 x*=R.transpose();
		 return x.norm();
     } // end bondAggregateVolatility


	 


 	 
// STRING MESSAGE
    

    string DriftlessLMM::toString()
    {
		vector<Real> vols(n); vols[0]=0;
		for(int i=1;i<n;i++) vols[i]=vol(i); 

		ostringstream os;
		os << "\nDriftless Libor Market Model, random dynamics: " << SG->toString() << endl
		   << "State variables: Gaussian forward transported Libors" << endl 
		   << "U_j=X_j(1+X_{j+1})...(1+X_{n-1})" << endl
		   << "These are driftless, very fast exact simulation." << endl
		   << factorLoading->toString() 
		   << "\n\nLibor volatilities:\n" << vols;      
        return os.str();
     }
	 
             



