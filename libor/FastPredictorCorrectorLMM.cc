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


/*
 * 
 *
 * Created on April 10, 2003, 9:00 PM
 */


#include "FastPredictorCorrectorLMM.h"
using namespace Martingale;



/*******************************************************************************
 *
 *                                 PredictorCorrectorLMM
 *
 ******************************************************************************/ 

// (X_p(T_t),...,X_{n-1}(T_t))\f$, 
vector<Real>& FastPredictorCorrectorLMM::XLvect(int t, int p) 
{ 
	XVec.setDimension(n-p);
	XVec.setIndexBase(p);
	for(int j=p;j<n;j++) XVec[j]=X(t,j);
	return XVec;
}


	
//  CONSTRUCTOR

    FastPredictorCorrectorLMM::
    FastPredictorCorrectorLMM(LiborFactorLoading* fl) : LiborMarketModel(fl),
    Z(n), X(n), Y(n), m0(n),
	m(new Real[n]),
	V(new Real[n]),
	G(new Real[n]),
	logLiborCovariationMatrices(n-1),
	logLiborCovariationMatrixRoots(n-1),
	SG(new MonteCarloLiborDriver(n)),
	XVec(*(new vector<Real>(n)))
	{        
        // initialize path arrays
        for(int j=0;j<n;j++){ X(0,j)=x[j]; Y(0,j)=log(x[j]); }
        
		// set pointers to matrices for time step simulation and initialize
		// the deterministic drift steps
		for(int t=0;t<n-1;t++){
			
			UTRMatrix<Real>& cvm_t=factorLoading->logLiborCovariationMatrix(t);
			UTRMatrix<Real>& cvmr_t=factorLoading->logLiborCovariationMatrixRoot(t);
			UTRMatrix<Real> *ptr_cvm_t=&cvm_t, *ptr_cvmr_t=&cvmr_t;
			logLiborCovariationMatrices.setMatrix(t,ptr_cvm_t);
            logLiborCovariationMatrixRoots.setMatrix(t,ptr_cvmr_t);
			// deterministic drift steps
			for(int j=t+1;j<n;j++){
				
				m0(t,j)=-0.5*cvm_t(j,j);
				for(int k=j+1;k<n;k++) m0(t,j)-=(x[k]/(1+x[k]))*cvm_t(j,k);
			}
		} // end for t
	} // end constructor

	

	LiborMarketModel* 
	FastPredictorCorrectorLMM::sample(int n, Real delta)
    {
		LiborFactorLoading* 
		fl=LiborMarketModel::sampleFactorLoading(n,delta,LiborMarketModel::CS);
		return new FastPredictorCorrectorLMM(fl);
	}

    

// WIENER INCREMENTS	
	

	void FastPredictorCorrectorLMM::
	printWienerIncrements(int t, int T)
    {
         // recall that Z is an upper triangular array
		 for(int s=t;s<T;s++){
			 
             for(int k=s+1;k<n;k++){ cout << Z(s,k) << " "; 
				                     if(k==n-1) cout << endl; }
		}
	} // end printWienerIncrements
	

// TIME STEP

     // time step T_t -> T_{t+1} for Libors X_j, j>=p.
     void FastPredictorCorrectorLMM::timeStep(int t, int p)
     {   
		 /* The matrices needed for the time step T_t->T_{t+1}.
          */
         UTRMatrix<Real>& C=logLiborCovariationMatrices.getMatrix(t); 
         UTRMatrix<Real>& R=logLiborCovariationMatrixRoots.getMatrix(t); 
                    
         int q=max(t+1,p);   
         // only Libors L_j, j>=q make the step. To check the index shifts 
         // to zero based array indices below consider the case p=t+1 
         // (all Libors), then q=t+1.
         
         
         // volatility step vector V
         for(int j=q;j<n;j++)
         {
             V[j]=0;
             for(int k=j;k<n;k++) V[j]+=R(j,k)*Z(t,k);  
         }
    
		 // the drift step
		 
         // compute predicted Libors using the cached deterministic drift
		 // steps m0(t,j)
         for(int j=q;j<n;j++){ Y(t+1,j)=Y(t,j)+m0(t,j)+V[j];
                               X(t+1,j)=exp(Y(t+1,j)); }
    
         // actual drift step using the predicted Libors
         for(int k=q+1;k<n;k++) G[k]=0.5*(2-1/(1+X(t,k))-1/(1+X(t+1,k)));
         for(int j=q;j<n;j++){ 
            
            m[j]=-0.5*C(j,j);
            for(int k=j+1;k<n;k++) m[j]-=G[k]*C(j,k);                
         } 

         // recompute the Libors with new drift step
         for(int j=q;j<n;j++){ Y(t+1,j)=Y(t,j)+m[j]+V[j];
                               X(t+1,j)=exp(Y(t+1,j)); }  
      
   }  // end timeStep
     


	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 

     Real FastPredictorCorrectorLMM::
     capletAggregateVolatility(int i)
     { 
         Real volsqr=factorLoading->integral_sgi_sgj_rhoij(i,i,0,Tc[i]);
		 return sqrt(volsqr);
     } 
	 

     Real FastPredictorCorrectorLMM::
	 swaptionAggregateVolatility(int p, int q, Real T)
     { 
          UTRMatrix<Real>& 
		  Q=factorLoading->logLiborCovariationMatrix(p,q,0,T).utrRoot();
		  vector<Real> x_pq(q-p,p);
		  for(int j=p;j<q;j++) x_pq[j]=(B0(j)-B0(j+1))/B_pq(p,q);
		  x_pq*=Q.transpose();
		  return x_pq.norm();
     } 
    
     
     
     
// STRING MESSAGE
    

    string FastPredictorCorrectorLMM::toString()
    {
		vector<Real> vols(n);
		for(int i=0;i<n;i++) vols[i]=vol(i); 

		ostringstream os;
		os << "\nLibor Market Model, random dynamics: " << SG->toString() << endl
		   << "fast predictor corrector simulation of true Libor dynamics." << endl
		   << factorLoading->toString() 
		   << "\n\nLibor volatilities:\n" << vols;
        
        return os.str();
     }
	 

