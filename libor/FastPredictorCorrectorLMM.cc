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

#include "FastPredictorCorrectorLMM.h"
#include "Array.h"
#include "Matrices.h"
#include "StochasticGenerator.h"
#include "LiborFactorLoading.h"
#include <cmath>
//#include <math.h>
#include <algorithm>
using namespace Martingale;



/*******************************************************************************
 *
 *                                 PredictorCorrectorLMM
 *
 ******************************************************************************/ 
 
 
void 
FastPredictorCorrectorLMM::
switchToQMC() 
{  
	if(SG) delete SG;
	SG = new SobolLiborDriver(n);
}
	
	
void 
FastPredictorCorrectorLMM::
switchToMC() 
{ 
	if(SG) delete SG;
	SG = new MonteCarloLiborDriver(n);
}


// LIBORS

Real 
FastPredictorCorrectorLMM::
L(int j, int t) const 
{ 
	return X(t,j)/delta[j]; 
}


Real 
FastPredictorCorrectorLMM::
XL(int j, int t) const 
{ 
	return X(t,j); 
}
 
 
// (X_p(T_t),...,X_{n-1}(T_t))\f$, 
const RealVector& 
FastPredictorCorrectorLMM::
XLvect(int t, int p) 
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
m(n), V(n), G(n),
logLiborCovariationMatrices(n-1),
logLiborCovariationMatrixRoots(n-1),
SG(new MonteCarloLiborDriver(n)),
XVec(n)
{        
    // set typeID
	type=FPC;    
	// initialize path arrays
    for(int j=0;j<n;j++){ X(0,j)=x[j]; Y(0,j)=log(x[j]); }
        
	// set pointers to matrices for time step simulation and initialize
	// the deterministic drift steps
	for(int t=0;t<n-1;t++){
			
		const UTRRealMatrix& cvm_t=factorLoading->logLiborCovariationMatrix(t);
		const UTRRealMatrix& cvmr_t=factorLoading->logLiborCovariationMatrixRoot(t);
		logLiborCovariationMatrices.setMatrix(t,cvm_t);
        logLiborCovariationMatrixRoots.setMatrix(t,cvmr_t);
		// deterministic drift steps
		for(int j=t+1;j<n;j++){
				
			m0(t,j)=-0.5*cvm_t(j,j);
			for(int k=j+1;k<n;k++) m0(t,j)-=(x[k]/(1+x[k]))*cvm_t(j,k);
		}
	} // end for t
} // end constructor

	

LiborMarketModel* 
FastPredictorCorrectorLMM::
sample(int n, int volType, int corrType)
{
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(n,volType,corrType);
	return new FastPredictorCorrectorLMM(fl);
}
 

    

// WIENER INCREMENTS	
	

void 
FastPredictorCorrectorLMM::
printWienerIncrements(int t, int s) const
{
     // recall that Z is an upper triangular array
	 for(int u=t;u<s;u++)
	 for(int k=u+1;k<n;k++){ 
			
			std::cout << Z(u,k) << " "; 
            if(k==n-1) std::cout << endl; 
	}
} // end printWienerIncrements
	


// TIME STEP

// time step T_t -> T_{t+1} for Libors X_j, j>=p.
void 
FastPredictorCorrectorLMM::
timeStep(int t, int p) 
{   
     /* The matrices needed for the time step T_t->T_{t+1}.
      */
     const UTRRealMatrix& C=logLiborCovariationMatrices.getMatrix(t); 
     const UTRRealMatrix& R=logLiborCovariationMatrixRoots.getMatrix(t); 
                    
     int q=max(t+1,p);   
     // only Libors L_j, j>=q make the step. To check the index shifts 
     // to zero based array indices below consider the case p=t+1 
     // (all Libors), then q=t+1.
         
         
     // volatility step vector V
     for(int j=q;j<n;j++){
		 
          V[j]=0;
          for(int k=j;k<n;k++) V[j]+=R(j,k)*Z(t,k);  
     }
    
     // the drift step
		 
     // compute predicted Libors using the cached deterministic drift
     // steps m0(t,j)
     for(int j=q;j<n;j++){ Y(t+1,j)=Y(t,j)+m0(t,j)+V[j];
                           X(t+1,j)=std::exp(Y(t+1,j)); }
    
     // actual drift step using the predicted Libors
     for(int k=q+1;k<n;k++) G[k]=0.5*(2-1/(1+X(t,k))-1/(1+X(t+1,k)));
     for(int j=q;j<n;j++){ 
            
         m[j]=-0.5*C(j,j);
         for(int k=j+1;k<n;k++) m[j]-=G[k]*C(j,k);                
     } 

     // recompute the Libors with new drift step
     for(int j=q;j<n;j++){ Y(t+1,j)=Y(t,j)+m[j]+V[j];
                           X(t+1,j)=std::exp(Y(t+1,j)); }  
      
}  // end timeStep



// PATHS


void 
FastPredictorCorrectorLMM::
newPath()
{
    SG->newWienerIncrements(0,n-1,Z);
	for(int t=0;t<n-1;t++)timeStep(t);
}


void 
FastPredictorCorrectorLMM::
newPath(int t, int p)
{
    SG->newWienerIncrements(0,t,Z);
    for(int s=0;s<t;s++)timeStep(s,p);
}
     


	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 
Real 
FastPredictorCorrectorLMM::
capletAggregateVolatility(int i) const
{ 
     Real volsqr=factorLoading->integral_sgi_sgj_rhoij(i,i,0,T[i]);
     return sqrt(volsqr);
} 
	 

Real 
FastPredictorCorrectorLMM::
swaptionAggregateVolatility(int p, int q, int t) const
{ 
     const UTRRealMatrix& 
     Q=factorLoading->logLiborCovariationMatrix(p,q,0,T[t]).utrRoot();
     RealVector x_pq(q-p,p);
	 for(int j=p;j<q;j++) x_pq[j]=(B0(j)-B0(j+1))/B_pq(p,q);
     x_pq*=Q.transpose();
	 return x_pq.norm();
} 
    
     
std::ostream& 
FastPredictorCorrectorLMM::
printSelf(std::ostream& os) const
{
	os << "\nLibor Market Model, fast predictor-corrector type."
	   << "\nRandom dynamics: ";
	SG->printSelf(os);
	factorLoading->printSelf(os);
	return os;
}
	 

