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


#include "StochasticProcesses.h"
#include "Random.h"
#include "StochasticProcess.h"
#include "FactorLoading.h"


MTGL_BEGIN_NAMESPACE(Martingale)


/********************************************************************************************
 *
 *                   Brownian motion
 * 
 ********************************************************************************************/


ScalarBrownianMotion::
ScalarBrownianMotion(int T, Real dt, Real x0=0.0) : 
BrownianScalarProcess(T) 
{ 
	sqrtdt=sqrt(dt); path[0]=x0; 
}
	
	
void 
ScalarBrownianMotion::
timeStep(int t){ path[t+1]=path[t]+sqrtdt*Z[t]; }
		



VectorBrownianMotion::
VectorBrownianMotion(int dim, int T, Real dt) : 
BrownianVectorProcess(dim,T) 
{  
	sqrtdt=sqrt(dt); 
}
	

VectorBrownianMotion::
VectorBrownianMotion(int dim, int T, Real dt, RealVector& x0) : 
BrownianVectorProcess(dim,T) 
{ 
	sqrtdt=sqrt(dt); 
	RealVector& X0=currentPath(0);
	for(int i=0;i<dim;i++) X0[i]=x0[i]; 
}
	

void 
VectorBrownianMotion::
timeStep(int t)
{	
	RealVector& Xt=currentPath(t);
	RealVector& XT=currentPath(t+1);
	for(int i=0;i<dim;i++) XT[i]=Xt[i]+sqrtdt*Z(t,i); 
}
	

void 
VectorBrownianMotion::
testPathFunctional(int t, int T, int nPath)
{
	BrownianVectorProcess* X=new VectorBrownianMotion(2,T,0.01);
	X->newPathSegment(t);
    RealVector& x=X->currentPath(t);

    Real analyticConditionalMean=x[0]+x[1],
         analyticVariance=2*(T-t)*0.01;
		 		 
     // Functional F(X)=X_1(T)+X_2(T)
     SumFunctional* SF=new SumFunctional(X,T);
     Real* vtMC=SF->conditionalMeanAndVariance(t,nPath,true,"Monte Carlo");
     X->switchToQMC();
     Real* vtQMC=SF->conditionalMeanAndVariance(t,nPath,true,"Quasi Monte Carlo");
		 		 
     cout << endl << endl
	      << "Paths: " << nPath << endl
	      << "Analytic conditional mean = " << analyticConditionalMean << endl
	      << "Monte Carlo conditional mean = " << vtMC[0] << endl
	      << "Quasi Monte Carlo conditional mean = " << vtQMC[0] << endl
	      << "Analytic variance = " << analyticVariance << endl
	      << "Monte Carlo variance = " << vtMC[1] << endl
	      << "Quasi Monte Carlo variance = " << vtQMC[1];
		 
} // end testPathFunctional


/********************************************************************************************
 *
 *                    Gaussian martingales
 * 
 ********************************************************************************************/

GaussianMartingale::
GaussianMartingale(int dim, int T, int ds, RealVector& x0, FactorLoading* fl) :
BrownianVectorProcess(dim,T), 
dt(ds), 
covariationMatrixRoots(T)
{
	// path intitialization
    RealVector& X0=*(path[0]);
	for(int j=0;j<dim;j++) X0[j]=x0[j];
		
    for(int t=0;t<T;t++){			
	
		const UTRRealMatrix& cv_t=factorLoading->covariationMatrix(t*dt,(t+1)*dt);
        covariationMatrixRoots[t]=&(cv_t.utrRoot());
	}
} // end constructor

GaussianMartingale::
~GaussianMartingale() 
{
	for(int t=0;t<T;t++) delete covariationMatrixRoots[t];		
} 

	
/** Time step t->t+1 based on current Z-increments.*/
void 
GaussianMartingale::
timeStep(int t)
{
    const UTRRealMatrix& R=*(covariationMatrixRoots[t]);
    RealVector& Xt=*(path[t]);
    RealVector& X=*(path[t+1]);
		
    for(int i=0;i<dim;i++){
			
	    X[i]=Xt[i];
	    for(int j=i;j<dim;j++) X[i]+=R(i,j)*Z(t,j);
    }
} // end timeStep



MTGL_END_NAMESPACE(Martingale)

			

