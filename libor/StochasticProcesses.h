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

#ifndef martingale_stochasticprocesses_h
#define martingale_stochasticprocesses_h

#include "Random.h"
#include "StochasticProcess.h"
#include "FactorLoading.h"

/*
 * StochasticProcesses.h
 *
 * Created on April 16, 2003, 8:00 AM
 */
 

// some examples of stochastic processes 


MTGL_BEGIN_NAMESPACE(Martingale)

/*! \file StochasticProcesses.h
 *  A collection of stochastic processes.
 */



/********************************************************************************************
 *
 *                   Brownian motion
 * 
 ********************************************************************************************/


 /** Brownian motion starting at the origin in dimension one.
  */
class ScalarBrownianMotion : public BrownianScalarProcess {
	
	Real sqrtdt; 
	
public:
	
	/** @param T number of time steps to horizon.
	 *  @param dt size of time step.
	 *  @param x0 starting point of the Brownian motion (default: the origin).
	 */
	ScalarBrownianMotion(int T, Real dt, Real x0=0.0) : BrownianScalarProcess(T) 
	{ sqrtdt=sqrt(dt); path[0]=x0; }
	
	/** Time step t->t+1 based on current Z-increments.*/
	void timeStep(int t){ path[t+1]=path[t]+sqrtdt*Z[t]; }
		
}; // end ScalarBrownianMotion



 /** Brownian motion in any dimension.
  */
class VectorBrownianMotion : public BrownianVectorProcess {
	
	Real sqrtdt; 
	
public:
	
    /** Brownian motion starts at zero.
	 *  @param dim dimension.
	 *  @param T number of time steps to horizon.
	 *  @param dt size of time step.
	 */
	VectorBrownianMotion(int dim, int T, Real dt) : 
	BrownianVectorProcess(dim,T) 
	{  sqrtdt=sqrt(dt); }
	
	
	/** @param dim dimension.
	 *  @param T number of time steps to horizon.
	 *  @param dt size of time step.
	 *  @param x starting point of Brownian motion.
	 */
	VectorBrownianMotion(int dim, int T, Real dt, vector<Real>& x0) : 
	BrownianVectorProcess(dim,T) 
	{ 
		sqrtdt=sqrt(dt); 
		vector<Real>& X0=currentPath(0);
		for(int i=0;i<dim;i++) X0[i]=x0[i]; 
	}
	
	/** Time step t->t+1 based on current Z-increments. */
	void timeStep(int t)
	{	
		vector<Real>& Xt=currentPath(t);
		vector<Real>& XT=currentPath(t+1);
		for(int i=0;i<dim;i++) XT[i]=Xt[i]+sqrtdt*Z[t][i]; 
	}
	
// A SMALL TEST PROGRAM
   
	 /** Allocates a two dimensional Brownian motion X with T time steps of
	  *  size dt=0.01 to the horizon. The path functional f is defined to be
	  *  f=X_1(T)+X_2(T). A path of X is computed to time t and the functional 
	  *  conditioned on the state at time t. The conditional distribution is
	  *  normal with mean X_1(t)+X_2(t) and variance 2*(T-t)*dt.
	  *  The conditional mean and variance are computed from a sample of size 
	  *  nPath and compared with the analytic values.
      */
     static void testPathFunctional(int t, int T, int nPath)
     {
         BrownianVectorProcess* X=new VectorBrownianMotion(2,T,0.01);
         X->newPathSegment(t);
         vector<Real>& x=X->currentPath(t);

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

}; // end VectorBrownianMotion



/********************************************************************************************
 *
 *                    Gaussian martingales
 * 
 ********************************************************************************************/

/** Driftless vectorial Ito process \f$Y(t)=\int_0^t\nu(s)dW(s)\f$ with 
 *  deterministic matrix \f$\nu(s)\f$ as in {@link FactorLoading}. 
 *  See book, 3.11. The process Y satisfies
 *  \f[dY_i(t)=\nu_i(t)dW(t)=\sigma_i(t)dV_i(t),\f]
 *  where the \f$V_i(t)\f$ are correlated Brownian motions with covariation process
 *  \f[d\langle V_i,V_j\rangle_t=\rho_{ij}dt.\f]
 */
class GaussianMartingale : public BrownianVectorProcess {
	
	// size of time step
	Real dt; 
	// covariationMatrixRoots[t] is the Cholesky root of the covariation matric C, 
	// C_ij=integral_{t*dt}^{(t+1)*dt}nu_i(s).nu_j(s)ds, i,j>=0; t=0,...,T-1.
	// this drives the time step t->t+1
	UTRMatrix<Real>** covariationMatrixRoots;
	
	FactorLoading* factorLoading;   // the factor loading
	
public:
	
	/** dim dimension.
	 *  T number of time steps to horizon.
	 *  ds size of time step.
	 *  x0 state at time t=0.
	 *  fl the factor loading.
	 */
	GaussianMartingale(int dim, int T, int ds, vector<Real>& x0, FactorLoading* fl) :
    BrownianVectorProcess(dim,T), 
    dt(ds), 
    covariationMatrixRoots(new UTRMatrix<Real>*[T])
    {
        // path intitialization
	    vector<Real>& X0=*(path[0]);
	    for(int j=0;j<dim;j++) X0[j]=x0[j];
		
	    for(int t=0;t<T;t++){			
	        UTRMatrix<Real>& cv_t=factorLoading->covariationMatrix(t*dt,(t+1)*dt);
	        covariationMatrixRoots[t]=&(cv_t.utrRoot());
	    }
    } // end constructor


	~GaussianMartingale() 
    {
		for(int t=0;t<T;t++) delete covariationMatrixRoots[t];		
		delete covariationMatrixRoots;	
	} 

	/** Time step t->t+1 based on current Z-increments.*/
	void timeStep(int t)
    {
	    UTRMatrix<Real>& R=*(covariationMatrixRoots[t]);
	    vector<Real>& Xt=*(path[t]);
	    vector<Real>& X=*(path[t+1]);
		
	    for(int i=0;i<dim;i++){
			
		    X[i]=Xt[i];
		    for(int j=i;j<dim;j++) X[i]+=R(i,j)*Z[t][j];
	    }
    } // end timeStep
			
			
}; // end GaussianMartingale




MTGL_END_NAMESPACE(Martingale)

#endif

