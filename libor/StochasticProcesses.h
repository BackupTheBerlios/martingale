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

#include "TypedefsMacros.h"
#include "StochasticProcess.h"          // base class
#include "Array.h"

MTGL_BEGIN_NAMESPACE(Martingale)


// we are using
class FactorLoading;                     // FactorLoading.h
class UTRRealMatrix;
// class BrownianScalarProcess



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
	ScalarBrownianMotion(int T, Real dt, Real x0=0.0);
	
	/** Time step t->t+1 based on current Z-increments.*/
	void timeStep(int t);
		
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
	VectorBrownianMotion(int dim, int T, Real dt);
	
	
	/** @param dim dimension.
	 *  @param T number of time steps to horizon.
	 *  @param dt size of time step.
	 *  @param x starting point of Brownian motion.
	 */
	VectorBrownianMotion(int dim, int T, Real dt, RealVector& x0);
	
	/** Time step t->t+1 based on current Z-increments. */
	void timeStep(int t);
	
// A SMALL TEST PROGRAM
   
	 /** Allocates a two dimensional Brownian motion X with T time steps of
	  *  size dt=0.01 to the horizon. The path functional f is defined to be
	  *  f=X_1(T)+X_2(T). A path of X is computed to time t and the functional 
	  *  conditioned on the state at time t. The conditional distribution is
	  *  normal with mean X_1(t)+X_2(t) and variance 2*(T-t)*dt.
	  *  The conditional mean and variance are computed from a sample of size 
	  *  nPath and compared with the analytic values.
      */
     static void testPathFunctional(int t, int T, int nPath);
	 
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
	Array1D<const UTRRealMatrix*> covariationMatrixRoots;
	
	FactorLoading* factorLoading;   // the factor loading
	
public:
	
	/** dim dimension.
	 *  T number of time steps to horizon.
	 *  ds size of time step.
	 *  x0 state at time t=0.
	 *  fl the factor loading.
	 */
	GaussianMartingale(int dim, int T, int ds, RealVector& x0, FactorLoading* fl);

	~GaussianMartingale();

	/** Time step t->t+1 based on current Z-increments.*/
	void timeStep(int t);
			
			
}; // end GaussianMartingale




MTGL_END_NAMESPACE(Martingale)

#endif

