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
 * LognormalLMM.h
 *
 * Created on April 11, 2003, 12:00 PM
 */

#ifndef martingale_lognormallmm_h    
#define martingale_lognormallmm_h

#include "Random.h"
#include "QuasiMonteCarlo.h"
#include "LiborMarketModel.h"
#include "StochasticGenerator.h"

MTGL_BEGIN_NAMESPACE(Martingale)




/*******************************************************************************
 *
 *                                 LognormalLMM
 *
 ******************************************************************************/ 

/** <p>Lognormal Approximation to the true Libor dynamics. Proof of the futility
 *  of lognormal approximations even clever ones. <i>Disregard.</i> 
 *  Slower than {@link LightSpeedLMM} and not accurate enough to be useful. 
 *  The preferreed Libor simulation is <code>LightSpeedLMM</code>
 *  since this is the only model which can be simulated without approximation
 *  and is by far the fastest also.</p>
 *
 * <p>Paths are generated using the log-Gaussian X1-approximation to the true Libor
 * dynamics, see book, 6.5.2. 
 * Time steps move directly from one point \f$T_i\f$ on the Libor
 * tenor structure to the following point \f$T_{i+1}\f$ and are simulated in 
 * the forward martingale measure \f$P_n\f$ at the terminal time \f$t=T_n\f$.
 * The stochastic dynamics is either pseudo random (MC) based on the Mersenne Twister or
 * quasi random (QMC) based on the Sobol sequence.</p>
 *
 * <p>Log-Gaussian X1-paths track true Libor paths quite closely. Nonetheless the 
 * accumulation factors 
 * \f[f=(1+X_j(T_t))(1+X_{j+1}(T_t))...(1+X_{n-1}(T_t))\f]
 * are very unforgiving of even the smallest inaccuracies if the number n-j of 
 * compounding periods is large. A good test of this is the price of a caplet halfway 
 * out to the horizon. The payoff has to be transported forward to the horizon over 
 * a substantial number of compounding periods. The lognormal approximation can only be 
 * used if the number of periods by which a cashflow has to be compounded forward is quite 
 * small.</p>
 * 
 *
 * @author  Michael J. Meyer
 */
class LognormalLMM : public LiborMarketModel {

	
	// Row t is used to drive the time step T_t->T_{t+1}
	// This is allocated as an n-dimensional matrix to preserve natural indexation
	// Row index starts at zero, first column not needed, see newWienerIncrements.
	UTRMatrix<Real> Z;
        
    // Row t is the vector eBY(T_t)=exp(B(T_t))Y(T_t), see LMM.ps.
    UTRMatrix<Real> eBY;

    // row t is the drift step vector for the time step T_t->T_{t+1},
	// drifts(t,*)=-integral_{T_t}^{T_{t+1}}eBu(s)ds. See LMM.ps. 
	// Allocated as an n-dimensional matrix to have natural indexation.
    UTRMatrix<Real> drifts;
    
    // volatility step vector
    Real* V;
    
    // eBInverse.matrix(t)=exp(B(T_t))^{-1}, t=1,...,n-1, 
	// not initialized for t=0 (where it would be the identity matrix);
    UTRMatrixSequence eBInverse;
                
    // covariationMatrixRoots.matrix(t) is the upper triangular root of the 
	// covariation matrix Jt of the volatility increment on [T_t,T_{t+1}],
	// Jt=integral_{T_t}^{T_{t+1}}eBnu(s)eBnu(s)'ds, see LMM.ps.
    UTRMatrixSequence covariationMatrixRoots;
	
	
    StochasticGenerator* SG;    // generates the Wiener increments driving the paths      
	
	vector<Real>& XVec;         // cache for fast returning of X-Libor vectors.



public:

	/** Switches to quasi random dynamics based on Sobol sequence.
     */
    void switchToQMC() 
	{  
		if(SG) delete SG;
		SG = new SobolLiborDriver(n);
	}
	
	
    /** Switches to pseudo random dynamics based on Mersenne Twister.
     */
    void switchToMC() 
	{ 
		if(SG) delete SG;
		SG = new MonteCarloLiborDriver(n);
	}
	
    UTRMatrix<Real>& getCovariationMatrixRoot(int t) const
	{ return covariationMatrixRoots.getMatrix(t); }
	
	UTRMatrix<Real>& getDrifts() { return drifts; }


// LIBORS


     /** Libor \f$L_j(t)\f$, value in current path.
      *
      * @param t time.
      */
     Real L(int j, int t) const { return XL(j,t)/delta[j]; }

     
     /** X-Libor \f$X_j(T_t)=\delta_jL_j(T_t)\f$, value in current path.
      *
      * @param j Libor index.
      * @param t discrete time.
      */
     Real XL(int j, int t) const 
	 { 
		  UTRMatrix<Real>& eBtInv=eBInverse.getMatrix(t);
		  // component j of matrix-vector product eBtInv*eBY(t,*)
		  Real sum=0;
		  for(int k=j;k<n;k++) sum+=eBtInv(j,k)*eBY(t,k);
		  return exp(sum);
	 }
	 
	 
	 /** X-Libor vector \f$(X_p(T_t),...,X_{n-1}(T_t))\f$, 
	  *  value in current path. Index base p, natural indices.
      *
      * @param p index of first Libor.
      * @param t discrete time.
      */
     vector<Real>& XLvect(int t, int p);
		 

	
//  CONSTRUCTOR
    
    
    /** Constructor, sets default stochastic generator: Mersenne Twister.
	 *
     * @param fl volatility and correlation structure, see
     * {@link FactorLoading}.
     */
    LognormalLMM(LiborFactorLoading* fl);			
    
	
	/** Sample LMM based on {@link CS_FactorLoading}.
	 *
	 * @param n dimension (number of accrual intervals).
	 * @param delta length of each accrual interval.
	 */
	static LiborMarketModel* sample(int n, Real delta);

    

// WIENER INCREMENTS	
	
	
	/** Prints the matrix Z of current Wiener increments,
	 *  method is used for testing only.
	 */
	void printWienerIncrements(int t, int T);
	
	
	/** <p>The effective dimension of the simulation, that is, the number of 
	 *  standard normal deviates needed to compute one path. 
	 *  This is the dimension of the low discrepancy sequence generator
	 *  needed for QMC simulation.
	 *
	 *  <p>Note: if this is too high (>623 for Mersenne Twister) a random number
	 *  generator may not be able to ensure equidistribution in this dimension
	 *  calling into question the mathematical basis for Monte Carlo expectations.
	 */ 
	int effectiveDimension(int t, int T) const { return (T-t)*(2*n-(1+t+T))/2; }


		
		
// TIME STEP

    

     /** <p>Evolves the X-Libors X_j, j>=p, from time T_t to time T_{t+1}.
      *  (discrete time t to discrete time t=1) in a single step
      *  using a predictor-corrector algorithm and the current matrix Z of 
	  *  standard normal increments. See document <i>LMM.ps</i>.</p>
      * 
      * @param t current discrete time (continuous time <code>T_t</code>).
      * @param p Libors evolved are X_j, j=p,p+1,...,n-1.
      */
     void timeStep(int t, int p);
		
     
     /** <p>Evolves the full set of Libors from discrete time 
      *  <code>t</code> to time <code>t+1</code> in a single time step
      *  using the current matrix Z of standard normal increments. </p>
      *
      * @param t current discrete time (continuous time <code>T_t</code>).
      */
     void timeStep(int t){ timeStep(t,t+1); }
     
     
	 
     
	 
// PATH GENERATION
	 
     
     /** Computes a full Libor path from time zero to the horizon.
      *  Moves only the X-Libor path.
      */
     void newPath()
     {
         SG->newWienerIncrements(0,n-1,Z);
         for(int t=0;t<n-1;t++)timeStep(t);
     }
     
    
     /** Path of Libors
      *  \f[t\in[0,T]\mapsto (L_p(t),L_{p+1}(t),...,L_{n-1}(t))\f]
      *  ie. the Libors \f$L_j(t), j>=p\f$, are computed from discrete
      *  time t=0 to discrete time t=T.
      *
      * @param T discrete time up to which Libors are computed.
      * @param p Libors evolved are \f$L_j, j=p,p+1,\dots,n-1\f$.
      */
     void newPath(int T, int p)
     {
         SG->newWienerIncrements(0,T,Z);
         for(int t=0;t<T;t++)timeStep(t,p);
     }
	 
	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 
   /** <p>Analytic approximation for the aggregate caplet volatility
	*  \f[\Sigma_i(0,T_i)=\sqrt{\langle Y_i\rangle_0^{T_i}},\f] 
	*  where \f$Y_i=log(X_i)\f$. Exact if Libors have deterministic volatility. 
	*  Quantity is needed for Black caplet fromula.</p>
	*
    * @param i caplet on \f$[T_i,T_{i+1}]\f$.
    */ 
     virtual Real capletAggregateVolatility(int i);
	 
	 
   /** <p>Analytic approximation for the aggregate caplet volatility
	*  \f[\Sigma_i(0,T)=\sqrt{\langle Y_i\rangle_0^T},\f] 
	*  where \f$Y_{p,q}=log(S_{p,q})\f$ is the logarithm of the swap rate. 
	*  Quantity is needed for Black approximate swaption formula.</p>
	*
    * @param p,q swap along on \f$[T_p,T_q]\f$.
    * @param T continuous time \f$T\leq T_p\f$.
    */ 
     Real swaptionAggregateVolatility(int p, int q, Real T);
    
     
     
     
// STRING MESSAGE
    
    
    /** A message what type of factor loading it is, all the parameter values.
     */
    string toString();
	 

	 
	 
// TEST PROGRAM
	 
   
/** Sets up a sample CS_FactorLoading in dimension n and a 
 *  LognormalLMM based on that FactorLoading. Then prints
 *  the terminal n/2 Libors along 20 sample paths.
 */
static void test(int n);
	
	 
             


}; // end LognormalLMM



MTGL_END_NAMESPACE(Martingale)	

#endif


