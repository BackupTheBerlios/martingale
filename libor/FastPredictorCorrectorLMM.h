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
 * FastPredictorCorrectorLMM.h
 *
 * Created on April 10, 2003, 9:00 PM
 */

#ifndef martingale_fastpclmm_h    
#define martingale_fastpclmm_h

#include "TypedefsMacros.h"
#include "LiborMarketModel.h"               // base class
#include "Matrix.h"                         // direct members

MTGL_BEGIN_NAMESPACE(Martingale)



// dependenies
class StochasticGenerator;       // StochasticGenerator.h
class SobolLiborDriver;
class MonteCarloLiborDriver;
class LiborFactorLoading;        // LiborFactorLoading.h
// class UTRRealMatrix;
// class RealArray1D;
// class RealVector;




/*******************************************************************************
 *
 *                                 PredictorCorrectorLMM
 *
 ******************************************************************************/ 

/**<p>Same as {@link PredictorCorrectorLMM} but with the following optimizations:
 *  A deterministic approximation to the drift step is used to compute the predicted 
 *  values. Once these are obtained the actual drift step is computed as usual.
 *  This cuts the computational effort in computing the drift step in half.
 *  All quantities involved are cached to speed up the computation.
 *
 * @author  Michael J. Meyer
 */
class FastPredictorCorrectorLMM : public LiborMarketModel {

	
	// Row t is used to drive the time step T_t->T_{t+1}
	// This is allocated as an n-dimensional matrix to preserve natural indexation
	// Row index starts at zero, first column not needed, see newWienerIncrements.
	UTRRealMatrix Z;
    
    // The following are allocated as lower triangular arrays:
    
    // Array containing the X-Libors X(t,j)=X_j(T_t)=delta_jL_j(T_t), t<=j.
    // The rows are the vectors X(T_t)=(X_t(T_t),...,X_{n-1}(T_t)).
    UTRRealMatrix X;
    
    // Array containing the log-Libors Y(t,j)=Y_j(T_t)=log(X_j(T_t)), t<=j.
    // The rows are the vectors Y(T_t)=(Y_t(T_t),...,Y_{n-1}(T_t)).
    UTRRealMatrix Y;
	
	// Array containing the deterministic drift step approximations
    // m0(t,j)=integral_{T_t}^{T_{t+1}}mu_j(s)ds, where mu_j(s) is the 
	// instantaneous drift of X_j computed approximating X_k(s)/(1+X_k(s)) 
	// with X_k(0)/(1+X_k(0)).
    UTRRealMatrix m0;

    // drift step vector
    RealArray1D m;
    
    // volatility step vector
    RealArray1D V;
    
    // vector used to store factors G_j(t)=0.5*[F_j(t)+F_j(t+1)] where 
	// F_j(t)=X_j(t)/(1+X_j(t)).
    RealArray1D G;
    
    // the log-Libor covariation matrices needed for the time steps
    UTRMatrixSequence logLiborCovariationMatrices;
                
    // the upper triangular roots of the above log-Libor covariation matrices
    UTRMatrixSequence logLiborCovariationMatrixRoots;
	
	
    StochasticGenerator* SG;    // generates the Wiener increments driving the paths     
	
	RealVector XVec;         // cache for fast returning of X-Libor vectors.



public:

	/** Switches to quasi random dynamics based on Sobol sequence.
     */
    void switchToQMC();
	
	
    /** Switches to pseudo random dynamics based on Mersenne Twister.
     */
    void switchToMC();



// LIBORS


     /** Libor \f$L_j(t)\f$, value in current path.
      *
      * @param t time.
      */
     Real L(int j, int t) const;

     
     /** X-Libor \f$X_j(T_t)=\delta_jL_j(T_t)\f$, value in current path.
      *
      * @param j Libor index.
      * @param t discrete time.
      */
     Real XL(int j, int t) const;
	 
	 
	 /** X-Libor vector \f$(X_p(T_t),...,X_{n-1}(T_t))\f$, 
	  *  value in current path. Index base p, natural indices.
      *
      * @param p index of first Libor.
      * @param t discrete time.
      */
     const RealVector& XLvect(int t, int p);
		 

	
//  CONSTRUCTOR

    
    /** Constructor, sets default stochastic generator: Mersenne Twister.
	 *
     * @param fl volatility and correlation structure, see
     * {@link FactorLoading}.
     */
    FastPredictorCorrectorLMM(LiborFactorLoading* fl);
	
	
	/** Sample LMM, quarterly accrual.
	 *
	 * @param n dimension (number of accrual intervals).
	 * @param volType type of volatility surface (VolSurface::JR,M,CONST).
	 * @param corrType type of correlations (Correlations::JR,CS).
	 */
	static LiborMarketModel* sample(int n, int volType, int corrType);

    

// WIENER INCREMENTS	
	
	
	/** Prints the matrix Z of current Wiener increments,
	 *  method is used for testing only.
	 */
	void printWienerIncrements(int t, int s) const;
	
	
	/** <p>The effective dimension of the simulation, that is, the number of 
	 *  standard normal deviates needed to compute one path. 
	 *  This is the dimension of the low discrepancy sequence generator
	 *  needed for QMC simulation.
	 *
	 *  <p>Note: if this is too high (>623 for Mersenne Twister) a random number
	 *  generator may not be able to ensure equidistribution in this dimension
	 *  calling into question the mathematical basis for Monte Carlo expectations.
	 */ 
	int effectiveDimension(int t, int s) const { return (s-t)*(2*n-(1+t+s))/2; }


		
		
// TIME STEP

    

     /** Time step \f$T_t\rightarrow T_{t+1}\f$ of the X-Libors \f$X_j, j>=p\f$
      *  (discrete time t to t+1) using a fast predictor-corrector algorithm 
	  *  and the current matrix Z of standard normal increments. See book, 6.5.1.
      * 
      * @param t current discrete time (continuous time <code>T_t</code>).
      * @param p Libors evolved are \f$X_j, j=p,p+1,\dots,n-1\f$.
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
     void newPath();
     
    
     /** Path of Libors
      *  \f[s\in[0,t]\mapsto (L_p(s),L_{p+1}(s),...,L_{n-1}(s))\f]
      *  ie. the Libors \f$L_j(s), j>=p\f$, are computed from discrete
      *  time s=0 to discrete time s=t.
      *
      * @param t discrete time up to which Libors are computed.
      * @param p Libors evolved are \f$L_j, j=p,p+1,\dots,n-1\f$.
      */
     void newPath(int t, int p);
	 
	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 
   /** <p>Analytic approximation for the aggregate caplet volatility
	*  \f[\Sigma_i(0,T_i)=\sqrt{\langle Y_i\rangle_0^{T_i}},\f] 
	*  where \f$Y_i=log(X_i)\f$. Exact if Libors have deterministic volatility. 
	*  Quantity is needed for Black caplet fromula.</p>
	*
    * @param i caplet on \f$[T_i,T_{i+1}]\f$.
    */ 
     Real capletAggregateVolatility(int i) const;
	 
	 
   /** <p>Analytic approximation for the aggregate swap rate volatility
	*  \f[\Sigma_{p,q}(0,T_t)=\sqrt{\langle Y\rangle_0^{T_t},\f] 
	*  where \f$Y=log(S_{p,q})\f$ is the logarithm of the swap rate. 
	*  Quantity is needed for Black approximate swaption formula.</p>
	*
    * @param p,q swap along on \f$[T_p,T_q]\f$.
    * @param t swaption exercise at time \f$T_t\leq T_p\f$.
    */ 
     Real swaptionAggregateVolatility(int p, int q, int t) const;
    
     
     
     
// STRING MESSAGE
    
    
    /** Message identifying the object (parameter values, etc)
     */
    std::ostream& printSelf(std::ostream& os) const;
	 

}; // end FastPredictorCorrectorLMM






MTGL_END_NAMESPACE(Martingale)	

#endif


