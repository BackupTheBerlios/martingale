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


#ifndef martingale_lowfactordriftlesslmm_h    
#define martingale_lowfactordriftlesslmm_h

#include "TypedefsMacros.h"
#include "LiborMarketModel.h"
//#include "Matrices.h"
#include "Matrix.h"

MTGL_BEGIN_NAMESPACE(Martingale)



// forward declarations
class UTRMatrixSequence;
class StochasticGenerator;       // StochasticGenerator.h
class SobolLiborDriver;
class MonteCarloLiborDriver;
class LiborFactorLoading;        // LiborFactorLoading.h



/*******************************************************************************
 *
 *                                 LowFactorDriftlessLMM
 *
 ******************************************************************************/ 

/**<p>A driftless Libor market model {@link DriftlessLMM} where
 * the covariation matrices C(t) which control the time steps \f$T_t\to T_{t+1}\f$
 * are approximated as products \f$C(t)\simeq R(t)R(t)'\f$ with R(t) of rank r. The idea
 * is to speed up the simulation by reducing matrix sizes.
 *
 * @author  Michael J. Meyer
 */
class LowFactorDriftlessLMM : public LiborMarketModel {
	
	
	int r;    // number of factors (rank to which covariation matrices are reduced).

	
	// n-1 by r matrix of standard normal increments.
	// Row t is used to drive the time step T_t->T_{t+1}
	RealMatrix Z;
    
    // Array containing the X-Libors X(t,j)=X_j(T_t)=delta_jL_j(T_t), t<=j.
    // The rows are the vectors X(T_t)=(X_t(T_t),...,X_{n-1}(T_t)).
    UTRRealMatrix U;
    
    // Array containing the log-Libors Y(t,j)=Y_j(T_t)=log(X_j(T_t)), t<=j.
    // The rows are the vectors Y(T_t)=(Y_t(T_t),...,Y_{n-1}(T_t)).
    UTRRealMatrix Y;
	
	// accrual factors H(t,j)=H_j(T_t)=B_j(T_t)/B_n(T_t), t<=j.
	// we take this up to j=n where H(t,n)=1 to avoid messy code.
	UTRRealMatrix H;         
 
	
	// Array containing the deterministic drift steps
    // m(t,j)=-0.5*integral_{T_t}^{T_{t+1}}sigma^2_j(s)ds of the log(U_j).
    UTRRealMatrix m;
    
    // volatility step vector
    RealArray1D V;
    
 
    // the log(U_j) covariation matrices of the needed for the time steps
    UTRMatrixSequence logLiborCovariationMatrices;
                
    // the upper triangular roots of the above covariation matrices
    MatrixSequence lowRankCovariationMatrixRoots;
	
	
    StochasticGenerator* SG;    // generates the Wiener increments driving the paths     
	
	RealVector XVec;          // cache for fast returning of X-Libor vectors.
  

public:

	/** Switches to quasi random dynamics based on Sobol sequence.
     */
    void switchToQMC();
	
	
    /** Switches to pseudo random dynamics based on Mersenne Twister.
     */
    void switchToMC();

	
// FACTOR ANALYSIS
	
	/** Prints the relative error in the trace norm of the rank r approximation
	 *  of the covariation matrices C(t) which control the time steps.
	 *  See {@link LiborFactorLoading#factorizationTest}.
	 */
   void factorizationTest() const;


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
	 *  Covariation matrices controlling the time steps approximated to rank r.
	 *
     * @param fl volatility and correlation structure, see {@link LiborFactorLoading}.
	 * @param r number of factors.
     */
    LowFactorDriftlessLMM(LiborFactorLoading* fl, int r);
	
	
	/** Sample LMM, quarterly accrual.
	 *
	 * @param n dimension (number of accrual intervals).
	 * @param volType type of volatility surface (VolSurface::JR,M,CONST).
	 * @param corrType type of correlations (Correlations::JR,CS).
	 * @param r number of factors.
	 */
	static LiborMarketModel* sample(int n, int volType, int corrType, int r=3);
    

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
	int effectiveDimension(int t, int s) const { return (s-t)*r; }


		
		
// TIME STEP

    

     /** <p>Time step \f$T_t\rightarrow T_{t+1}\f$ of X-Libors \f$X_j, j>=p\f$.
      *  (discrete time t to discrete time t+1) in a single step
      *  using a predictor-corrector algorithm and the current matrix Z of 
	  *  standard normal increments. See book 6.5.1.</p>
      * 
      * @param t current discrete time (continuous time \f$T_t\f$).
      * @param p Libors evolved are \f$X_j, j=p,p+1,...,n-1\f$.
      */
     void timeStep(int t, int p);
		
     
     /** Evolves the full set of Libors from discrete time t
      *  to time t+1 in a single time step
      *  using the current matrix Z of standard normal increments.
      *
      * @param t current discrete time (continuous time <code>T_t</code>).
      */
     void timeStep(int t){ timeStep(t,t+1); }
     
     
	 
	 
// PATH GENERATION
	 
     
     /** Computes a full Libor path from time zero to the horizon.
      */
     void newPath();
     
    
     /** Path of Libors
      *  \f[s\in[0,t]\mapsto (L_p(s),L_{p+1}(s),...,L_{n-1}(s))\f]
      *  ie. the Libors \f$L_j(s), j>=p\f$, are computed from discrete
      *  time s=0 to discrete time s=t.
      *
      * @param t discrete time up to which Libors are computed.
      * @param p Libors evolved are \f$L_j, j=p,p+1,...,n-1\f$.
      */
     void newPath(int t, int p);
	 
	 

//  OVERRIDE lIBOR VOLATLTIES    
//  FASTER VERSIONS OF BONDS, SWAPRATES,...
     

     /** Annualized volatility of \f$L_i\f$.
      *
      * @param i Libor index.
      */
     Real vol(int i) const;
	 
     


// FORWARD TRANSPORTING FACTORS
	 
	 
	 	 
	 /** Accrual factor \f$H_i(0)=B_i(0)/B_n(0)\f$, forward transport  
      *  \f$T_i\rightarrow T_n\f$ with Libors in state at time t=0.
      */
     Real H_i0(int i) const;
     

    /** <p>Accrual factor \f$H_i(T_t)=B_i(T_t)/B_n(T_t)\f$. This factor shifts
     *  a cashflow from time \f$T_i\rightarrow T_n\f$
     *  with Libors in state at time \f$T_t\f$ (discrete time t).
     *  In other words the shift is carried out at discrete time t.
	 *  Needs all \f$U_j, j\geq i\f$ up to time \f$T_t\f$.
     *
     * @param t current discrete time (continuous time <code>T_t</code>).
     * @param i cashflow shifted from time <code>T_i</code> to horizon.
     */
     Real H_it(int i, int t);


     
    /** \f$H_i(T_i)\f$, that is the shift \f$T_i\rightarrow T_n\f$ 
	 *  is carried out with Libors in state at time \f$T_i\f$ (discrete time t=i).
     *
     * @param i cashflow shifted from time \f$T_i\f$ to horizon.
     */
     Real H_ii(int i) const;
	 
	 
	/** \f$H_{p,q}(T_t)=B_{p,q}(T_t)/B_n(T_t)\f$, the forward price at horizon of the 
	 *  annuity along \f$[T_p,T_q]\f$ at time \f$T_t\f$.
     *
     * @param p,q annuity along \f$[T_p,T_q]\f$.
	 * @param t price at time \f$T_t\f$ (accrued forward to time \f$T_n\f$).
     */
     Real H_pq(int p, int q, int t);


 
	 
              
// ZERO COUPON BONDS 


    /** The zero coupon bond \f$B_i(0)=B(0,T_i)\f$.
     *  
     * @param i bond matures at time \f$T_i\f$.
     */
     Real B0(int i) const;

 
    /** The zero coupon bond \f$B_i(T_t)=B(s,T)\f$ with \f$s=T_t, T=T_i\f$. 
     *  Assumes that all Libors \f$U_k, k=t,t+1,...,n-1\f$ 
     *  have been computed up to time \f$s=T_t\f$ (discrete time t).
     *
     *  @param t bond evaluated at time \f$T_t\leq T_i\f$
     *  @param i bond matures at time \f$T_i\f$.
     */
    Real B(int i, int t);

	

     

	 
// FORWARD SWAP RATES                      
             

     
     /** <p>The forward swap rate \f$S_{pq}(t)=k(t,[T_p,T_q])\f$ at discrete 
      *  time t (continuous time \f$T_t\f$). 
      *  Needs all Libors \f$U_j, j\geq p\f$ up to time t.</p>
      *
      * @param p,q swap along \f$[T_p,T_q]\f$.
      * @param t discrete time.
      */ 
     Real swapRate(int p, int q, int t);


     /** <p>The forward swap rate \f$S_{pq}(t)=k(t,[T_p,T_q])\f$ at time t=0.
      *
      * @param p,q swap along \f$[T_p,T_q]\f$.
      */ 
     Real swapRate(int p, int q) { return swapRate(p,q,0); } 

 
	 
	 
// ANNUITY NUMERAIRE
                  
     
     
     /** The annuity \f$B_{pq}(t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$.
      *  Needs all Libors \f$U_j, j\geq p+1\f$ up to time t.
      *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      * @param t discrete time.
      */ 
     Real B_pq(int p, int q, int t);

     
     /** The annuity \f$B_{pq}(t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$
      *  at time t=0.
      *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      */ 
     Real B_pq(int p, int q) { return B_pq(p,q,0); }
	 
	 
	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 
   /** <p>Forecast for the total caplet volatility 
	*  \f[\Sigma_i(0,T_i)=\sqrt{\langle Y_i\rangle_0^{T_i}},\f]
	*  on the interval \f$[0,T_i]\f$. 
	*  Here \f$Y_i=log(X_i)\f$. Exact if Libors have deterministic volatility. 
	*  Quantity is needed for Black caplet formula.</p>
	*
    * @param i caplet on \f$[T_i,T_{i+1}]\f$.
    */ 
     Real capletAggregateVolatility(int i) const;
	 
	 
   /** Forecast for the total swap rate volatility
	*  \f[\Sigma_{p,q}(0,T_t)=\sqrt{\langle Y\rangle_0^{T_t}},\f]
	*  on the interval \f$[0,T_t]\f$. 
	*  Here \f$Y=log(S_{p,q})\f$ is the logarithm of the swap rate. 
	*  Quantity is needed for Black approximate swaption formula.
    *
    * @param p,q swap along on \f$[T_p,T_q]\f$.
    * @param t swaption exercise at time \f$T_t\leq T_p\f$.
    */ 
     Real swaptionAggregateVolatility(int p, int q, int t) const;


   /** Forecast for the total volatility \f$\Sigma(0,T_t)\f$ (of the
	*  relevant quantity) on the interval \f$[0,T_t]\f$ 
    *  needed for pricing the call on a bond. See book, 6.9.4.
    *
    * @param B the bond
    * @param t call on bond exercises at time \f$T_t\f$.
    */ 
     Real bondAggregateVolatility(Bond* B, int t) const;


	 
// STRING MESSAGE
    
    
    /** A message what type of factor loading it is, all the parameter values.
     */
    std::ostream& printSelf(std::ostream& os) const;
	 
	 
             


}; // end DriftlessLMM



MTGL_END_NAMESPACE(Martingale)	

#endif


