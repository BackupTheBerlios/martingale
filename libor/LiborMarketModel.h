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
 * LiborMarketModel.h
 *
 * Created on April 11, 2003, 12;45 AM
 */

#ifndef martingale_libormarketmodel_h    
#define martingale_libormarketmodel_h

#include "TypedefsMacros.h"
#include "Utils.h"
#include "Matrix.h"                // problem with typedefs in forward declarations

MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file LiborMarketModel.h
 * Interface and partial implementation for all Libor market models based on
 * factor loadings {@link LiborFactorLoading} with a deterministic volatility 
 * {@link VolSurface} and constant correlations {@link Correlations}.
 */


// dependencies
struct LiborFactorLoadingType;
class LiborFactorLoading;
class Bond;  // defined below
// class RealArray1D;



/*******************************************************************************
 *
 *                     LiborMarektModelType
 * 
 ******************************************************************************/
 
 /** Object encapsulating the following information:<br> 
  *  type of simulation (driftless, low factor driftless, predictor-corrector, 
  *  fast predictor-corrector) and {@link LiborFactorLoadingType}.
  */
 struct LiborMarketModelType{
	  
	  /** DL (driftless), LFDL (low factor drifless)
	   *  PC (predictor-corrector), FPC (fast predictor-corrector).
	   */
	  static const int DL=0, LFDL=1, PC=2, FPC=3;
	 
	  /** Type of factor loading. */
	  const LiborFactorLoadingType* flType;
	 
	  /** type flag: DL, LFDL, PC, FPC. */
	  const int type;
	 	 
	  /** @param lflType factorloading type (vols and corrs).
	   *  @param lmmType LiborMarketModel::DL,LFDL,PC,FPC
	   */
	  LiborMarketModelType(const LiborFactorLoadingType* lflType, int lmmType) :
	  flType(lflType), type(lmmType)  {    }                                                                         
	  
	  std::ostream& printSelf(std::ostream&) const;
};





/*******************************************************************************
 *
 *                    LIBOR MARKET MODEL
 *
 ******************************************************************************/

/**<p>Libor Market Model based on state variables with Gaussian {@link FactorLoading}.
 * See book, Chapter 6. The state variables need not be the Libors \f$X_j\f$ themselves. 
 * They could in fact be the coterminal 
 * swap rates (all swaps terminate at the final date \f$T_n\f$, swap market model). In 
 * fact the fastest model uses the forward transported Libors
 *\f[U_j=X_j(1+X_{j+1})\dots(1+X_{n-1})\f] 
 * as state variables since these are local martingales and hence driftless in the 
 * forward martingale measure at \f$t=T_n\f$, see {@link DriftlessLMM}.</p>
 *
 * <p> If the state variables have a state dependent drift they will not follow a 
 * log-Gaussian process despite the Gaussian factor loadings. In fact the distribution of 
 * the process will most likely be unknown. This motivates the attempts to get rid of the
 * drift term ({@link DriftlessLMM}).
 * </p>
 *
 * <p>Provides interface to Libor path generation and Libors and convenience
 * implementations of a  variety of methods such as zero coupon bonds,
 * accumulation factors, swap rates and annuity numeraire. Override these if 
 * faster computations are possible in your model.</p>
 *
 * @author  Michael J. Meyer
 */
class LiborMarketModel {
	
protected:
	
	LiborMarketModelType type;
    
	// initialized from the factorloading
    int n;                            // number of forward Libors
    RealArray1D delta;                // delta[t]=T_{t+1}-T_t,length of compounding intervals
    RealArray1D T;                    // tenor structure, Tc[t]=T_t
    RealArray1D l;                    // initial Libors, l[j]=L_j(0)
    RealArray1D x;                    // initial X-Libors x[j]=X_j(0)=delta_jL_j(0)
    
    // volatility and correlation structure
    LiborFactorLoading* factorLoading;
	

public:
	
	/** Flags for the type of LMM. 
	 *  DL {@link DriftlessLMM}, 
	 *  LFDL {@link LowFactorDriftlessLMM}
	 *  PC {@link PredictorCorrectorLMM},
	 *  FPC {@link FastPredictorCorrectorLMM}.
	 */
	static const int DL=0, LFDL=1, PC=2, FPC=3;
	
	/** Type object (integer and string type flags). 
	 */
	const LiborMarketModelType* getType() const { return &type; }
	
	/** Convert integer type flag to string "DL", "LFDL", "PC", "FPC".
	 */
	static std::string lmmType(int type);
	
    /** The number n of accrual periods.
     */
    int getDimension() const { return n; }
    
    /** The array of continuous times \f$T_j\f$.
     */
    const RealArray1D& getTenorStructure() const { return T; }
    
    
    /** The array of accrual periods \f$delta_j\f$.
     */
    const RealArray1D& getDeltas() const { return delta; }
    
    /** Initial Libors \f$L_j(0)\f$.
     */
    const RealArray1D& getInitialLibors() const { return l; }
    
    
    /** The array of initial X-Libors \f$X_j(0)=\delta_jL_j(0)\f$.
     */
    const RealArray1D& getInitialXLibors() const { return x; }
    
    
    /** The {@link FactorLoading} of the Libor process
     */
    LiborFactorLoading* getFactorLoading() const { return factorLoading; }
	
		
    /** Switch to quasi random dynamics based on low discrepancy sequence.
     *  Default: nothing happens.
     */
    virtual void switchToQMC() const {  }
	
	
    /** Switches to pseudo random dynamics based on random number generator.
     *  Default: nothing happens.
     */
    virtual void switchToMC() const {  }
	
	
	/** Restarts the underlying Sobol sequence generator.
     *  Default: nothing happens.
     */
    virtual void restartSobolGenerator() const {  }



// INTITIAL LIBORS
	
   /** \f$L_i(0)\f$.*/
   Real L_i0(int i) const { return l[i]; }
   
   /** \f$X_i(0)=\delta_jL_j(0)\f$. */
   Real X_i0(int i) const { return x[i]; }


   
// CONSTRUCTION

    
    /** Constructor
	 *
     * @param fl volatility and correlation structure, see
     * {@link LiborFactorLoading}.
     * @param lmmType type of Libor market model: {@link #DL}, LFDL, PC, FPC.
     */
    LiborMarketModel(LiborFactorLoading* fl, int lmmType=DL);
    virtual ~LiborMarketModel(){ }
   
   
    /** Sample Libor market model, quarterly accrual.
	 *
	 * @param n dimension (number of Libor accrual intervals).
	 * @param lmmType type of Libor market model: {@link #DL}, LFDL, PC, FPC.
	 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
	 * @param corrType type of correlations, Correlations::CS, JR.
	 */
	static LiborMarketModel* sample(int n, int lmmType, int volType, int corrType);
 
    

// ABSTRACT METHODS
		
	
     /** X-Libor \f$X_j(T_t)\f$ from current path. 
      */
     virtual Real XL(int j, int t) const = 0;


	/** The vector \f$(X_p(T_t),...,X_{n-1}(T_t))\f$ from current path. 
      */
     virtual const RealVector& XLvect(int t, int p) = 0;
	 

// PATH COMPUTATION

     /** Write a new set of standard normal increments driving a full
	  *  Libor path to time t. Empty implementation to accommodate LMMs
	  *  that are not based on such a mechanism. 
	  *
	  *  <p>All the concrete LMMs implemented here are based on this
	  *  mechanism and forward the request to a {@link StochasticGenerator}
	  */
     virtual void newWienerIncrements(int t){   }

     /** Time step t-&gt;t+1 for Libors \f$X_j\f$ or 
	  *  \f$U_j\f$ with \f$j\geq p\f$.
	  */
	 virtual void timeStep(int t, int p) = 0;
	 
	 
	 /** A full Libor path to the horizon. Default implementation
	  *  using {@link timeStep(int,int)}.
	  */
	 virtual void newPath();
 
     
     /** Path of Libors
      *  \f[s\in[0,t]\mapsto (L_p(s),L_{p+1}(s),...,L_{n-1}(s))\f]
      *  ie. the Libors \f$L_j(s), j>=p\f$, are computed from discrete
      *  time s=0 to discrete time s=t. Default implementation
	  *  using {@link timeStep(int,int)}.
      *
      * @param t discrete time up to which Libors are computed.
      * @param p Libors evolved are \f$L_j, j=p,p+1,...,n-1\f$.
      */
     virtual void newPath(int t, int p);
	 
	 
	 	
	 /** <p>The effective dimension of the simulation, that is, the number of 
	  *  independent deviates of the distribution of the stochastic driver
	  *  (typically standard normal) needed to compute one path from discrete 
	  *  time t to discrete time s. See {@link StochasticGenerator}.
	  *
	  *  <p>Note: if this is too high (>623 for Mersenne Twister) a random number
	  *  generator may not be able to ensure equidistribution in this dimension
      *  calling into question the mathematical basis for Monte Carlo expectations.
	  *
	  *  <p>Default implementation, returns zero.
	  */ 
	 virtual int effectiveDimension(int t, int s) const { return 0; }
	 
	 
	 /** Diagnostic. Print the the standard normal increments driving the current path.
	  *  Empty implementation.
	  */
	 virtual void printWienerIncrements(int s, int t) const { }


     
     
//  LIBOR VOLATILITIES
     

     /** Annualized volatility of \f$L_i(t)\f$.
      *
      * @param i Libor index.
      */
     Real vol(int i) const;
	 
     


	
// FORWARD TRANSPORTING FACTORS
	 
	 
	 /** Libor \f$L_j(T_t)\f$.*/
	 Real L(int j, int t) const;
	 
	 
	 /** \f$H_0(0)=1/B_n(0)\f$, accrual factor \f$T_0\rightarrow T_n\f$ with Libors 
      *  in state at time t=0.
      */
     virtual Real H0() const;
	 
	 
	 	 
	 /** \f$H_i(0)=B_i(0)/B_n(0)\f$, accrual factor \f$T_i\rightarrow T_n\f$ with Libors 
      *  in state at time t=0.
      */
     virtual Real H_i0(int i) const;
     

    /** Accrual factor \f$H_i(T_t)=B_i(T_t)/B_n(T_t)\f$. This factor shifts
     *  a cashflow from time \f$T_i\rightarrow T_n\f$ with Libors in state 
     *  at time \f$T_t\f$ (discrete time t). In other words the shift is carried 
	 *  out at discrete time t. Needs all Libors \f$L_j, j\geq i\f$
     *  up to time \f$T_t\f$.
     *
     * @param t discrete time (continuous time \f$T_t\f$).
     * @param i cashflow shifted from time \f$T_i\f$ to horizon.
     */
     virtual Real H_it(int i, int t);
	 

     
    /** Accrual factor \f$H_t(T_t)=1/B_n(T_t)\f$.
	 *  Cashflow is shifted \f$T_t\rightarrow T_n\f$ with Libors in state
	 *  at time \f$T_t\f$.
     *
     * @param t cashflow shifted from time \f$T_t\f$ to horizon.
     */
     virtual Real H_ii(int t);

 
	 
              
// ZERO COUPON BONDS 


    /** The zero coupon bond \f$B_i(0)=B(0,T_i)\f$.
     *  
     * @param i bond matures at time \f$T_i\f$.
     */
     virtual Real B0(int i) const;
 
    /** The zero coupon bond \f$B_i(T_t)=B(s,T)\f$ with \f$s=T_t, T=T_i\f$. 
     *  Assumes that all Libors \f$X_k, k=t,t+1,...,n-1\f$ have been computed up 
     *  to time \f$s=T_t\f$ (discrete time t).
     *
     *  @param t bond evaluated at time \f$T_t\leq T_i\f$.
     *  @param i bond matures at time \f$T_i\f$. 
     */
    virtual Real B(int i, int t);
	

     

	 
// FORWARD SWAP RATES                      
             

   /** <p>The forward swap rate \f$S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)\f$.
    *  Needs all Libors \f$X_j, j\geq p\f$ up to time \f$T_t\f$.</p>
    *
    * <p>Straightforward implementation following the definition.
    * Extremely inefficient. Used only to test the correctness of the 
    * streamlined implementation.</p>
    *
    * @param p,q swap along \f$[T_p,T_q]\f$.
    * @param t discrete time.
    */ 
     virtual Real swRate(int p, int q, int t);
     
     
     /** The forward swap rate \f$S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)\f$.
      *  Needs all Libors \f$X_j, j\geq p\f$ up to time \f$T_t\f$.
      *
      * @param p,q swap along \f$[T_p,T_q]\f$.
      * @param t discrete time.
      */ 
     virtual Real swapRate(int p, int q, int t);


     /** The forward swap rate \f$S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)\f$ 
	  *  at time t=0.
      *
      * @param p,q swap along \f$[T_p,T_q]\f$.    
	  */ 
     virtual Real swapRate(int p, int q) const;
 
	 
	 
// ANNUITY NUMERAIRE
                  
     
    /** <p>The annuity \f$B_{pq}(t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$.
      *  Needs all Libors \f$L_j, j\geq p+1\f$ up to time t.</p>
      *
      * <p>Straightforward implementation following the definition.
      * Extremely inefficient. Used only to test the correctness of the 
      * streamlined implementation.</p>
      *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      * @param t current (discrete) time.
      */ 
     virtual Real b_pq(int p, int q, int t);
     
     
    /** <p>The annuity \f$B_{pq}(t)=\sum_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$.
      *  Needs all Libors \f$L_j, j\geq p+1\f$ up to time t.</p>
	  *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      * @param t current (discrete) time.
      */      
	 virtual Real B_pq(int p, int q, int t);

     
    /** The annuity \f$B_{pq}(t)=\sum_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$
      * at time t=0.
	  *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      */     
     virtual Real B_pq(int p, int q) const;
	 
	 
	 
	 /** <p>The forward transported annuity \f$H_{p,q}(T_t)=B_{pq}(T_t)/B_n(T_t)\f$.
      *  Needs all Libors \f$L_j, j\geq t\f$ up to time \f$T_t\f$.</p>
	  *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      * @param t current (discrete) time.
      */      
	 virtual Real H_pq(int p, int q, int t);

     
    /** The annuity \f$B_{pq}(t)=\sum_{k=p}^{q-1}\delta_kB_{k+1}(t)\f$
      * at time t=0.
	  *
      * @param p,q annuity along \f$[T_p,T_q]\f$.
      */     
     virtual Real H_pq(int p, int q) const;



	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 
   /** <p>The forecast of caplet volatility \f$\Sigma(0,T_i)\f$ to expiration.
	*  Quantity is needed for Black caplet formula. See book, 6.9.2.</p>
    *  Default implementation: error message and abort.</p>
    *
    * @param i caplet on \f$[T_i,T_{i+1}]\f$.
    */ 
     virtual Real capletAggregateVolatility(int i) const;
	 
	 
   /** <p>The forecast \f$\Sigma_{p,q}(0,T_t)\f$ of swap rate volatility to 
	*  expiration, see book 6.7.</p>
	*
    *  <p>Quantity is needed for approximate analytic swaption price.
	*  Default implementation: error message and abort.</p>
    *
    * @param p,q swap along on \f$[T_p,T_q]\f$.
    * @param t swaption exercise at time \f$T_t\leq T_p\f$.
    */ 
     virtual Real swaptionAggregateVolatility(int p, int q, int t) const;

     
	 /** <p>The volatility forecast \f$\Sigma(0,T_t)\f$ needed
	  *  to price the call on the bond \f$B\f$. See book, 6.9.2.</p>
	  *
      *  <p>Default implementation: error message and abort.
	  *  This is implemnted only in {@link DriftlessLMM}.</p>
	  *
	  * @param B the bond.
	  * @param t aggregate volatility until time \f$T_t\f$.
	  */
	 virtual Real bondAggregateVolatility(Bond* B, int t) const;
	 
	 

// STRING MESSAGE	 
	    
   /** Message and fields.*/
   virtual std::ostream& printSelf(std::ostream& os) const = 0;
    


}; // end LiborMarketModel




MTGL_END_NAMESPACE(Martingale)	

#endif


