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
 * Derivatives.h
 *
 * Created on March 21, 2003, 12:45 PM
 */
 
 
#ifndef martingale_derivatives_h    
#define martingale_derivatives_h

#include "TypedefsMacros.h"
#include "Matrices.h"
#include "RandomObject.h"
#include "ControlledRandomVariable.h"


MTGL_BEGIN_NAMESPACE(Martingale)


// we are using
class LiborMarketModel;
class Bond;



/*! \file LiborDerivatives.h
 *  <p>The class {@link Derivative} is the general interface to derivatives based 
 *  on forward prices (rather than discounted prices). The file then declares
 *  several Libor derivatives.</p>
 */
 
 // WHICH LIBORS ARE NEEDED:
 // Generally if a forward price is computed at discrete time t (continuous time T_t)
 // we need all Libors X_j, j>=t, and until time t for the forward accrual from
 // time t to the horizon n and the same is true for the U_j in the LightSpeedLMM.
 

/*******************************************************************************
 *
 *                     DERIVATIVES
 *
 ******************************************************************************/

/** General derivative. Terminology uses forward prices
 *  (as is usual for Libor derivatives) not discounted prices (as is more common
 *  for most other derivatives).
 */ 
class Derivative {
    
	     
     /* A control variate is either implemented or not, deault: false; 
	  * reset this from concrete subclass.
	  */
     bool  controlVariateIsDefined;     
       

public:								 

// ACCESSORS
	
	/** Is a control variate defined? */
	bool hasControlVariate() const { return controlVariateIsDefined; }
	
	/** Sets flag indicating a control variate to true. */
	void armControlVariate(){ controlVariateIsDefined=true; }
	

// CONSTRUCTOR


     /** Default constructor.
      */
     Derivative() : controlVariateIsDefined(false) {   }


 
	 
// FORWARD TRANSPORTED PAYOFF

     
     /** <p>Forward transported payoff computed from a new path
	  *  of the underlying. This method computes a new path of the underlying 
	  *  and then derives the payoff from this new path.</p>
      *
      *  <p>If the option has payoffs at several points along the path these 
      *  must all be transported forward and aggregated at the horizon.</p>
      *
      * <p>Why the forward transporting? We operate under the assumption
	  * that the forard martingale measure at the horizon is used for the price
	  * dynamics of the underlying and in the martingale pricing formula.
	  * This is the cleanest approach which applies to Libor derivatives and also
	  * to all other derivatives.</p>
      */
     virtual Real nextForwardPayoff() = 0;
     

     /** The forward transported payoff as a random variable.
      */
	 class ForwardPayoff : public RandomVariable {
		 
		 Derivative* LD;    // pointer to containing class
	 public:	 
		 ForwardPayoff(Derivative* LD0) : LD(LD0) {  }
		 Real nextValue();

	 }; // end ForwardPayoff
		 
		 
     /** The forward transported payoff as a random variable.
      */
     RandomVariable* forwardPayoff();
     

    
// FORWARD TRANSPORTED PAYOFF WITH CONTROL VARIATE

           
     /** Mean of the control variate.
      *  This is a default implementation (error message and program abort)
      *  intended to be overidden in concrete subclasses which do implement
      *  a control variate.
      */
     virtual Real controlVariateMean();
      
     
     /** <p>Payoff-Controlvariate pair computed from current path of the 
	  *  underlying. This is a default implementation (error message and program 
	  *  abort) intended to be overidden in concrete subclasses which do implement
      *  a control variate.</p>
      */
     virtual RealVector nextControlledForwardPayoff();
      
      
	 /** The controlled forward transported payoff as a random variable.
      */
	 class ControlledForwardPayoff : public ControlledRandomVariable {
		 
		 Derivative* LD;    // pointer to containing class
	 public:	 
		 ControlledForwardPayoff(Derivative* LD0) : LD(LD0) 
		 { 
			 // subclass must set beta coefficient
			 setBeta(); 
		 }
		 RealVector nextValue(); 
		 Real getControlVariateMean();
         
     }; // end ControlledForwardPayoff
	 

     /** The controlled forward transported payoff as a random variable.
      */
     ControlledRandomVariable* controlledForwardPayoff();
     

	 
// PRICING
	  
     
      /** The value of the forward price at time <code>t=0</code>.
       *  This is a default implementation (error message and program abort)
       *  intended to be overidden in concrete subclasses.</p>
       *
       */
      virtual Real analyticForwardPrice() const;
	 
     
      /** The value of the forward price at time <code>t=0</code>.
       *
       * @param nPath number of Libor paths simulated.
       */
      Real monteCarloForwardPrice(int nPath);
	  
	  
	  /** The value of the time forward price at time <code>t=0</code>. 
       *  Progress is reported during computation.
       *
       * @param nPath number of Libor paths simulated.
	   * @param message string descriptive of computation.
       */
      Real monteCarloForwardPrice(int nPath, string message);
             
      
      
      /** The value of the time forward price at time 
       *  <code>t=0</code> using the control variate.
       *
       * @param nPath number of Libor paths simulated.
       */
      Real controlledMonteCarloForwardPrice(int nPath);
	  
	  
	  /** Message and fields.*/
      virtual std::ostream& printSelf(std::ostream& os) const = 0;
      
     
}; // end Derivative



// GLOBAL INSERTION

std::ostream& operator << (std::ostream& os, const Derivative& fl);




/*******************************************************************************
 *
 *                          LIBOR DERIVATIVE
 *
 ******************************************************************************/

/** Class which factors out some common features of all Libor derivatives.
 */
class LiborDerivative : public Derivative {
	
protected:
	
	LiborMarketModel* LMM;  // the object generating the Libors
	int n;                  // number of Libors including L_0
	int horizon;            // time until which Libors have to be computed
	                        // to compute a payoff of the derivative

public:
	
	/** @param lmm underlying Libor market model.
	 *  @param t horizon (time until which Libors are needed to computed derivative payoff).
	 */
	LiborDerivative(LiborMarketModel* lmm, int t);
	
	/** The underlying Libor market model. */
	LiborMarketModel* getLMM() const { return LMM; }
	
	/** Resets the pointer to the underlying Libor market model. Be careful to reset 
	 *  all pointers to the LMM and all relevant parameters when implementing 
	 *  another LiborDerivative.
	 */
	virtual void setLMM(LiborMarketModel* lmm) { LMM=lmm; }
	
	/** Dimension n of the underlying Libor market model. */
	int  getDimension() const { return n; }
	
	/** Discrete time t until which Libors have to be evolved for the 
	 *  computation of the derivative payoff.
	 */
	int  getHorizon() const { return horizon; }
	
	/** The effective dimension of computing one payoff. This is the number
	 *  of independent uniform deviates generated to compute the necessary 
	 *  Libors. Default implementation returns zero.
	 */
	int effectiveDimension() const;
	
	/** Message announcing a generic Libor derivative. */
    std::ostream& printSelf(std::ostream& os) const;
	
		
		
// TESTING PRICING ROUTINES

/** Tests the analytic price against Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates from
 *  a sample of 20000 paths of the underlying Libor process.
 */
virtual void testPrice();


/** Tests the analytic price against Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates from
 *  a sample of 20000 paths of the underlying Libor process.
 *  Choice of Libor market models: PredictorCorrector,
 *  DriftlessLMM. Quarterly accrual. Changes the underlying Libor
 *  market model but does not restore the original state (model).
 *
 * @param lmmType type of LMM: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::CONST, JR, M.
 * @param corrType type of correlations: {@link Correlations}::JR, CS.
 */
void priceTest(int lmmType, int volType, int corrType);

		
}; // end LiborDerivative




/*******************************************************************************
 *
 *                          CAPLETS
 *
 ******************************************************************************/


/** Caplet \f$Cplt([T_i,T_{i+1}],\kappa)\f$ with strike rate \f$\kappa\f$ 
 *  pays off \f$h=delta_i*(L_i(T_i)-\kappa)^+=(X_i(T_i)-\kappa_i)^+\f$ at time
 *  \f$T_{i+1}\f$, where \f$\kappa_i=\delta_i\kappa\f$. It is based on a 
 *  {@link LiborMarketModel}.
 */
class Caplet : public LiborDerivative {
   
	int  i;                  // accrual interval [T_i,T_{i+1}]
	
    Real kappa,             // strike rate
	     delta_i;           // T_{i+1}-T_i
	

	     
public: 
	
    /** Resets the underlying Libor market model. 
	 */
	void setLMM(LiborMarketModel* lmm);
	

// CONSTRUCTOR

    /** Caplet on \f$[T_k,T_{k+1}]\f$ based on a Libor process 
	 *  with strike rate <code>strike</code>.
     *
     * @param k caps Libor on \f$[T_k,T_{k+1}]\f$.
     * @param strike strike rate.
	 * @param lmm underlying Libor market model.
     */
    Caplet(int k, Real strike, LiborMarketModel* lmm);
	

    /** Sample at the money caplet with i=n/3, ie. one third of the way out to the horizon.
	 *  Based on a {@link DriftlessLMM}.
	 *
	 * @param n dimension (number of Libor accrual intervals).
	 * @param int lmmType type of Libor market model: {@link LiborMarketModel::DL,PC,FPC}.
	 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
	 * @param corrType type of correlations, Correlations::CS, JR.
	 */
	static Caplet* sample(int n, int lmmType, int volType, int corrType);
		
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff computed from current state of the Libor generator and transported 
	  *  forward from time \f$T_{i+1}\f$ to time \f$T_n\f$.
      */
     Real nextForwardPayoff();

	 
     /** Mean of the control variate. This is \f$X_i(0)*H_{i+1}(0)\f$,
	  *  see book, 6.9.2.
      */
     Real controlVariateMean();
 
	 
    /** Control variate is forward transported Libor 
     *  \f$X_i(T_i)H_{i+1}(T_i)\f$, a \f$P_n\f$-martingale. See book, 6.9.2.
     */
    RealVector nextControlledForwardPayoff();
    
    
   
      
// ANALYTIC PRICE

	 
   /** Black caplet price (this is the martingale price in the LMM).
	*  Approximate price if Libor volatility is stochastic as in the
	*  {@link DriftlessLMM}.
    */
   Real analyticForwardPrice() const;
   

   
// PRINTING

	/** Message identifying the derivative. */
    std::ostream& printSelf(std::ostream& os) const;


}; // end Caplet






/*******************************************************************************
 *
 *                          SWAPTIONS
 *
 ******************************************************************************/


/** The payer swaption \f$Swpn(T,[T_p,T_q],\kappa)\f$ with strike rate \f$\kappa\f$ 
 *  exercisable at time \f$T_t\leq T_p\f$ pays off \f$h=B_{p,q}(T_t)*(S_{p,q}(T_t)-\kappa)^+$ 
 *  at time \f$T_t\f$. Here \f$B_{p,q}\f$ and \f$S_{p,q}\f$ are the annuity and swap 
 *  rate along \f$[T_p,T_q]\f$ as usual. It is based on a {@link LiborMarketModel}.
 */
class Swaption : public LiborDerivative {
   
	int  p,q,                // swap period [T_p,T_q]
	     t;                  // Swaption exercisable at time T_t
	                        
    Real kappa;           // strike rate
	

	     
public:
	
	

// CONSTRUCTOR

    /** <p>The payer swaption \f$Swpn(T,[T_p,T_q],\kappa)\f$ with strike rate 
     *  \f$\kappa\f$ exercisable at time \f$T_t\leq T_p\f$.</p>
     *
     * @param p,q period of swap \f$[T_p,T_q]\f$.
     * @param strike strike rate.
	 * @param t swaption exercises at time \f$T_t\leq T_p\f$.
	 * @param lmm underlying Libor market model.
     */
    Swaption(int _p, int _q, int _t, Real strike, LiborMarketModel* lmm);
	

    /** Sample at the money swaption with p=n/3, q=n, t=p.
	 *  Based on a {DriftlessLMM}.
	 *	 
 	 * @param n dimension (number of Libor accrual intervals).
	 * @param int lmmType type of Libor market model: {@link LiborMarketModel::DL,PC,FPC}
	 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
	 * @param corrType type of correlations, Correlations::CS, JR.
	 */
	static Swaption* sample(int n, int lmmType, int volType, int corrType);
	
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff computed from current state of the Libor generator and transported 
	  *  forward from time \f$T_t\f$ to time \f$T_n\f$.
      */
     Real nextForwardPayoff();

	 
     /** Mean of the control variate. This is \f$(B_p(0)-B_q(0))/B_n(0)\f$, see book, 6.9.2.
      */
     Real controlVariateMean();
	 
	 
    /** Control variate is \f$H_p(T_t)-H_q(T_t)\f$, a \f$P_n\f$-martingale.
	 *  See book, 6.9.2.
     */
    RealVector nextControlledForwardPayoff();
    
    
         
// ANALYTIC PRICE

	 
   /** Black caplet price (this is the martingale price in the LMM).
	*  Approximate price if Libor volatility is stochastic as in the
	*  {@link DriftlessLMM}.
    */
   Real analyticForwardPrice() const;
   


// PRINTING

	/** Message identifying the derivative. */
    std::ostream& printSelf(std::ostream& os) const;



}; // end swaption




/*******************************************************************************
 *
 *                   CALL ON BONDS
 *
 ******************************************************************************/


/** Call on a general bond \f$B(s)=\sum\nolimits_j=p^{q-1}c_jB_j(s)\f$
 *  with strike rate \f$K\f$ based on a Libor market model and exercisable at a
 *  Libor reset point $\f$T_t\f$. It is assumed that \f$c_j\geq0\f$.
 */
class BondCall : public LiborDerivative {
   
	Bond* B;                 // the underlying bond
	int p,q;                 // coupon periods [T_p,T_q]

	Real K;                  // strike rate
	int t;                   // call exercises at time T_t

	     
public:
	
	/** Resets the underlying Libor market model. Be careful to reset all
	 *  all relevant parameters when implementing another LiborDerivative.
	 */
	void setLMM(LiborMarketModel* lmm);
	

// CONSTRUCTOR

    /** <p>Call on general bond portfolio \f$D(t)=\sum\nolimits_j=p^{q-1}c_jB_j(t)\f$
     *
     * @param D the underlying bond.
	 * @apram strike the strike price.
	 * @param s call exercisable at time \f$T_s\f$.
     */
    BondCall(Bond* D, Real strike, int s);
	
	
	/** Sample bond with p=n/3, q=2*n/3, all coupons random in [0.5,1.5].
	 *  Call on this bond with strike rate = cash price of the bond
	 *  at the horizon exercisable at time \f$T_p\f$. 
	 *	 
	 * @param n dimension (number of Libor accrual intervals).
	 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
	 * @param volType type of volatility surface, {@link VolSurface}::CONST, JR, M.
	 * @param corrType type of correlations, {@link Correlations}::CS, JR.
	 */
	static BondCall* sample(int n, int lmmType, int volType, int corrType);
	
	 
	/** Sample call on zero coupon bond with i=n/2, strike price = cash price 
	 *  of the bond at the horizon exercisable at time \f$T_{i-1}\f$. 
	 *  This is a worst case for the assumptions of the analytic price formulas.
	 *	 
	 * @param n dimension (number of Libor accrual intervals).
	 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
	 * @param volType type of volatility surface, {@link VolSurface}::CONST, JR, M.
	 * @param corrType type of correlations, {@link Correlations}::CS, JR.
	 */
	static BondCall* sampleCallOnZeroCouponBond(int n, int lmmType, int volType, int corrType);

     
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff \f$h\f$ at time \f$T_t\f$ accrued forward to time \f$T_n\f$:
	  *  \f[h/B_n(T_t)=h*H_t(T_t)=\left(\sum_{j=p}^{q-1}c_jH_j(T_t)-KH_t(T_t)\right)^+\f]
	  * computed from a new Libor path.
      */
     Real nextForwardPayoff();
	 
	 
     /** Mean of the control variate, the forward price at time zero.
      */
     Real controlVariateMean();
 
	 
    /** Control variate is the forward bond price.
     */
    RealVector nextControlledForwardPayoff();
		 
   
      
// ANALYTIC PRICE

	 
   /** Black-Scholes price, assumes bond price at expiration is lognormal.
	*  See book, 6.9.4.
    */
   Real analyticForwardPrice() const;
   
   
  
// PRINTING

	/** Message identifying the derivative. */
    std::ostream& printSelf(std::ostream& os) const;

   

}; // end Bond





MTGL_END_NAMESPACE(Martingale)

#endif
