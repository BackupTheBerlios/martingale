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

#include <string>
#include <sstream>
#include <iostream>
#include "FinMath.h"
#include "RandomObject.h"
#include "ControlledRandomVariable.h"
#include "PredictorCorrectorLMM.h"
#include "FastPredictorCorrectorLMM.h"
#include "LognormalLMM.h"
#include "DriftlessLMM.h"
#include "LowFactorDriftlessLMM.h"

MTGL_BEGIN_NAMESPACE(Martingale)

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
     virtual Real nextForwardPayoff() const =0;
     

     /** The forward transported payoff as a random variable.
      */
	 class ForwardPayoff : public RandomVariable {
		 
		 Derivative* LD;    // pointer to containing class
	 public:	 
		 ForwardPayoff(Derivative* LD0) : LD(LD0) {  }
		 Real nextValue() { return LD->nextForwardPayoff(); }

	 }; // end ForwardPayoff
		 
		 
     /** The forward transported payoff as a random variable.
      */
     RandomVariable* forwardPayoff() 
     {
         return new ForwardPayoff(this);
     } 
     

    
// FORWARD TRANSPORTED PAYOFF WITH CONTROL VARIATE

           
     /** Mean of the control variate.
      *  This is a default implementation (error message and program abort)
      *  intended to be overidden in concrete subclasses which do implement
      *  a control variate.
      */
     virtual Real controlVariateMean() const
     {
         cout << "Derivative.controlVariateMean():"
		      << endl << "no control variate implemented, aborting.";
         exit(0);
         return 0.0;    // keeps the compiler happy
     } 
      
     
     /** <p>Payoff-Controlvariate pair computed from current path of the 
	  *  underlying. This is a default implementation (error message and program 
	  *  abort) intended to be overidden in concrete subclasses which do implement
      *  a control variate.</p>
      */
     virtual vector<Real> nextControlledForwardPayoff() const
     {
         cout << "Derivative.nextControlledForwardPayoff():"
		      << endl << "no control variate implemented, aborting.";
         exit(0);
		 // keep the compiler happy
         vector<Real> v(2); v[0]=0.0; v[1]=0.0;
         return v;    
     } 
      
      
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
		 vector<Real> nextValue() { return LD->nextControlledForwardPayoff(); }   
		 Real getControlVariateMean() { return LD->controlVariateMean(); }
         
     }; // end ControlledForwardPayoff
	 

     /** The controlled forward transported payoff as a random variable.
      */
     ControlledRandomVariable* controlledForwardPayoff() 
     { return new ControlledForwardPayoff(this); } 

     

	 
// PRICING
	  
     
      /** The value of the forward price at time <code>t=0</code>.
       *  This is a default implementation (error message and program abort)
       *  intended to be overidden in concrete subclasses.</p>
       *
       */
      virtual Real analyticForwardPrice() const
      {
         cout << endl 
              << "Derivative.analyticForwardPrice():" << endl
              << "no analytic price implemented, aborting.";
         exit(0);
         return 0.0;    // keeps the compiler happy
      }
      
     
      /** The value of the forward price at time <code>t=0</code>.
       *
       * @param nPath number of Libor paths simulated.
       */
      Real monteCarloForwardPrice(int nPath) 
      {
          return forwardPayoff()->expectation(nPath);
      }
	  
	  
	  /** The value of the time forward price at time <code>t=0</code>. 
       *  Progress is reported during computation.
       *
       * @param nPath number of Libor paths simulated.
	   * @param message string descriptive of computation.
       */
      Real monteCarloForwardPrice(int nPath, string message) 
      {
          return forwardPayoff()->expectation(nPath,message);
      }
             
      
      
      /** The value of the time forward price at time 
       *  <code>t=0</code> using the control variate.
       *
       * @param nPath number of Libor paths simulated.
       */
      Real controlledMonteCarloForwardPrice(int nPath) 
      {
          return controlledForwardPayoff()->expectation(nPath);
      }
      
     
}; // end Derivative




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
	LiborDerivative(LiborMarketModel* lmm, int t) : 
	LMM(lmm), n(lmm->getDimension()), horizon(t) 
	{ }
	
	/** The underlying Libor market model. */
	LiborMarketModel* getLMM(){ return LMM; }
	
	/** Resets the pointer to the underlying Libor market model. Be careful to reset 
	 *  all pointers to the LMM and all relevant parameters when implementing 
	 *  another LiborDerivative.
	 */
	virtual void setLMM(LiborMarketModel* lmm) { LMM=lmm; }
	
	/** Dimension n of the underlying Libor market model. */
	int  getDimension(){ return n; }
	
	/** Discrete time t until which Libors have to be evolved for the 
	 *  computation of the derivative payoff.
	 */
	int  getHorizon(){ return horizon; }
	
	/** The effective dimension of computing one payoff. This is the number
	 *  of independent uniform deviates generated to compute the necessary 
	 *  Libors. Default implementation returns zero.
	 */
	int effectiveDimension(){ return LMM->effectiveDimension(0, horizon); }
	
	/** Message announcing a generic Libor derivative. */
	virtual string toString(){ return "Generic LiborDerivative"; }
	
		
		
	
// TESTING PRICING ROUTINES

/** Tests the analytic price against Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates from
 *  a sample of 20000 paths of the underlying Libor process.
 */
virtual void testPrice()
{ 
	 // all prices forward prices at time T_n
     Real aprice,               // analytic price
          mcprice,              // Monte carlo price
          cvmcprice,            // Monte carlo price with control variate
	      epsilon=0.00000001;   // replacement for zero denominator in relError
         
	 // correlation with control variate
     cout << toString() << endl
	      << LMM->toString() << endl
	      << endl << "Effective dimension of the simulation = " 
	      << effectiveDimension()
	      << endl;
		
	  int nPath=20000;
      // prices
      aprice=analyticForwardPrice();
      mcprice=monteCarloForwardPrice(nPath);
	  cvmcprice=controlledMonteCarloForwardPrice(nPath);

      cout << "\nForward prices, " << nPath << " paths:" << endl
		   << "Analytic: " << aprice << endl
		   << "Monte Carlo: " << mcprice << " , relative error: " 
		   << relativeError(aprice,mcprice,epsilon) << "%" << endl
		   << "Monte Carlo with control variate: " << cvmcprice << ", relative error: " 
		   << relativeError(aprice,cvmcprice,epsilon) << "%"
		   << endl;
		   

      // QMC dynamics
	  LMM->switchToQMC();
		 
      // prices
      mcprice=monteCarloForwardPrice(nPath);
	  LMM->restartSobolGenerator();
	  cvmcprice=controlledMonteCarloForwardPrice(nPath);
 
      cout << "Quasi Monte Carlo, relative error: " 
		   << relativeError(aprice,mcprice,epsilon) << "%" << endl
		   << "Quasi Monte Carlo with control variate, relative error: " 
		   << relativeError(aprice,cvmcprice,epsilon) << "%"
		   << endl;

       
} // end testPrice


/** Tests the analytic price against Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates from
 *  a sample of 20000 paths of the underlying Libor process.
 *  Choice of Libor market models: PredictorCorrector,
 *  FastPredictorCorrector, DriftlessLMM. Changes the underlying Libor
 *  market model but does not restore the original state (model).
 *
 * @param delta length of each accrual period (will be reset).
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link DriftlessLMM} Libor market model.
 */
void priceTest
(Real delta, bool LS=true, bool PC=false, bool FPC=false) 
{
     Timer watch;

	 if(LS){
		 watch.start();
		 setLMM(DriftlessLMM::sample(n,delta)); 
         testPrice();
	     watch.stop();
	     watch.report("DriftlessLMM");
	} // end if
	 
	 if(PC){
		 watch.start();
	     setLMM(PredictorCorrectorLMM::sample(n,delta));   
         testPrice();
    	 watch.stop();
	     watch.report("PredictorCorrectorLMM");
	 } // end if
	
	 if(FPC){
		 watch.start();
	     setLMM(FastPredictorCorrectorLMM::sample(n,delta));   
         testPrice();
    	 watch.stop();
	     watch.report("FastPredictorCorrectorLMM");
	 } // end if

} // end testPrice
	

		
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
	void setLMM(LiborMarketModel* lmm)
	{ 
		LiborDerivative::setLMM(lmm); 
		// we must now adjust delta_i
		delta_i=(lmm->getDeltas())[i];
	}
	

// CONSTRUCTOR

    /** Caplet on \f$[T_k,T_{k+1}]\f$ based on a Libor process 
	 *  with strike rate <code>strike</code>.
     *
     * @param k caps Libor on \f$[T_k,T_{k+1}]\f$.
     * @param strike strike rate.
	 * @param lmm underlying Libor market model.
     */
    Caplet(int k, Real strike, LiborMarketModel* lmm) :
    LiborDerivative(lmm,min<int>(k+1,lmm->getDimension()-1)),               // Libors needed until time t=min(k+1,n-1)
	i(k), 
	kappa(strike),
	delta_i((lmm->getDeltas())[i])
    {   
		armControlVariate();
	} 

    /** Sample at the money caplet with i=n/3, ie. one third of the way out to the horizon.
	 *  Based on a {@link DriftlessLMM} with a {@link CS_FactorLoading}.
	 *  Temporarily based on a {@link PredictorCorrectorLMM} with a {@link JR_FactorLoading},
	 *  to study failure of the predictor-corrector algorithm.
	 */
	static Caplet* sample(int n, Real delta)
    {
		LiborMarketModel* lmm=DriftlessLMM::sample(n,delta);
		int i=n/3;
		Real kappa=(lmm->getInitialLibors())[i];
		return new Caplet(i,kappa,lmm);
	}
		
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff computed from current state of the Libor generator and transported 
	  *  forward from time \f$T_{i+1}\f$ to time \f$T_n\f$.
      */
     Real nextForwardPayoff() const 
	 { 
		 LMM->newPath(horizon,i);                      // needed Libors L_j, j>=i.
		 Real X_iT_i,h,f;
		 X_iT_i=LMM->XL(i,i);                          // Libor X_i(T_i)
		 h=max<Real>(X_iT_i-delta_i*kappa,0.0);        // payoff at time T_{i+1}
		 f=LMM->H_ii(i+1);                             // H_{i+1}(T_{i+1}) accrual T_{i+1}->T_n at time T_{i+1}
		 // move this from time T_{i+1} to time T_n
         return f*h;	
	 }

	 
     /** Mean of the control variate. This is \f$X_i(0)*H_{i+1}(0)\f$,
	  *  see book, 6.9.2.
      */
     Real controlVariateMean() const
     { 
		 return (LMM->X_i0(i))*(LMM->H_i0(i+1));
     } 
 
	 
    /** Control variate is forward transported Libor 
     *  \f$X_i(T_i)H_{i+1}(T_i)\f$, a \f$P_n\f$-martingale. See book, 6.9.2.
     */
    vector<Real> nextControlledForwardPayoff() const
    {
		 LMM->newPath(horizon,i);            // needed Libors L_j, j>=i
		 
         Real X_iT_i=LMM->XL(i,i),                        // Libor X_i(T_i)
		      h=max<Real>(X_iT_i-delta_i*kappa,0.0),      // payoff at time T_{i+1}
		      f=LMM->H_ii(i+1);                           // H_{i+1}(T_{i+1}), accrual T_{i+1}->T_n at time T_{i+1}      
        
         vector<Real> v(2); 
		 v[0]=h*f;                                        // forward transported payoff
		 v[1]=X_iT_i*(LMM->H_it(i+1,i));                  // control variate
		 return v;
     } 
    
    
   
      
// ANALYTIC PRICE

	 
   /** Black caplet price (this is the martingale price in the LMM).
	*  Approximate price if Libor volatility is stochastic as in the
	*  {@link DriftlessLMM}.
    */
   Real analyticForwardPrice() const
   {
       Real Li0=LMM->L_i0(i),                                           // L_i(0)
            cSigma=LMM->capletAggregateVolatility(i),                   // aggregate volatility to expiry    
	        Nplus=FinMath::N(FinMath::d_plus(Li0,kappa,cSigma)),
	        Nminus=FinMath::N(FinMath::d_minus(Li0,kappa,cSigma)),
	        f=LMM->H_i0(i+1);                                           // H_{i+1}(0)=B_{i+1}(0)/B_n(0)
 
       return delta_i*(Li0*Nplus-kappa*Nminus)*f;  
   
   } //end analyticForwardPrice
   

   
// TO STRING

/** Message identifying the object (parameter values, etc)
 */
std::string toString()
{
     ostringstream os;
	 os << "Caplet Cplt([T_"<<i<<",T_"<<i+1<<"],"<<kappa<<")" << endl;

	return os.str();
}


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
    Swaption(int _p, int _q, int _t, Real strike, LiborMarketModel* lmm) :
    LiborDerivative(lmm,_t),                           // Libors needed until time T_s                  
	p(_p), q(_q), t(_t),                                  
	kappa(strike)
    {   
		armControlVariate();
	} 

    /** Sample at the money swaption with p=n/3, q=n, t=p.
	 *  Based on a {DriftlessLMM} with a {@link CS_FactorLoading}.
	 */
	static Swaption* sample(int n, Real delta)
    {
		LiborMarketModel* lmm=DriftlessLMM::sample(n,delta);
		int p=n/3, q=n, t=p;
		Real kappa=lmm->swapRate(p,q);
		return new Swaption(p,q,t,kappa,lmm);
	} 
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff computed from current state of the Libor generator and transported 
	  *  forward from time \f$T_t\f$ to time \f$T_n\f$.
      */
     Real nextForwardPayoff() const 
	 {
		 LMM->newPath(horizon,t);                      // needed Libors L_j, j>=t
		 Real S_pqT,B_pqT,h,f;
		 S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
		 B_pqT=LMM->B_pq(p,q,t);
		 h=B_pqT*max<Real>(S_pqT-kappa,0.0);           // payoff at time T_t
		 f=LMM->H_ii(t);                               // H_t(T_t) accrual T_t->T_n at time T_t
		 // move this from time T_t to time T_n
         return f*h;	
	 }

	 
     /** Mean of the control variate. This is \f$(B_p(0)-B_q(0))/B_n(0)\f$, see book, 6.9.2.
      */
     Real controlVariateMean() const
     { 
		 Real fp=LMM->H_i0(p),        // H_p(0)
		      fq=LMM->H_i0(q);        // H_q(0)
		 return fp-fq;
     } 
 
	 
    /** Control variate is \f$H_p(T_t)-H_q(T_t)\f$, a \f$P_n\f$-martingale.
	 *  See book, 6.9.2.
     */
    vector<Real> nextControlledForwardPayoff() const
    {
		 LMM->newPath(horizon,t);            // needed Libors L_j, j>=t
		 
		 Real S_pqT,B_pqT,h,ft,fp,fq;
		 S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
		 B_pqT=LMM->B_pq(p,q,t);
		 h=B_pqT*max<Real>(S_pqT-kappa,0.0);           // payoff at time T_t
		 ft=LMM->H_ii(t),                              // H_t(T_t) accrual to horizon from time t
		 fp=LMM->H_it(p,t),                            // H_p(T_t)
		 fq=LMM->H_it(q,t);                            // H_q(T_t)
        
         vector<Real> v(2); 
		 v[0]=h*ft;                                     // forward transported payoff
		 v[1]=fp-fq;                                    // control variate
		 return v;
     } 
    
    
   
      
// ANALYTIC PRICE

	 
   /** Black caplet price (this is the martingale price in the LMM).
	*  Approximate price if Libor volatility is stochastic as in the
	*  {@link DriftlessLMM}.
    */
   Real analyticForwardPrice() const
   {
       Real S_pq=LMM->swapRate(p,q),                                    // S_{p,q}(0)
            swpnSigma=LMM->swaptionAggregateVolatility(p,q,t),          // aggregate volatility to T_t    
	        Nplus=FinMath::N(FinMath::d_plus(S_pq,kappa,swpnSigma)),
	        Nminus=FinMath::N(FinMath::d_minus(S_pq,kappa,swpnSigma)),
	        H_pq=LMM->H_pq(p,q);                                        // forward B_{p,q}(0)
  
       return H_pq*(S_pq*Nplus-kappa*Nminus);  
   
   } //end analyticForwardPrice
   


// TO STRING

/** Message identifying the object (parameter values, etc)
 */
std::string toString()
{
     ostringstream os;
	 os << "Swaption Swpn(T_"<<t<<",[T_"<<p<<",T_"<<q<<"],"<<kappa<<")" << endl;

	return os.str();
}



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
	void setLMM(LiborMarketModel* lmm)
	{ 
		LiborDerivative::setLMM(lmm); 
		// we must also reset the libor marke model for the underlying bond
		B->setLMM(lmm);      
	}
	

// CONSTRUCTOR

    /** <p>Call on general bond portfolio \f$D(t)=\sum\nolimits_j=p^{q-1}c_jB_j(t)\f$
     *
     * @param D the underlying bond.
	 * @apram strike the strike price.
	 * @param s call exercisable at time \f$T_s\f$.
     */
    BondCall(Bond* D, Real strike, int s) : 
	LiborDerivative(D->getLMM(),s),      // Libors needed until time T_s
	B(D),
	p(D->get_p()), q(D->get_q()), 
	K(strike), t(s)
    {   
		armControlVariate();
    }
	
	
	/** Sample bond with p=n/3, q=2*n/3, all coupons random in [0.5,1.5].
	 *  Call on this bond with strike rate = cash price of the bond
	 *  at the horizon exercisable at time \f$T_p\f$. 
	 *  Based on a {DriftlessLMM} with a {@link CS_FactorLoading}.
	 */
	static BondCall* sample(int n, Real delta)
    {
		LiborMarketModel* lmm=DriftlessLMM::sample(n,delta);
		int p=n/3, q=2*n/3, t=p;
		vector<Real> c(q-p,p);
		for(int j=p;j<q;j++) c[j]=0.5+Random::U01();
		Bond* B=new Bond(p,q,c,lmm);
		Real K=B->cashPrice();
		
		return new BondCall(B,K,t);
	}
	 
	/** Sample call on zero coupon bond with i=n/2, strike price = cash price 
	 *  of the bond at the horizon exercisable at time \f$T_{i-1}\f$. 
	 *  Based on a {DriftlessLMM} with a {@link CS_FactorLoading}.
	 *  This is a worst case for the assumptions of the analytic price formulas.
	 */
	static BondCall* sampleCallOnZeroCouponBond(int n, Real delta)
    {
		LiborMarketModel* lmm=DriftlessLMM::sample(n,delta);
		int i=n/2, t=i-1;
		Bond* B=new Bond(i,lmm);
		Real K=B->cashPrice();
		
		return new BondCall(B,K,t);
	}

     
	 
    
// FORWARD TRANSPORTED PAYOFF
	
	
	    
     /** Payoff \f$h\f$ at time \f$T_t\f$ accrued forward to time \f$T_n\f$:
	  *  \f[h/B_n(T_t)=h*H_t(T_t)=\left(\sum_{j=p}^{q-1}c_jH_j(T_t)-KH_t(T_t)\right)^+\f]
	  * computed from a new Libor path.
      */
     Real nextForwardPayoff() const 
	 {  
		 LiborMarketModel* lmm=B->getLMM();
		 lmm->newPath(horizon,t);             // needed Libors L_j, j>=t
		 Real F=B->forwardPrice(t),
		      Ht=lmm->H_ii(t);                // H_t(T_t)

		 return max<Real>(F-K*Ht,0.0);
	 }
//<----------- INVESTIGATE ----------->
// if we use LMM instead of lmm=B->getLMM(), then the paths of the Libor model underlying B
// don't move even with our initilization of LiborDerivative::LMM as B->getLMM(), why is that???		

	 
     /** Mean of the control variate, the forward price at time zero.
      */
     Real controlVariateMean() const
     { 
		 return B->forwardPrice();
     } 
 
	 
    /** Control variate is the forward bond price.
     */
    vector<Real> nextControlledForwardPayoff() const
    {
		 LiborMarketModel* lmm=B->getLMM();
		 lmm->newPath(horizon,t);              // needed Libors L_j, j>=t
		 Real F=B->forwardPrice(t),
		      Ht=lmm->H_ii(t);                 // H_t(T_t)
        
         vector<Real> v(2); 
		 v[0]=max<Real>(F-K*Ht,0.0);	       // forward transported payoff
		 v[1]=B->forwardPrice(t);              // control variate
		 return v;
     } 
		 
   
      
// ANALYTIC PRICE

	 
   /** Black-Scholes price, assumes bond price at expiration is lognormal.
	*  See book, 6.9.4.
    */
   Real analyticForwardPrice() const
   {
       LiborMarketModel* lmm=B->getLMM();
	   Real F=B->forwardPrice(),
	        H_t0=lmm->H_i0(t),                                          // H_t(0)
	        Q=F/H_t0,
            bondSigma=LMM->bondAggregateVolatility(B,t),                // aggregate volatility to T_t    
	        Nplus=FinMath::N(FinMath::d_plus(Q,K,bondSigma)),
	        Nminus=FinMath::N(FinMath::d_minus(Q,K,bondSigma));
	                                             
        return F*Nplus-K*H_t0*Nminus;  

   } //end analyticForwardPrice
   
   
  
   

// TO STRING

/** Message identifying the object (parameter values, etc)
 */
std::string toString()
{
     ostringstream os;
	 os << "Call on bond B with strike K =: " << K << endl
		<< "Bond B:" << B->toString() << endl;

	return os.str();
}

   

}; // end Bond





MTGL_END_NAMESPACE(Martingale)

#endif
