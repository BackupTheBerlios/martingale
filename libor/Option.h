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
 
 
#ifndef martingale_option_h    
#define martingale_option_h

#include "TypedefsMacros.h"
#include "PathGenerator.h"           // base class
#include "Pricing.h"                 // needed in template function
#include <cstdlib>                   // exit(int)




MTGL_BEGIN_NAMESPACE(Martingale)


// dependencies
class LiborMarketModel;
class Bond;
// class RealVector;
class LmmLattice;
class LmmNode;


/*! \file Option.h
 *  <p>The class {@link Option} is the general interface to all European or 
 *  Bermudan options. In fact Bermudan options allow exercise at an arbitrary 
 *  sequence of times and European Options are treated as a special case where 
 *  exercise is possible only at expiry.
 */
 

/*******************************************************************************
 *
 *                     Option
 *
 ******************************************************************************/

/** General Bermudan option. Exercise can be allowed at an arbitrary sequence of
 *  times. European and American options are special cases (exercise only at 
 *  expiry respectively at each time step). Terminology uses forward prices at the
 *  horizon (as is usual for Libor derivatives) not discounted prices (as is more 
 *  common for most other derivatives). This works in all cases. The horizon is
 *  usually the date of option expiry. In the case of Libor derivatives the 
 *  horizon is the end of the last Libor accrual period regardless of the
 *  date of option expiry.
 */ 
class Option {
	
protected:
    
	 Real T_;                 // time to expiry in years
	 bool hasLattice_,        // evaluation in lattice implemented
	      hasMonteCarlo_,     // Monte Carlo path pricing implemented
	      hasControlVariate_; // control variate is implemented
       

public:								 
	
    /** Number of time steps to expiration.*/
	Real getExpiration(){ return T_; }

    /** @param T time to expiry in years.
	 *  @param lt lattice pricing is implemented.
	 *  @param mc Monte Carlo pricing is implemented.
	 *  @param cv control variate is implemented.
	 */
    Option(Real T, bool lt, bool mc, bool cv) : 
	T_(T), hasLattice_(lt), hasMonteCarlo_(mc), hasControlVariate_(cv)
	{   }
	
	
// PAYOFF

/** Wether the option can be exercised at continuous time t < expiration.
 *  Returns false defaulting to European exercise. Override as appropriate.*/
virtual bool isExercisable(Real t){ return false; }	
		    
/** Returns -1.0. Override meaningfully.*/
virtual Real forwardPayoffAlongCurrentPath(){ return -1.0; }
	 
/** Returns default 0.0. Override. */
virtual Real controlVariateMean(){ return 0.0; }
	 
/** Returns 0.0. Override.*/
virtual Real controlVariateAlongCurrentPath(){ return 0.0; }

/** Returns -1.0. Override meaningfully.*/
virtual Real forwardPayoff(LmmNode* node){ return -1.0; }
	
	 
/** The PathGenerator associated with the option.
 *  Error message that nothing is implemented in this generality.
 *  Aborts.
 */
virtual PathGenerator* getPathGenerator();

/** Error message that nothing is implemented in this generality. 
 *  Returns -1.
 */
virtual Real analyticForwardPrice() const;
	
/** Forward price computed from nPaths paths of the PathGenerator
 *  associated with the option.
 */
Real monteCarloForwardPrice(int nPaths);

/** Forward price computed from nPaths paths of the PathGenerator
 *  associated with the option using the control variate implemented 
 *  with the option.
 */
Real controlledMonteCarloForwardPrice(int nPaths);

/** Forward price computed in theLattice.*/
// must be a template since Lattice is a class template not a class.
template<class LatticeType>
Real latticeForwardPrice(LatticeType* theLattice)
{
	return Pricing::latticeForwardPrice(theLattice,this);
}

/** Forward price computed in the default lattice.
 *  Error message that nothing is implemented in this generality. 
 *  Returns -1.
 */
virtual Real latticeForwardPrice();

/** Message descriptive of derivative. */
virtual std::ostream& printSelf(std::ostream& os) const;
      
     
}; // end Derivative






/*******************************************************************************
 *
 *                          LIBOR DERIVATIVE
 *
 ******************************************************************************/

/** PathGenerator for Libor derivatives. Maintains reference to underlying 
 *  LiborMarketModel and forwards request for a newPath() to the LMM while
 *  controlling which Libors are computed until which time.
 */
struct LiborPathGenerator : public PathGenerator {
	
	LiborMarketModel* LMM;
	/** Libors are computed until time T_t.*/   int t;                     
	/** Libors evolved are L_j, j>=i.*/         int i;                     
	
public:
	
	/** @param lmm the underlying LiborMarketModel.
	 *  @param s time until which Libors are computed.
	 *  @param k Libors evolved are L_j, j>=k.
	 */
	LiborPathGenerator(LiborMarketModel* lmm, int s, int k) :
	LMM(lmm), t(s), i(k)
    {   }
	
	void newPath();
	
};
	
	

/** Class which factors out some common features of all Libor derivatives.
 *  PathGenerators for Libor derivatives see to it that only those Libors
 *  are computed which are really needed for the option payoff. The standing
 *  assumption is that the underlying LMM makes one time step per Libor
 *  compounding period.
 */
class LiborDerivative : public Option {
	
protected:
	
	LiborMarketModel* LMM;     // the underlying LMM
	LiborPathGenerator* LPG;   // path generator, maintains reference to underlying LMM
	int n;                     // number of Libors including L_0
	int t;                     // Libors needed until time T_t for option payoff


public:
	
/** @param lmm underlying Libor market model.
 *  @param s Libors are needed until time \f$T_s\f$. 
 *  @param i Libors needed are L_j, j>=i.
 *  @param lt lattice pricing is implemented.
 *  @param mc Monte Carlo pricing is implemented.
 *  @param cv control variate is implemented.
 */
// T_t is handed to Option as "time to expiration".
LiborDerivative
(LiborMarketModel* lmm, int s, int i, bool lt, bool mc, bool cv);
	
/** Number t of Libor compounding periods to expiry (at T_t).*/
// Usually this is t but not always, see Caplet where it is t-1.
virtual int getPeriodsToExpiry(){ return t; }
	
/** Dimension n of the underlying Libor market model. */
int  getDimension() const { return n; }
	
/** The effective dimension of computing one payoff. This is the number
 *  of independent uniform deviates generated to compute the necessary 
 *  Libors. Default implementation returns zero.
 */
int effectiveDimension() const;
	
/** The PathGenerator. */
PathGenerator* getPathGenerator(){ return LPG; }
	
/** The default lattice used for pricing: a two factor
 *  LmmLattice based on the underlying LMM. Constructed on the heap
 *  and must be deallocated by the user.
 */
LmmLattice* getDefaultLattice();

/** Forward price computed in the default lattice.
 *  See {@link getDefaultLattice()}.
 */
Real latticeForwardPrice();
	
/** Message announcing a generic Libor derivative. */
std::ostream& printSelf(std::ostream& os) const;
	

// TESTING PRICING ROUTINES

/** Tests the analytic price against Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates from
 *  a sample of 20000 paths of the underlying Libor process.
 */
virtual void testPrice();

		
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
    Real kappa,              // strike rate
	     delta_i;            // T_{i+1}-T_i
	

	     
public: 
	
/** Number t of Libor compounding periods to expiry (at T_t).*/
int getPeriodsToExpiry(){ return t-1; }
	

// CONSTRUCTOR

/** Caplet on \f$[T_j,T_{j+1}]\f$ based on a Libor process 
 *  with strike rate <code>strike</code>.
 *
 * @param j caps Libor on \f$[T_j,T_{j+1}]\f$.
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
		

// PAYOFF

/** Forward payoff along the current Libor path.*/
Real forwardPayoffAlongCurrentPath();

/** Mean of the control variate. This is \f$(B_p(0)-B_q(0))/B_n(0)\f$, see book, 6.9.2.*/
Real controlVariateMean();

/** Control variate along the current Libor path.*/
Real controlVariateAlongCurrentPath();


// Payoffs at nodes not implemented. Delay problem: payoff determined at time T_i
// but received at time T_{i+1}. So forward transporting involves (one step)
// path history. Not amenable to lattice treatment.
      
// ANALYTIC PRICE

	 
/** Black caplet price. See book 6.6, 6.8.2.
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
	                        
    Real kappa;              // strike rate
	

	     
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
    Swaption(int p_, int q_, int t_, Real strike, LiborMarketModel* lmm);
	

    /** Sample at the money swaption with t=p.
	 *	 
 	 * @param n dimension (number of Libor accrual intervals).
	 * @param int lmmType type of Libor market model: {@link LiborMarketModel::DL,PC,FPC}
	 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
	 * @param corrType type of correlations, Correlations::CS, JR.
	 */
	static Swaption* sample
	(int p, int q,
	 int lmmType = 0,    // LiborMarketModel::DL, we use literals to avoid includes 
	 int volType = 2,    // VolSurface::CONST
	 int corrType =1     // Correlations::CS
	);
	
	 
    
// FORWARD TRANSPORTED PAYOFF AND CONTROL VARIATE
		    
/** Swaption payoff compounded forward from time \f$T_t\f$ to time \f$T_n\f$.*/
Real forwardPayoffAlongCurrentPath();
	 
/** This is \f$H_p(0)-H_q(0)\f$, see book, 6.9.2.*/
Real controlVariateMean();
	 
/** Control variate is \f$H_p(T_t)-H_q(T_t)\f$, a \f$P_n\f$-martingale. See book, 6.9.2.*/
Real controlVariateAlongCurrentPath();

/** Payoff at LmmNode compounded forward to time \f$T_n\f$.*/
Real forwardPayoff(LmmNode* node);

// PRICING
	 
/** See book 6.7, 6.8.3.*/
Real analyticForwardPrice() const;
   
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
	
/** The underlying bond.*/
Bond* getBond(){ return B; }

/** The strike price. */
Real getStrike(){ return K; }

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

        
// FORWARD TRANSPORTED PAYOFF AND CONTROL VARIATE
		    
/** Swaption payoff compounded forward from time \f$T_t\f$ to time \f$T_n\f$.*/
Real forwardPayoffAlongCurrentPath();
	 
/** This is \f$H_p(0)-H_q(0)\f$, see book, 6.9.2.*/
Real controlVariateMean();
	 
/** Control variate is \f$H_p(T_t)-H_q(T_t)\f$, a \f$P_n\f$-martingale. See book, 6.9.2.*/
Real controlVariateAlongCurrentPath();

/** Payoff at LmmNode compounded forward to time \f$T_n\f$.*/
Real forwardPayoff(LmmNode* node);
   
/** See book, 6.9.4.*/
Real analyticForwardPrice() const;
   
/** Message identifying the derivative. */
std::ostream& printSelf(std::ostream& os) const;

   
}; // end BondCall





MTGL_END_NAMESPACE(Martingale)

#endif
