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
#include "Pricing.h"                   // needed in template function
#include "Node.h"                      // problem with template forward declarations
#include <cstdlib>                     // exit(int)




MTGL_BEGIN_NAMESPACE(Martingale)


// dependencies
class LiborMarketModel;
class Bond;
// class RealVector;
class StandardBrownianNode;
class LmmLattice; 




/*! \file Option.h
 *  <p>The class {@link Option} is the general interface to all European or 
 *  Bermudan options. In fact Bermudan options allow exercise at an arbitrary 
 *  sequence of times and European respectively American options are treated as 
 *  a special case where exercise is possible only at expiry respectively at 
 *  each time step. See {@link Option}, {@link LiborDerivative}, 
 *  {@link Caplet}, {@link Swaption} and {@link BondCall}.
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
	 bool hasAnalytic_,       // analytic pricing formula implemented
	      hasLattice_,        // evaluation in lattice implemented
	      hasMonteCarlo_,     // Monte Carlo path pricing implemented
	      hasControlVariate_; // control variate is implemented
       

public:								 
	
    /** Number of time steps to expiration.*/
	Real getExpiration(){ return T_; }

    /** @param T time to expiry in years.
	 *  @param an analytic pricing formula is implemented.
	 *  @param lt lattice pricing is implemented.
	 *  @param mc Monte Carlo pricing is implemented.
	 *  @param cv control variate is implemented.
	 */
    Option(Real T, bool an, bool lt, bool mc, bool cv) : 
	T_(T), 
	hasAnalytic_(an), hasLattice_(lt), hasMonteCarlo_(mc), hasControlVariate_(cv)
	{   }
	
	
// PAYOFF

/** Wether the option can be exercised at continuous time t < expiration.
 *  Returns false defaulting to European exercise. Override as appropriate.*/
virtual bool isExercisable(Real t){ return false; }	
		    
/** The next forward payoff random sample.*/
virtual Real nextForwardPayoff(){ return -1.0; }
	 
/** Returns default 0.0. Override. */
virtual Real controlVariateMean(){ return 0.0; }
	 
/** Vector Z with Z[0] the next forward payoff random sample 
 *  and Z[1] the corresponding control variate. Default returns 
 *  (-1.0,0.0). Return by value since it is so small.
 */
virtual const RealArray1D nextControlledForwardPayoff();

/** Forward compounded payoff of the option if exercised at node node
 *  living at time step t in lattice theLattice. Returns -1.0. 
 *  To be overriddden by appropriate specializations.
 */
template<typename LatticeType, int t>
Real forwardPayoff
(typename LatticeType::NodeType* node, LatticeType* theLattice, int t){ return -1.0; }
	
	 
/** Error message that nothing is implemented in this generality. 
 *  Returns -1.
 */
virtual Real analyticForwardPrice() const;
	
/** Forward price computed from N sample payoffs of the option.*/
Real monteCarloForwardPrice(int N);

/** Forward price computed from from N sample payoffs of the option
 *  using the control variate implemented with the option.
 */
Real controlledMonteCarloForwardPrice(int nPaths);

/** The correlation of the forward payoff with its control variate
 *  computed from N sample payoff - controlVariate pairs of the option.
 */
Real correlationWithControlVariate(int N);

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

/** Class which factors out some common features of all Libor derivatives.
 *  PathGenerators for Libor derivatives see to it that only those Libors
 *  are computed which are really needed for the option payoff. 
 *
 * <p>To compound the payoff forward from time the \f$T_t\f$ of exercise to 
 *  the horizon we need all Libors \f$X_j,\ j\geq t\f$. The standing
 *  assumption is that the underlying LMM makes one time step per Libor
 *  compounding period. This class does not take ownership of the underlying 
 *  LMM (delete separately).
 */
class LiborDerivative : public Option {
	
protected:
	
	LiborMarketModel* LMM;     // the underlying LMM
	int n;                     // number of Libors including L_0
	int t;                     // Libors are needed until time T_t for option payoff

public:
	
/** @param lmm underlying {@link LiborMarketModel}.
 *  @param s option expires at \f$T_s\f$. 
 *  @param an analytic pricing formula is implemented.
 *  @param lt lattice pricing is implemented.
 *  @param mc Monte Carlo pricing is implemented.
 *  @param cv control variate is implemented.
 */
// T_t is handed to Option as "time to expiration".
LiborDerivative
(LiborMarketModel* lmm, int s, bool an, bool lt, bool mc, bool cv);

virtual ~LiborDerivative(){ }
	
/** Number t of Libor compounding periods to expiry (at T_t).*/
// Usually this is t but not always, see Caplet where it is t-1.
virtual int getPeriodsToExpiry(){ return t; }
	
/** Dimension n of the underlying Libor market model. */
int  getDimension() const { return n; }
	
/** The effective dimension of computing one payoff. This is the number
 *  of independent uniform deviates generated to compute one sample of
 *  the necessary Libors.
 */
int effectiveDimension() const;


/** The default lattice used for pricing: a two factor
 *  LmmLattice based on the underlying LMM. Constructed on the heap
 *  and must be deallocated by the user. This lattice is built with
 *  steps=6 time steps per Libor compounding period unless the constant
 *  <code>LATTICE_MAX_STEPS</code> defined in TypedefsMacros.h forces
 *  reduction in the number steps. This is adequate up to semiannual 
 *  compounding. The maximum number of possible time steps depends on
 *  the available RAM. We assume 1GB. Otherwise adjust
 *  <code>LATTICE_MAX_STEPS</code> by trial and error.
 */
LmmLattice* getDefaultLattice();


/** Payoff at a node living at time step s in theLattice. 
 *  Empty default implementation: aborts with error message.
 */
virtual Real forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s);


/** Price computed in theLattice. */
Real latticeForwardPrice(LmmLattice* theLattice);
	
	
/** Forward price computed in the default lattice.
 *  See {@link getDefaultLattice()}.
 */
Real latticeForwardPrice(){ return latticeForwardPrice(getDefaultLattice()); }

	
/** Message announcing a generic Libor derivative. */
std::ostream& printSelf(std::ostream& os) const;
	

// TESTING PRICING ROUTINES

/** Tests the analytic price against Lattice, Monte Carlo
 *  and Quasi Monte Carlo prices with and without control variates.
 *  Reports error relative to the analytic price if this price is defined. 
 *  Note however that the analytic price itself may be an approximation.
 *
 * @param nPath number of payoff samples (usually computed from Libor paths).
 */
virtual void testPrice(int nPath);

		
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
int getPeriodsToExpiry(){ return i; }


// CONSTRUCTOR

/** Caplet on \f$[T_j,T_{j+1}]\f$ based on a Libor process 
 *  with strike rate <code>strike</code>.
 *
 * @param j caps Libor on \f$[T_j,T_{j+1}]\f$.
 * @param strike strike rate.
 * @param lmm underlying Libor market model.
 */
Caplet(int j, Real strike, LiborMarketModel* lmm);
	

/** Sample at the money caplet with i=n/3, ie. one third of the way out to the horizon.
 *  Based on a {@link DriftlessLMM}.
 *
 * @param n dimension (number of Libor accrual intervals).
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,LFDL,PC,FPC.
 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
 * @param corrType type of correlations, Correlations::CS, JR.
 */
static Caplet* sample
(int n, 
 int lmmType = 0,    // LiborMarketModel::DL, we use literals to avoid includes 
 int volType = 2,    // VolSurface::CONST
 int corrType =1     // Correlations::CS
);
		

// PAYOFF

/** Next random sample of the payoff compounded forward 
 *  from time \f$T_{i+1}\f$ to time \f$T_n\f$.
 */
Real nextForwardPayoff();

/** Mean of the control variate. 
 *  This is \f$(B_p(0)-B_q(0))/B_n(0)\f$, see book, 6.9.2.*/
Real controlVariateMean();

/** Vector Z with Z[0] the next forward payoff random sample and Z[1] 
 *  the corresponding control variate. Return by value since it is so small.*/
const RealArray1D nextControlledForwardPayoff();

/** Payoff at node compounded forward to time \f$T_n\f$.
 *  Specialization overrides base class template.
 *  Slightly inaccurate (see {@link LiborFunctional#capletForwardPayoff})
 *  but we can do no better in a lattice.
 *
 *  @param s lattice time step at which the node lives.
 */
Real forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s);



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
 *  exercisable at time \f$T_t\leq T_p\f$ pays off 
 *  \f$h=B_{p,q}(T_t)*(S_{p,q}(T_t)-\kappa)^+\f$ 
 *  at time \f$T_t\f$. Here \f$B_{p,q}\f$ and \f$S_{p,q}\f$ are the annuity and swap 
 *  rate along \f$[T_p,T_q]\f$ as usual. It is based on a {@link LiborMarketModel}.
 */
class Swaption : public LiborDerivative {
   
	int  p,q,                // swap period [T_p,T_q]	
	     t;                  // swaption exercises at T_t
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
Swaption(int p, int q, int t, Real strike, LiborMarketModel* lmm);
	

/** Sample at the money swaption with t=p.
 *	 
 * @param p,q period of swap \f$[T_p,T_q]\f$.
 * @param int lmmType type of Libor market model: {@link LiborMarketModel}::DL,LFDL,PC,FPC
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

/** Next random sample of the payoff payoff compounded forward 
 *  from time \f$T_t\f$ to time \f$T_n\f$.
 */
Real nextForwardPayoff();

/** Mean of the control variate. 
 *  This is \f$H_p(0)-H_q(0)\f$, see book, 6.9.2.*/
Real controlVariateMean();

/** Vector Z with Z[0] the next forward payoff random sample and Z[1] 
 *  the corresponding control variate. Return by value since it is so small.*/
const RealArray1D nextControlledForwardPayoff();


/** Payoff at node compounded forward to time \f$T_n\f$.
 *  Specialization overrides base class template.
 *
 *  @param s lattice time step at which the node lives.
 */
Real forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s);


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


/** Call on a general bond \f$B(s)=\sum\nolimits_{j=p}^{q-1}c_jB_j(s)\f$
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
 * @param strike the strike price.
 * @param s call exercisable at time \f$T_s\f$.
 */
BondCall(Bond* D, Real strike, int s);
		
/** Sample call on bond with random coupons in [0.5,1.5] received at 
 *  \f$T_j,\ p\leq j<q\f$. Strike price = cash price of the bond,
 *  expiration at \f$T_p\f$. 
 *	 
 * @param p coupons begin at \f$T_p\f$.
 * @param q coupons end at \f$T_q\f$.
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,LFDL,PC,FPC.
 * @param volType type of volatility surface, {@link VolSurface}::CONST, JR, M.
 * @param corrType type of correlations, {@link Correlations}::CS, JR.
 */
static BondCall* sample
(int p, int q, 
 int lmmType = 0,    // LiborMarketModel::DL, we use literals to avoid includes 
 int volType = 2,    // VolSurface::CONST
 int corrType =1     // Correlations::CS
);
	
	 
/** Sample call on zero coupon bond maturing at \f$T_p\f$ 
 *  exercisable at time \f$T_{p-1}\f$, strike price = cash price  
 *  This is a worst case for the assumptions of the analytic price formulas.
 *  Dimension of LMM is chosen to be p+3.
 *	 
 * @param p bond matures at \f$T_p\f$
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,LFDL,PC,FPC.
 * @param volType type of volatility surface, {@link VolSurface}::CONST, JR, M.
 * @param corrType type of correlations, {@link Correlations}::CS, JR.
 */
static BondCall* sampleCallOnZeroCouponBond
(int p, 
 int lmmType = 0,    // LiborMarketModel::DL, we use literals to avoid includes 
 int volType = 2,    // VolSurface::CONST
 int corrType =1     // Correlations::CS
);

        
// FORWARD TRANSPORTED PAYOFF AND CONTROL VARIATE

/** Next random sample of the payoff payoff compounded forward 
 *  from time \f$T_t\f$ to time \f$T_n\f$.
 */
Real nextForwardPayoff();

/** Mean of the control variate. 
 *  This is the forward price of the bond at time zero.*/
Real controlVariateMean();

/** Vector Z with Z[0] the next forward payoff random sample and Z[1] 
 *  the corresponding control variate. Return by value since it is so small.*/
const RealArray1D nextControlledForwardPayoff();

/** Payoff at node compounded forward to time \f$T_n\f$.
 *  Specialization overrides base class template.
 *  
 *  @param s lattice time step at which the node lives.
 */
Real forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s);
   
/** See book, 6.9.4.*/
Real analyticForwardPrice() const;
   
/** Message identifying the derivative. */
std::ostream& printSelf(std::ostream& os) const;

   
}; // end BondCall





MTGL_END_NAMESPACE(Martingale)

#endif
