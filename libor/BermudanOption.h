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
spyqqqdia@yahoo.com

*/


 
 
#ifndef martingale_bermudanoption_h    
#define martingale_bermudanoption_h

#include "Option.h"
#include "Trigger.h"



MTGL_BEGIN_NAMESPACE(Martingale)

// dependencies
class LiborMarketModel;
class LmmNode;



/*! \file BermudanOption.h
 *  <p>Some Bermudan options. Exercise is allowed at an arbitrary sequence of
 *  times. American options are treated as the special case where exercise is
 *  possible at each time step.
 */
 

/*******************************************************************************
 *
 *                     BERMUDAN SWAPTION
 *
 ******************************************************************************/
 

/** Bermudan payer swaption on \f$[T_p,T_q]\f$. Exercise is possible at each Libor
 *  reset date \f$T_j,\ p\leq j<q\f$. No analytic pricing formula. We use the
 *  exercise trigger described in P. Jaeckel, "Monte Carlo Methods in Finance",
 *  12.7. This trigger is much faster than general convex exercise (book, 4.6.3)
 *  and produces slightly better results. No control variate is implemented.
 */ 
class BermudanSwaption : public LiborDerivative {
   
	int  p,q;                // initial swap period [T_p,T_q]
	int nPath;               // number of training paths for the exercise trigger
    Real kappa;              // strike rate
	
	PjTrigger trigger;       // the exercise trigger
	friend class PjTrigger;
		     
public:

/** The srike rate. */
Real getStrike(){ return kappa; }

/** Wether or not the option can be exercised at continuous time t.*/
bool isExercisable(Real t);

// CONSTRUCTOR

/** 
 * @param p,q period of swap \f$[T_p,T_q]\f$.
 * @param paths number of training paths for the exercise trigger.
 * @param strike strike rate.
 * @param lmm underlying Libor market model.
 * @param verbose messages during trigger optimization
 */
BermudanSwaption
(int p_, int q_, int paths, Real strike, LiborMarketModel* lmm, bool verbose=false);
	

/** Sample at the money Bermudan swaption.
 *	 
 * @param p,q swap interval \f$[T_p,T_q]\f$.
 * @param paths number of training paths for the exercise trigger.
 * @param verbose messages during trigger optimization.
 * @param int lmmType type of Libor market model: {@link LiborMarketModel::DL,PC,FPC}
 * @param volType type of volatility surface, VolSurface::CONST, JR, M.
 * @param corrType type of correlations, Correlations::CS, JR.
 */
static BermudanSwaption* sample
(int p, int q, int paths, 
 bool verbose=false,
 int lmmType = 0,    // LiborMarketModel::DL, we use literals to avoid includes 
 int volType = 2,    // VolSurface::CONST
 int corrType =1     // Correlations::CS
);
	
	 
    
// FORWARD TRANSPORTED PAYOFF AND CONTROL VARIATE

		    
/** Next random sample of the swaption payoff (resulting from the exercise
 *  strategy implemented by the trigger) compounded forward to time \f$T_n\f$.*/
Real nextForwardPayoff();

/** Swaption payoff if exercised at time \f$T_s\f$ in the current Libor path
 *  compounded forward to time \f$T_n\f$.
 */
Real currentForwardPayoff(int s);

/** Payoff at LmmNode compounded forward to time \f$T_n\f$.*/
Real forwardPayoff(LmmNode* node);
  
/** Message identifying the derivative. */
std::ostream& printSelf(std::ostream& os) const;

}; // end swaption



MTGL_END_NAMESPACE(Martingale)

#endif
