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


#ifndef martingale_liborfunctional_h    
#define martingale_liborfunctional_h

#include "TypedefsMacros.h"
#include "Array.h"


MTGL_BEGIN_NAMESPACE(Martingale)

// dependencies
class BondCall;


/*! \file LiborFunctional.h
 *  Free standing functions to computing various functionals F of the
 *  Libor process as deterministic functions \f$F=F(H)\f$ of a vector 
 *  \f[H=(H_p(t),\dots,H_n(t))\f] 
 *  of accrual factors at time t.  
 */

	
/** <p>Free standing functions to computing various functionals F of the
 *  Libor process as deterministic functions \f$F=F(H)\f$ of a vector 
 *  \f[H=(H_p(t),\dots,H_n(t))\f] 
 *  of accrual factors at time t. It is assumed that the vector
 *  H has natural indexation \f$H[j]=H_j(t),\ j=p,\dots,n\f$ and of course
 *  that it contains all the accrual factors which are needed for the computation
 *  of F. 
 *
 *  <p><b>Assumption:</b> all Libor compounding periods have the same length.
 *  This class has to be rewritten if we work with LMM with irregular tenor 
 *  structure.
 *
 *  <p>Note that the time t does not explicitly enter in any of the
 *  computations. This information is in the vector H.
 */
namespace LiborFunctional {
	
	/** Libor X_j.
	 * @param H vector of accrual factors.
	 */
	inline Real X(int j, const RealArray1D& H){ return (H[j]-H[j+1])/H[j+1]; }
	
	/** Forward annuity (forward price of a basis point) on [T_p,T_q].
	 * @param H vector of accrual factors.
	 * @param delta length of all Libor accrual periods.
	 */
	Real H_pq(int p, int q, const RealArray1D& H, Real delta);	
		
	/** Swaprate for swap on [T_p,T_q].
	 * @param H vector of accrual factors.
	 * @param delta length of all Libor accrual periods.
	 */
	Real swapRate(int p, int q, const RealArray1D& H, Real delta);
		
		
	/** Payoff of a forward swaption with strike rate kappa exercising into a 
     *  swap on the interval [T_p,T_q]. Payoff is accrued forward to the horizon T_n. 
	 *
	 * @param H vector of accrual factors.
	 * @param delta length of all Libor accrual periods.
     */
    Real forwardSwaptionPayoff
	(int p, int q, Real kappa, const RealArray1D& H, Real delta);
	
	
	/** Payoff of a caplet with strike rate kappa on the interval [T_i,T_{i+1}]. 
	 *  Payoff is accrued forward to the horizon T_n. It is assumed that the
	 *  vector H is computed at time T_i.
	 *
	 * @param H vector of accrual factors.
	 * @param delta length of all Libor accrual periods.
     */
    Real forwardCapletPayoff
	(int i, Real kappa,const RealArray1D& H, Real delta);
	
	
	/** Payoff of the {@link BondCall} accrued forward to the horizon T_n. 
	 *
	 * @param bc the BondCall.
	 * @param H vector of accrual factors.
     */
    Real forwardBondCallPayoff(BondCall* bc, const RealArray1D& H);

};



MTGL_END_NAMESPACE(Martingale)

#endif


