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

	
/** Free standing functions to computing various functionals F of the
 *  Libor process as deterministic functions 
 *  \f[F=F(H)\hbox{ of a vector } H=(H_p(t),\dots,H_n(t))\f] 
 *  of accrual factors at time \f$t\f$. It is assumed that the vector
 *  H has natural indexation \f$H[j]=H_j(t),\ j=p,\dots,n\f$ and of course
 *  that it contains all the accrual factors whih are needed.
 */
namespace LiborFunctional {
	
	/** Libor X_j.
	 * @param H vector of accrual factors.
	 */
	inline Real X(int j, const RealArray1D& H){ return (H[j]-H[j+1])/H[j+1]; }
	
	/** Forward annuity (forward price of a basis point) on [T_p,T_q].
	 * @param H vector of accrual factors.
	 * @param delta vector of accrual periods.
	 */
	Real Hpq(int p, int q, const RealArray1D& H, const RealArray1D& delta);	
		
	/** Swaprate for swap on [T_p,T_q].
	 * @param H vector of accrual factors.
	 * @param delta vector of accrual periods.
	 */
	Real S_pq(int p, int q, const RealArray1D& H, const RealArray1D& delta);

};
