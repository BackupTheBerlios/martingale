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

#ifndef martingale_bond_h    
#define martingale_bond_h

#include "TypedefsMacros.h"
#include "Array.h"                // problem with typedefs in forward declarations

MTGL_BEGIN_NAMESPACE(Martingale)


// dependencies
class LiborMarketModel;


/*! \file Bond.h
 * Bonds in an LMM with an arbitrary sequence of coupon payments.
 * There is no principal, the principal is part of the last coupon.
 */


/*******************************************************************************
 *
 *                          BONDS
 *
 ******************************************************************************/



/** General <a name="bond">bond</a>
 *  \f[B(t)=\sum\nolimits_j=p^{q-1}c_jB_j(t),\f]
 *  that is, a linear combination of zero coupon bonds based on a 
 *  Libor market model. It is assumed that \f$c_j\geq0\f$.
 *  {@link BondCall}.
 *
 */
class Bond {
   
	LiborMarketModel* LMM;  // the object generating the Libors

	int n,                  // number of Libors including L_0
	    p,q;                // coupon periods [T_p,T_q]

	RealArray1D c;          // c_j units of B_j(t), j=p,...,q-1
	RealArray1D b;          // b_j=c_p+...+c_j, j=p,...,q-1

	     
public: 
	
/** The underlying Libor market model */
LiborMarketModel* getLMM() const { return LMM; }
	
/** Resets the underlying Libor market model.*/
void setLMM(LiborMarketModel* lmm){ LMM=lmm; }
	
/** See <a href="bond">bond definition</a>. */
int get_p() const { return p; }
	
/** See <a href="bond">bond definition</a>. */
int get_q() const { return q; }
	
/** See <a href="bond">bond definition</a>. */
const RealArray1D& get_c() const { return c; }
	
/** \f$b_j=c_p+\dots+c_j, p\leq j<q\f$.*/
const RealArray1D& get_b() const { return b; }

// CONSTRUCTOR

/** <p>General bond portfolio, ie. linear combination \f$B\f$ of zero coupon bonds
 *  \f[B(t)=\sum\nolimits_j=p^{q-1}d_jB_j(t)\f]
 *
 * @param k,m = p,q.
 * @param d coupons.
 * @param lmm underlying Libor market model.
 */
Bond(int k, int m, const RealArray1D& d, LiborMarketModel* lmm);
	
	
/** <p>Zero coupon bonds maturing at time \f$T_i\f$.
 *
 * @param i bond matures at time \f$T_i\f$.
 * @param lmm underlying Libor market model.
 */
Bond(int i, LiborMarketModel* lmm);


// PPRICE
	
/** Forward price \f$F(T_t)=\sum\nolimits_{j=p}^{q-1}c_jH_j(T_t)\f$ at time
 *  \f$T_t\leq T_p\f$ from current Libor path.
 */
Real forwardPrice(int t) const;

/** Forward price \f$F(0)=\sum\nolimits_{j=p}^{q-1}c_jH_j(0)\f$ at time \f$t=0\f$.*/
Real forwardPrice() const;
	 
/** Cash price \f$B(T_t)=\sum\nolimits_{j=p}^{q-1}c_jB_j(T_t)\f$ at time
 *  \f$T_t\leq T_p\f$ from current Libor path.
 */
Real cashPrice(int t) const;

/** Cash price \f$B(0)=\sum\nolimits_{j=p}^{q-1}c_jB_j(0)\f$ at time \f$t=0\f$.*/
Real cashPrice() const;
	 
/** Message identifying the object.*/
std::ostream& printSelf(std::ostream& os) const;
	
}; // end Bond



MTGL_END_NAMESPACE(Martingale)	

#endif


