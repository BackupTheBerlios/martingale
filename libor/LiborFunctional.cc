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



#include "LiborFunctional.h"
#include "Option.h"
#include "Bond.h"


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(LiborFunctional)
	
	
Real 
H_pq(int p, int q, const RealArray1D& H, Real delta)
{
	Real sum=0.0;
	for(int j=p;j<q;j++) sum+=H[j+1];
			
	return delta*sum;
}


Real 
swapRate(int p, int q, const RealArray1D& H, Real delta)
{
	Real Hpq=0.0;
	for(int j=p;j<q;j++) Hpq+=H[j+1];
	Hpq*=delta;
			
	return (H[p]-H[q])/Hpq;
}


Real 
forwardSwaptionPayoff(int p, int q, Real kappa, const RealArray1D& H, Real delta)
{
	Real Hpq=H_pq(p,q,H,delta);
    Real Spq=swapRate(p,q,H,delta);
	return Spq>kappa ? (Spq-kappa)*Hpq : 0.0;
}


Real 
forwardCapletPayoff(int i, Real kappa, const RealArray1D& H, Real delta)
{
	Real XiTi=X(i,H);                 // X_i(T_i)
	return XiTi>delta*kappa ? (XiTi-delta*kappa)*H[i+1] : 0.0;
}


Real 
forwardBondCallPayoff(BondCall* bc, const RealArray1D& H)
{
	Real  K = bc->getStrike();
	Bond* B = bc->getBond();
	
	int t = bc->getPeriodsToExpiry(),     // option expires at T_t 
	    p = B->get_p(),                   // coupon period [T_p,T_q)
	    q = B->get_q();                   
	const RealArray1D& c = B->get_c();    // coupons
	
	Real F_B = 0.0,                         // forward bond price
	     F_K = K*H[t];                      // forward strike price
	for(int i=p;i<q;i++) F_B+=c[i]*H[i];
		
	return F_B-F_K>0? F_B-F_K : 0.0;
}
		
	
	

MTGL_END_NAMESPACE(LiborFunctional)
MTGL_END_NAMESPACE(Martingale)
