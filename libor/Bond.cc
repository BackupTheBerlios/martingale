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

#include "Bond.h"
#include "LiborMarketModel.h"


MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;



// CONSTRUCTOR

Bond::
Bond(int k, int m, const RealArray1D& d, LiborMarketModel* lmm) :
LMM(lmm),
n(lmm->getDimension()), p(k), q(m), 
c(d),
b(q-p,p)
{   
	c.setIndexBase(p);
	b[p]=c[p];
	for(int j=p;j<q-1;j++) b[j+1]=b[j]+c[j+1];
}
	
	
// zero coupon bond maturing at T_i
Bond::
Bond(int i, LiborMarketModel* lmm) :
LMM(lmm),
n(lmm->getDimension()), p(i), q(i+1), 
c(1,p),
b(1,p)
{   
	c[p]=b[p]=1.0;
}

    
// PPRICE
	
//Forward price F(T_t)=\sum_{j=p}^{q-1}c_jH_j(T_t)
Real 
Bond::
forwardPrice(int t) const
{
     Real sum=0;
  	 for(int j=p;j<q;j++){ sum+=c[j]*(LMM->H_it(j,t)); }
     return sum;	
}   

	 
// Forward price F(0) at time 0
Real 
Bond::
forwardPrice() const 
{ 
	return forwardPrice(0); 
}

	 
// Cash price B(T_t)=\sum_{j=p}^{q-1}c_jB_j(T_t)
Real 
Bond::
cashPrice(int t) const 
{ 
	return forwardPrice(t)/(LMM->H_ii(t)); 
}

	 
// Cash price B(0) at time 0
Real 
Bond::
cashPrice() const 
{ 
	return forwardPrice(0)/(LMM->H0()); 
}
	 
	 
// PRINTING 

ostream& 
Bond::
printSelf(ostream& os) const
{
    return 
	os << "Bond along [T_"<<p<<",T_"<<q<<"], n="<<n<< endl
	   << "Coupons: " << c;
}



	
MTGL_END_NAMESPACE(Martingale)


