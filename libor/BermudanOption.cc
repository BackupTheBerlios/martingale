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
 


#include "BermudanOption.h"
#include "LmmLattice.h"
//#include "BasketLattice.h"
#include "FinMath.h"
#include "LiborMarketModel.h"          // Bond


MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::endl;



/*******************************************************************************
 *
 *                        BERMUDAN SWAPTION
 *
 ******************************************************************************/


BermudanSwaption::
BermudanSwaption
(int a, int b, int paths, Real strike, LiborMarketModel* lmm, bool verbose) :
// Libors L_j, j>=a needed until time b-1.
LiborDerivative
(lmm, new LiborPathsToTriggerTime(lmm,b-1,a),b-1,false,true,true,false),                                           
p(a), q(b), nPath(paths), kappa(strike)
{   
	LPG->setTrigger(new PjTrigger(lmm,this,verbose));
}


BermudanSwaption* 
BermudanSwaption::
sample(int p, int q, int paths, bool verbose, int lmmType, int volType, int corrType)
{
	LiborMarketModel* lmm=LiborMarketModel::sample(q,lmmType,volType,corrType);
	Real kappa=lmm->swapRate(p,q);
	return new BermudanSwaption(p,q,paths,kappa,lmm,verbose);
} 

	 
bool
BermudanSwaption::
isExercisable(Real t)
{
	const RealArray1D& T = LMM->getTenorStructure();
	for(int i=q-1;i>=p;--i) if(t==T[i]) return true;
	return false;
}


Real 
BermudanSwaption::
currentForwardPayoff(int s)
{
	if(s==q) return 0.0;                          // never exercised
	Real S_sqT,H_sqT,h;
	S_sqT=LMM->swapRate(s,q,s);                   // swaprate S_{s,q}(T_s)
	if(S_sqT<kappa) return 0.0;
    H_sqT=LMM->H_pq(s,q,s);
    return H_sqT*(S_sqT-kappa);                   // payoff at time T_s	
}


Real 
BermudanSwaption::
forwardPayoffAlongCurrentPath()
{
	int s=LPG->getTime();                         // exercise time	
    return currentForwardPayoff(s);
}


Real 
BermudanSwaption::
forwardPayoff(LmmNode* node)
{ 
	return node->forwardSwaptionPayoff(p,q,kappa); 
}
 

ostream& 
BermudanSwaption::
printSelf(ostream& os) const
{  return os << "\nBermudan swaption on [T_"<<p<<",T_"<<q<<"], strike rate "<< kappa; }

	


MTGL_END_NAMESPACE(Martingale)