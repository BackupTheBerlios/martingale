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

#include "Node.h"
#include "LatticeOption.h"
#include "LmmLattice.h"
#include "LowFactorDriftlessLMM.h"
#include "Derivatives.h"
#include <math.h>

using namespace Martingale;




/**********************************************************************************
 *
 *            LATTICE OPTION
 *
 *********************************************************************************/
 

Real 
LatticeOption::
forwardPrice() 
{
	computeConditionalExpectations();
	Node* root=theLattice->getRoot();
	return root->getPi();
}
		

	
void 
LatticeOption::
computeConditionalExpectations()
{
	// run through the list of nodes at time t=s and set the known value pi_s=h
	std::list<Node*>& nodes_s=theLattice->getNodeList(s);	    
    std::list<Node*>::const_iterator theNode;
	for(theNode=nodes_s.begin(); theNode!=nodes_s.end(); ++theNode)
	{
		// the list is a list of pointers to Node, so *theNode is pointer to Node
		Node* node=*theNode;      
		node->setPi(forwardPayoff(node));
	}
	
	// backward computation through earlier nodes
	for(int t=s-1; t>=0;t--){
		
		std::list<Node*>& nodes=theLattice->getNodeList(t);	    
    	for(theNode=nodes.begin(); theNode!=nodes.end(); ++theNode)
	    {
		     Node* node=*theNode;
			 std::list<Edge>& edges=node->getEdges();
				
			 Real E_th=0.0;
			 std::list<Edge>::const_iterator theEdge;
			 for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge)
	         {	
		         Real p=theEdge->probability;
				 // E_{t+1}(h) at the node the edge points to.    
				 Real E_t1h=theEdge->node->getPi(); 
				 E_th+=p*E_t1h;
	         }	
			 node->setPi(E_th);			 
	    } // end for nodes
			
	} // end for t
} // end computeConditionalExpectations
	


/**********************************************************************************
 *
 *                     SWAPTIONS
 *
 *********************************************************************************/


Real 
LatticeSwaption::
forwardPayoff(Node* node) 
{
	Real swapRate=node->swapRate(p,q);
	if(swapRate<=kappa) return 0.0;
	// else, forward price of the annuity B_{p,q}
    Real Hpq=node->Hpq(p,q); 
	return (swapRate-kappa)*Hpq;
}

	

void 
LatticeSwaption::
test(int s, int p, int q)
{
	Timer watch; watch.start();
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(q,VolSurface::CONST,Correlations::CS); 
	// number of time steps in each Libor accrual interval
	int nSteps=2;
	ConstVolLmmLattice3F theLattice(fl,s,nSteps);
	LmmNode* root=theLattice.getRoot();
	Real strike=root->swapRate(p,q);           // swap rate at time zero
				
	LiborMarketModel* lmm=new LowFactorDriftlessLMM(fl,3);
	Derivative* swpnLmm=new Swaption(p,q,s,strike,lmm);
	Real mcPrice=swpnLmm->monteCarloForwardPrice(20000),
	     aPrice =swpnLmm->analyticForwardPrice();
		
	cout << "\n\n\nSwaption forward price: "
	     << "\nAnalytic: " << aPrice
	     << "\nMonte Carlo, 3 factors: " << mcPrice;
		

	LatticeSwaption<LmmNode3F> swpn(&theLattice,s,p,q,strike,nSteps);
	Real treePrice=swpn.forwardPrice();
		
	cout << "\nTree: " << treePrice;
        
	watch.stop();
	watch.report("Time");
}
	

	
	

 
