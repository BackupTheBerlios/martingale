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

#ifndef martingale_latticeoption_h    
#define martingale_latticeoption_h

#include "LmmLattice.h"
#include  "LowFactorDriftlessLMM.h"
#include "Derivatives.h"
#include <math.h>

MTGL_BEGIN_NAMESPACE(Martingale)




/**********************************************************************************
 *
 *            MARTINGALE pi_t=E_t(h) in a LATTICE
 *
 *********************************************************************************/
 


/** <p>European option evaluated in a {@link Lattice}. We work with forward prices
 *  at the horizon T of the underlying asset price process instead of discounted prices
 *  since this works for all options including Libor derivatives.
 *
 *  <p>The option specifies a payoff h at discrete time t=s (continuous time \f$T_s\f$)
 *  accrued forward to the horizon T. With this the forward option price
 *  \f$\pi_t\f$ satisfies \f$\pi_s=h\f$ and \f$\pi_t=E_t(h)\f$ for \f$t<s\f$ and this 
 *  price can be computed recursively as
 *
 * \f[pi_t=E_t(\pi_{t+1}),\quad\quad \pi_s=h.\f]
 *
 * In this manner the price is computed backward from time t=s to time t=0 through all
 * nodes in the lattice. At each node the conditional expectation \f$E_t(\pi_{t+1})\f$
 * is computed by averaging the values of \f$\pi_{t+1}\f$ over all nodes which can be
 * reached from this node weighted according to the transition probabilities.
 * The class provides methods to roll back the price to the root node and report this value
 * as the forward price of the option.
 *
 * <p>We only need the type of the lattice nodes as a template parameter since this then 
 * dtermines the lattice type as Lattice<Node>.
 *
 * @param Node type of node used in the underlying lattice.
 */
template<typename Node>
class LatticeEuropeanOption {
	
	Lattice<Node>* theLattice; // underlying lattice
	
	int s;                     // time at which the values of the payoff h are known
	                           // typically time of option expiration
	
public:
	
	
	/** @param lattice the underlying {@link Lattice}.
	 *  @param t time at which the values of <code>this</code> are known.
	 */
	LatticeEuropeanOption(Lattice<Node>* lattice, int t) : theLattice(lattice), s(t) {  }
		
		
	/** The option payoff h at time t=s accrued forward to the horizon of the underlying
	 *  asset price process evaluated at all nodes at time t=s. Applies only to nodes at time s and this is 
     *  not checked. The time s when the option pays off is a property of this class.
	 */
	virtual Real forwardPayoff(Node* node) = 0;
		
		
	/** The forward price at time t=0.
	 */
	Real forwardPrice() 
	{
		computeConditionalExpectations();
		Node* root=theLattice->getRoot();
		return root->getPi();
	}
		

		
private:		
		
	/** This method writes the forward price \f$\pi_t=E_t(h)\f$ into all
	 *  nodes in the lattice at times \f$t<s\f$, where s is the time at which the values
	 *  of <code>this</code> are known. The value is written into the field <code>pi</code>
	 *  of each node.
	 */
	void computeConditionalExpectations()
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
	
	
}; // end LmmLatticeEuropeanOption 



/**********************************************************************************
 *
 *                     SWAPTIONS
 *
 *********************************************************************************/



/** <p>The forward payer swaption which can be exercised at time 
 *  \f$T_s\f$ into a payer swap on the interval \f$[T_p,T_q]\f$
 * evaluated in a {@link Libortree2F} lattice.
 *
 * @param Node must extend LmmNode.
 */
template<typename Node>
class LatticeSwaption : public LatticeEuropeanOption<Node> {
	
	int s,    // exercise time T_s
	    p,    // swap begins at time T_p
	    q;    // swap ends at time T_q
	
	Real kappa;  // the strike rate
	
public:
	
	/** @param lattice the underlying LmmLattice.
	 *  @param _s swaption is exercisable at time T_s, s=s0.
	 *  @param _p  swap begins at time T_p, p=p0.
	 *  @param _q  swap ends at time T_q, q=q0.
	 */
	LatticeSwaption(Lattice<Node>* lattice, int _s, int _p, int _q, Real strike) :
	LatticeEuropeanOption<Node>(lattice,_s),
	p(_p), q(_q),
	kappa(strike)
    {  }
	
	// the NodeType extends LmmNode
	Real forwardPayoff(Node* node) 
    {
		Real swapRate=node->swapRate(p,q);
		if(swapRate<=kappa) return 0.0;
		// else, forward price of the annuity B_{p,q}
	    Real Hpq=node->Hpq(p,q); 
		return (swapRate-kappa)*Hpq;
	}

	
// TEST
	
	/** Allocate a sample LMMTree2F of dimension q and compute the forward price
	 *  of the at the money payer swaption exercisable at time \f$T_s\f$ into a swap 
	 *  along \f$[T_p,T_q]\f$.
	 */
	static void test(int s, int p, int q)
    {
		Timer watch; watch.start();
		ConstVolLiborFactorLoading* fl=ConstVolLiborFactorLoading::sample(q);
		ConstVolLmmLattice2F theLattice(fl,s);
		LmmNode* root=theLattice.getRoot();
		
		Real strike=root->swapRate(p,q);
		LatticeSwaption<LmmNode2F> swpn(&theLattice,s,p,q,strike);
		Real treePrice=swpn.forwardPrice();
		
		cout << "\n\n\nSwaption forward price: " 
		     << "\nTree: " << treePrice;
		
		LiborMarketModel* lmm=new LowFactorDriftlessLMM(fl,2);
		Derivative* swpnLmm=new Swaption(p,q,s,strike,lmm);
		Real mcPrice=swpnLmm->monteCarloForwardPrice(10000);
		
		cout << "\nMonte Carlo: " << mcPrice;
        
		watch.stop();
		watch.report("Time");
	}
	
}; // end LatticeSwaption


typedef LatticeSwaption<LmmNode2F> LatticeSwaption2F;
	
	
	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
