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

#ifndef martingale_lmmlattice_europeanoption_h    
#define martingale_lmmlattice_europeanoption_h

#include "LiborTree2F.h"

MTGL_BEGIN_NAMESPACE(Martingale)




/**********************************************************************************
 *
 *            FUNCTIONAL OF A TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
 


/** <p>European Libor derivative evaluated in a {@link LiborTree2F} lattice.
 *  The option specifies a payoff h at discrete time t=s (continuous time \f$T_s\f$)
 *  accrued forward to the horizon of the underlying LMM. Thus the forward otpion price
 *  \f$\pi_t\f$ satisfies \f$\pi_s=h\f$ and \f$\pi_t=E_t(h)\f$ for \f$t<s\f$ and this 
 *  price can be computed recursively as
 *
 * \f[pi_t=E_t(\pi_{t+1}),\quad\quad \pi_s=h.\f]
 *
 * In this manner the price is computed backward from time t=s to time t=0 through all
 * nodes in the lattice. At each node the conditional expectation \f$E_t(\pi_{t+1})\f$
 * is computed by averaging the values of \f$\pi_{t+1}\f$ over all nodes which can be
 * reahed from this node weighted according to the transition probabilities.
 * The class provides methods to roll back the price to the root node and report this value
 * as the forward price of the option.
 */
class LmmLatticeEuropeanOption {
	
	LiborTree2F* theTree;    // underlying 2 factor LMM tree
	
	int s;                   // time at which the values h of the martingale are known
	
public:
	
	/** @param tree the underlying {@link LiborTree2F}.
	 *  @param t time at which the values of <code>this</code> are known.
	 */
	LmmLatticeEuropeanOption(LiborTree2F* tree, int t) : theTree(tree), s(t) {  }
		
		
	/** The option payoff h at time t=s accrued forward to the horizon of the underlying
	 *  LMM evaluated at all nodes at time t=s. Applies only to nodes at time s and this is 
     *  not checked. The time s when the option pays off is a property of this class.
	 */
	virtual Real forwardPayoff(LiborTree2F::Node* node) const = 0;
		
		
	/** The forward price at time t=0.
	 */
	Real forwardPrice()
	{
		computeConditionalExpectations();
		LiborTree2F::Node* root=theTree->getRoot();
		return root->pi;
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
		std::list<LiborTree2F::Node*>& nodes_s=theTree->getNodeList(s);	    
	    std::list<LiborTree2F::Node*>::const_iterator theNode;
		for(theNode=nodes_s.begin(); theNode!=nodes_s.end(); ++theNode)
		{
			// the list is a list of pointers to Node, so *theNode is pointer to Node
			LiborTree2F::Node* node=*theNode;      
			node->pi=forwardPayoff(node);
		}
		
		// backward computation through earlier nodes
		for(int t=s-1; t>=0;t--){
			
			std::list<LiborTree2F::Node*>& nodes_s=theTree->getNodeList(t);	    
	    	for(theNode=nodes_s.begin(); theNode!=nodes_s.end(); ++theNode)
		    {
			     LiborTree2F::Node* node=*theNode;
				 Matrix<LiborTree2F::Edge*> edges=node->edges;
				
				 Real E_th=0.0;
				 for(int k=-1;k<2;k++)
	         	 for(int l=-1;l<2;l++)
		         if((k!=0)||(l!=0))
		         {	
		         	 LiborTree2F::Edge* edge=edges(k,l);
			         Real p=edge->probability;
					 // E_{t+1}(h) at the node the edge points to.    
					 Real E_t1h=edge->node->pi;    
					 E_th+=p*E_t1h;
		         }	
				 node->pi=E_th;
		    } // end for nodes
		} // end for t
	} // computeConditionalExpectations
	
	
}; // end LmmLatticeEuropeanOption 



/**********************************************************************************
 *
 *                     SWAPTIONS
 *
 *********************************************************************************/



/** The forward payer swaption which can be exercised at time 
 *  \f$T_s\f$ into a payer swap on the interval \f$[T_p,T_q]\f$
 * evaluated in a {@link Libortree2F} lattice.
 */
class LatticeSwaption : public LmmLatticeEuropeanOption {
	
	int s,    // swaption is exercisable at time T_s
	    p,    // swap begins at time T_p
	    q;    // swap ends at time T_q
	
	Real kappa;  // the strike rate
	
public:
	
	/** @param tree the underlying LMMTree2F.
	 *  @param s0 swaption is exercisable at time T_s, s=s0.
	 *  @param p0  swap begins at time T_p, p=p0.
	 *  @param q0  swap ends at time T_q, q=q0.
	 */
	LatticeSwaption(LiborTree2F* tree, int s0, int p0, int q0, Real strike) :
	LmmLatticeEuropeanOption(tree,s),
	s(s0), 
	p(p0), q(q0),
	kappa(strike)
    {  }
	
	Real forwardPayoff(LiborTree2F::Node* node) const
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
		LiborFactorLoading* fl=JR_FactorLoading::sample(q);
		LiborTree2F theTree(fl,s);
		LiborTree2F::Node* root=theTree.getRoot();
		Real strike=root->swapRate(p,q);
		LatticeSwaption swpn(&theTree,s,p,q,strike);
		
		cout << "Swaption forward price = " << swpn.forwardPrice();
        
		watch.stop();
		watch.report("Time");
	}
	
}; // end LatticeSwaption
	
	
	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
