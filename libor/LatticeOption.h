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

// definitions are in header
#include "TypedefsMacros.h"
#include "LatticeOption.h"
#include "Node.h"
#include "LmmLattice.h"
#include "LowFactorDriftlessLMM.h"
#include "Derivatives.h"
#include <math.h>
#include <vector>

MTGL_BEGIN_NAMESPACE(Martingale)


using std::vector;


/**********************************************************************************
 *
 *            MARTINGALE pi_t=E_t(h) in a LATTICE
 *
 *********************************************************************************/
 


/** <p>European option evaluated in a {@link Lattice}. We work with forward prices
 *  at the horizon T of the underlying asset price process instead of discounted prices
 *  since this works for all options including Libor derivatives.
 *
 *  <p>The option specifies a payoff h at discrete time t=s (after s time steps)
 *  accrued forward to the horizon T. With this the forward option price
 *  \f$\pi_t\f$ satisfies \f$\pi_s=h\f$ and \f$\pi_t=E_t(h)\f$ for \f$t<s\f$ and this 
 *  price can be computed by backward recursion
 *
 * \f[pi_t=E_t(\pi_{t+1}),\quad\quad \pi_s=h.\f]
 *
 * starting from the time t=s of the option payoff to time t=0 through all
 * nodes in the lattice. Here discrete time t is an integer and denotes the time reached 
 * after time step number t as usual.
 * 
 * <p>At each node the conditional expectation \f$E_t(\pi_{t+1})\f$
 * is computed by averaging the values of \f$\pi_{t+1}\f$ over all nodes which can be
 * reached from this node weighted according to the transition probabilities.
 * The class provides methods to roll back the price to the root node and report this value
 * as the forward price of the option.
 *
 * <p>We only need the type of the lattice nodes as a template parameter since this then 
 * dtermines the lattice type as Lattice<Node>.
 *
 * @param Lattice type of the underlying lattice.
 */
template<typename Lattice>
class LatticeEuropeanOption {
	
protected:
	
	 /** the underlying lattice. */
	Lattice* theLattice;      
	
	/** the type of nodes in the lattice. */
	typedef typename Lattice::NodeType NodeType;
	
	/** time step at which the option payoff is known (typically expiration).*/
	int expiration;            
	                           
	
public:
	
	
	/** The underlying lattice. */
	Lattice* getLattice(){ return theLattice; }
	
	
	/** @param lattice the underlying {@link Lattice}.
	 *  @param t number of time steps to option expiration.
	 */
	LatticeEuropeanOption(Lattice* lattice, int t) : 
    theLattice(lattice), expiration(t) 
    {  }
	
	~LatticeEuropeanOption(){ delete theLattice; }
		
		
	/** The option payoff h at time t=s accrued forward to the horizon of the underlying
	 *  asset price process evaluated at all nodes at time t=s. Applies only to nodes at time s and this is 
     *  not checked. The time s when the option pays off is a property of this class.
	 */
	virtual Real forwardPayoff(NodeType* node) = 0;
		
		
	/** The forward price at time t=0.
	 */
	Real forwardPrice()
    {
	    computeConditionalExpectations();
	    NodeType* root=theLattice->getRoot();
	    return root->getPi();
    }	

		
private:		
		
/** Rolls back the forward price \f$\pi_t=E_t(h)\f$ starting from time t=s to time t=0
 *  through all nodes in the lattice. Here s is the discrete time of option expiration.
 *  The value \f$\pi_t=E_t(h)\f$ is written into the field <code>pi</code> of each node.
 */
void computeConditionalExpectations()
{
	// run through the list of nodes at time t=expiration and set the known value pi_t=h
	vector<NodeType*>* nodes_s=theLattice->getNodeList(expiration);
	// *theNode is pointer to Node
    vector<NodeType*>::const_iterator theNode=nodes_s->begin();
	while(theNode!=nodes_s->end()) {
	
		NodeType* node=*theNode;
		node->setPi(forwardPayoff(node));
		theNode++;
	}
	
	// backward computation through earlier nodes
	for(int t=expiration-1; t>=0;t--){
		
		vector<NodeType*>* nodes_t=theLattice->getNodeList(t);
	    theNode=nodes_t->begin();
    	while(theNode!=nodes_t->end()) {
	    
			NodeType* node = *theNode; 
			vector<Edge*>* edges=node->getEdges();
				
			 Real E_th=0.0;
			 vector<Edge*>::const_iterator theEdge=edges->begin();
			 while(theEdge!=edges->end()) {
	        
		         Edge* edge = *theEdge;
				 Real p=edge->probability;
				 // E_{t+1}(h) at the node the edge points to.    
				 Real E_t1h=edge->node->getPi(); 
				 E_th+=p*E_t1h;
				 ++theEdge;
	         }	
			 node->setPi(E_th);	
			 ++theNode;
	    } // end while(theNode)
			
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
 * evaluated in a lattice for the driftless Libor market model
 * {@link DriftlessLMM}.
 *
 * @param Lmm_Lattice the type of lattice used to price the swaption.
 */
template<typename Lmm_Lattice>
class LatticeSwaption : public LatticeEuropeanOption<Lmm_Lattice> {
	
	int s,          // swaption exercises at time T_s
	    p,          // swap begins at time T_p
	    q;          // swap ends at time T_q
	
	Real kappa;     // the strike rate
	
	
public:
	
	/** The strike rate.*/
	Real getStrike(){ return kappa; }
	
	/** @param lattice the underlying LmmLattice.
	 *  @param s swaption is exercisable at time T_s, s=s0.
	 *  @param p_ swap begins at time T_p, p=p0.
	 *  @param q_  swap ends at time T_q, q=q0.
	 */
	LatticeSwaption
	(Lmm_Lattice* lattice, int s, int p_, int q_, Real strike) :
	// m=s*nSteps time steps needed until T_s = LatticeEuropeanOption::expiration
	LatticeEuropeanOption<Lmm_Lattice>(lattice,s*(lattice->getData()->nSteps)),          
	p(p_), q(q_),
	kappa(strike)
    {  }
	
	/** Payoff of the swaption at a node living at the time of swaption exercise
	 *  accrued forward to the horizon \f$T_n\f$ of the underlying Libor process.
	 */ 
	//the NodeType extends LmmNode
	Real forwardPayoff(NodeType* node){ return node->forwardSwaptionPayoff(p,q,kappa); }
	
// SAMPLE AND TEST
	
/** Sample at the money swaption in dimension q on accrual interval [T_p,T_q]
 *  with nSteps time steps in each accrual period. 
 *  @apram verbose details on lattice and Libor factorloading.
 */	
static LatticeSwaption* 
sample(int p, int q, int nSteps, bool verbose)
{
    Lmm_Lattice* lattice = Lmm_Lattice::sample(q,p,nSteps,verbose);
	NodeType* root=lattice->getRoot();
	Real strike=root->swapRate(p,q);           // swap rate at time zero	
	
	return new LatticeSwaption(lattice,p,p,q,strike);
}
	
	
/** Allocates a sample at the money LatticeSwaption along \f$[T_p,T_q]\f$
 *  exercisable at time \f$T_p\f$ and computes the forward price both in
 *  the lattice and in a driftless Libor market model with the same number of 
 *  factors.
 * @param nSteps time steps in each Libor compounding period.
 * @apram verbose details on lattice and Libor factorloading.
 */
static void test(int p, int q, int nSteps, bool verbose=false)
{	
	cout << "\n\nComputing swaption price:";
	
	Timer watch; watch.start();
    LatticeSwaption* latticeSwaption = sample(p,q,nSteps,verbose);
	Real latticePrice=latticeSwaption->forwardPrice();
	cout << "\n\nLattice price without volatility rescaling: " << latticePrice;
	watch.stop();
	watch.report("Time");

	// lattice price with volatility rescaled
	latticeSwaption->getLattice()->rescaleVols();
	latticePrice=latticeSwaption->forwardPrice();
	cout << "\n\nLattice price with volatility rescaling: " << latticePrice;
	watch.stop();
	watch.report("Time");
		
	// LMM and Monte carlo
	int r = latticeSwaption->getLattice()->nFactors();    // number of factors
	LiborFactorLoading* fl = latticeSwaption->getLattice()->getFactorLoading();
	Real strike = latticeSwaption->getStrike();
		
	delete latticeSwaption;
	
	watch.start();
	LiborMarketModel* lmm=new LowFactorDriftlessLMM(fl,r);
	Derivative* swpnLmm=new Swaption(p,q,p,strike,lmm);
	Real aPrice =swpnLmm->analyticForwardPrice();
		
	cout << "\nAnalytic: " << aPrice
	     << "\nMonte Carlo price, 20000 paths: ";
		 
	Real mcPrice=swpnLmm->monteCarloForwardPrice(20000);	
	cout << mcPrice;
	
	watch.stop();
	watch.report("Time");

}	
	
	
}; // end LatticeSwaption



/** Swaption in two factor lightweight Lmmlattice.*/
typedef LatticeSwaption<LmmLattice2F> LatticeSwaption2F;
/** Swaption in three factor lightweight Lmmlattice.*/
typedef LatticeSwaption<LmmLattice3F> LatticeSwaption3F;




MTGL_END_NAMESPACE(Martingale)

#endif
 
