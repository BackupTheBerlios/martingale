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


MTGL_BEGIN_NAMESPACE(Martingale)



// forward declaration
template<typename Node> 
class Lattice;

class LmmNode2F;
class LmmNode3F;             // for typedef LatticeSwaption<LmmNode3F> below




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
 * @param Node type of node used in the underlying lattice.
 */
template<typename Node>
class LatticeEuropeanOption {
	
	Lattice<Node>* theLattice; // underlying lattice
	
	int s;                     // time at which the values of the payoff h are known
	                           // typically time of option expiration
	
public:
	
	
	/** @param lattice the underlying {@link Lattice}.
	 *  @param t number of time steps to option exercise.
	 */
	LatticeEuropeanOption(Lattice<Node>* lattice, int t) : theLattice(lattice), s(t) {  }
		
		
	/** The option payoff h at time t=s accrued forward to the horizon of the underlying
	 *  asset price process evaluated at all nodes at time t=s. Applies only to nodes at time s and this is 
     *  not checked. The time s when the option pays off is a property of this class.
	 */
	virtual Real forwardPayoff(Node* node) = 0;
		
		
	/** The forward price at time t=0.
	 */
	Real forwardPrice();
		

		
private:		
		
	/** Rolls back the forward price \f$\pi_t=E_t(h)\f$ starting from time t=s to time t=0
	 *  through all nodes in the lattice. Here s is the discrete time of option expiration.
	 *  The value \f$\pi_t=E_t(h)\f$ is written into the field <code>pi</code> of each node.
	 */
	void computeConditionalExpectations();	
	
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
 * @param Node must extend LmmNode.
 */
template<typename Node>
class LatticeSwaption : public LatticeEuropeanOption<Node> {
	
	int s,          // swaption exercises at time T_s
	    nSteps,     // number of time steps in each Libor accrual interval
	    p,          // swap begins at time T_p
	    q;          // swap ends at time T_q
	
	Real kappa;     // the strike rate
	
public:
	
	/** @param lattice the underlying LmmLattice.
	 *  @param s swaption is exercisable at time T_s, s=s0.
	 *  @param p_ swap begins at time T_p, p=p0.
	 *  @param q_  swap ends at time T_q, q=q0.
	 *  @param steps number of time steps in each Libor accrual interval.
	 */
	LatticeSwaption(Lattice<Node>* lattice, int s, int p_, int q_, Real strike, int steps=1) :
	LatticeEuropeanOption<Node>(lattice,s*steps),          // m=s*nSteps time steps needed until time T_s 
	nSteps(steps),
	p(p_), q(q_),
	kappa(strike)
    {  }
	
	/** Payoff of the swaption at a node living at the time of swaption exercise
	 *  accrued forward to the horizon \f$T_n\f$ of the underlying Libor process.
	 */ 
	//the NodeType extends LmmNode
	Real forwardPayoff(Node* node);

	
// TEST
	
	/** Allocate a sample {@link LMMlattice3F} of dimension q and computes the forward 
	 *  price of the at the money payer swaption exercisable at time \f$T_s\f$ into a swap 
	 *  along \f$[T_p,T_q]\f$.
	 */
	static void test(int s, int p, int q);
	
}; // end LatticeSwaption




/** Swaptions in two and three factor lattices for the {@link DriftlessLMM}.
 */
typedef LatticeSwaption<LmmNode2F> LatticeSwaption2F;
typedef LatticeSwaption<LmmNode3F> LatticeSwaption3F;	
	
	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
