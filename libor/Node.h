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


#ifndef martingale_node_h    
#define martingale_node_h

#include "TypedefsMacros.h"
#include "Utils.h"
#include "Edge.h"
#include "Array.h"                    // problem with typdefs in forward declarations
#include <vector>                     // direct member


MTGL_BEGIN_NAMESPACE(Martingale)

using std::vector;

// forward declarations
class LiborFactorLoading;             // LiborFactorLoading.h
class std::ostream;
class BondCall;
// class RealArray1D;


/*! \file Node.h
 *  Nodes in a stochastic lattice in general and in lattices for the 
 *  driftless Libor market model {@link DriftlessLMM} and baskets of assets 
 *  with constant volatility in particular. See book, chapter 8.
 *
 * <p><a name="lmm-node">Nodes for a general LMM lattice:</a>
 *  A node in a lattice for the Libor market model lives at a time in some
 *  accrual interval \f$(T_{t-1},T_t]\f$. At this time the surviving Libors are
 *  \f$X_t,X_{t+1},\dots,X_{n-1}\f$, where n is the number of accrual intervals.
 *  Nodes only store the state of the Brownian driver (see below) to obtain the 
 *  smallest possible memory footprint. This allows us to build lattices with
 *  millions of nodes. The vector of accrual factors
 *  \f[H=(H_t,...,H_n)\f]
 *  and all other functionals of the Libor process are computed from 
 *  this vector with minimal computational effort. See {@link LiborFunctional}.
 *
 * <a name="lmm-node-3f">Nodes for a 3 factor LMM lattice:</a>
 * <p>Nodes in a 3 factor lattice for the Libor market model {@link LmmLattice3F}
 *  compute the \f$H_j\f$ from the volatility parts \f$V_j\f$ of the forward transported
 *  Libors \f$U_j\f$ (book, 6.8) approximated as
 *  \f[V_j(s)\simeq\sigma_j\left[R_{j1}Z_1(s)+R_{j2}Z_2(s)+R_{j3}Z_3(s)\right],\f]
 *  See book, 8.1.2, for details and notation. 
 *
 * <p>The \f$Z_j(t)\f$ are independent standard Brownian motions which evolve from the state
 * \f$Z_j(0)=0\f$ and then tick up or down in ticks of size \f$a=\sqrt{dt}\f$ where dt is the 
 * size of the time step. The state at any node is then given by the triple of integers
 * (k0,k1,k2) with
 * \f[Z_1=k0*a,\quad Z_2=k1*a,\quad Z_3=k2*a.\f]
 * Each node stores the integer vector k as the basic state.
 *
 * <p>The case of nodes for an asset basket is similar except that the fundamental
 * vector which is computed is the vector S of asset prices at this node. 
 */
 


/**********************************************************************************
 *
 *            NODES A STOCHASTIC LATTICE
 *
 *********************************************************************************/



/** Node in a general stochastic lattice. Allocates empty list of edges.
 *  The edges have to be filled in in the derived node classes.
 *  More details at Node.h.
 */
class Node {
	
protected:
	
	/** number of time steps to reach the node. */
	int s; 

	/** Cache for the value \f$\pi_s=E_s(h)\f$ of some conditional expectation 
	 *  at this node. 
	 */
	Real  pi;               
	
	/** List of edges originating at this node. */
	vector<Edge*> edges;  
		                    
		
public:
	
// ACCESSORS
	
	/** Discrete time t at which the node lives. */
	int getTimeStep() const { return s; }
	Real getPi() const { return pi; }
	void setPi(Real x){ pi=x; }
	
	/** List of edges */
	vector<Edge*>* getEdges() { return &edges; }
	
	
// CONSTRUCTOR
	
	/** @param s_ number of time steps to reach the node from time zero.
	 */
	Node(int s_) : s(s_), edges() {  }
	
	// manual deallocation is necessary since we are dealing with
	// lists of pointers. Make sure the pointees are deallocated
	// before the pointers disappear
	~Node();
	
// DIAGNOSTIC
	
	/** Diagnostic. Prints the transition probabilities and corresponding Q_ij(t) values.
	 */
	void printTransitionProbabilities() const;
	 
	

}; // end Node
	


/**********************************************************************************
 *
 *            NODE IN GENERAL LMM LATTICE
 *
 *********************************************************************************/


// forward declaration
class LmmLatticeData;
	

/** <p>Node in an {@link LmmLattice}.
 *  Stores only the two state variables from which the accrual factors
 *  \f$H_j\f$ are computed. Smallest possible memory footprint.
 */
class LmmNode : public Node {
	
protected:
	
   LmmLatticeData* lattice;   // information about the LmmLattice the node lives in.

public:	
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a of the driving Brownian motion Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
LmmNode(int s, const IntArray1D& k, LmmLatticeData* latticeData);

virtual ~LmmNode(){ delete[] k_; }


/** The state Z_j=k[j]*a of the Brownian driver at this node.
 *  Here a=sqrt(dt) is the ticksize of a standard Brownian motion Z over a time step 
 *  of size dt.
 */
int* getIntegerTicks(){ return k_; }
	
/** <p>The vector \f$H=(H_p,\dots,H_n)\f$ of accrual factors at the node.
 *  Natural indices j=p,...,n.
 *  <p>This is a view of this vector in a static workspace (speed, no memory allocation).
 *  If the vector has to be preserved it must be copied before commands execute
 *  which overwrite the workspace.
 */
const RealArray1D& Hvect(int p);
	
/** The forward price \f$H_{p,q}\f$ of the annuity \f$B_{p,q}\f$ over the 
 *  interval [T_p,T_q] at this node.
 */
Real H_pq(int p, int q);
		
/** The swaprate \f$S_{p,q}\f$ for a swap on [T_p,T_q] at this node.
 */
Real swapRate(int p, int q);
		
/** Payoff of a forward swaption with strike rate kappa exercising into 
 *  a swap on the interval [T_p,T_q] if exercised at this node. 
 *  Payoff is compounded forward to the horizon T_n. 
 */
Real forwardSwaptionPayoff(int p, int q, Real kappa);

/** Forward compounded payoff of a caplet with strike rate kappa at this node. 
 *  Assumes that the node lives at the Libor reset point at which Libor for
 *  this caplet is set.
 */
Real forwardCapletPayoff(Real kappa);

/** Forward payoff of the BondCall bc if exercised at this node. */
Real forwardBondCallPayoff(BondCall* bc);

/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os);


private:

	static RealArray1D H_;           // workspace for accrual factors H_j
	static RealArray1D V_;           // workspace for volatility parts V_j of the log(U_j)
	                                 // book, 8.1

	int* k_;                         // Z_j=k[j]*a, state of the Brownian driver
	                                 // a=sqrt(dt) the tick size of a standard Brownian 
	                                 // motion over an interval of length dt.
	
/** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
 */
int get_t() const;


}; // end LmmNode
	
		
	
	

/**********************************************************************************
 *
 *            NODE IN GENERAL BASKET
 *
 *********************************************************************************/

// forward declaration
class BasketLatticeData;


/** <p>Node in a {@link BasketLattice}.
 *  Stores only the r state variables from which the assets
 *  \f$S_j\f$ are computed. Smallest possible memory footprint.
 *  The only service provided is the computation of the asset
 *  price vector S. 
 */
class BasketNode : public Node {

public:	
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a of the driving Brownian motion Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
BasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData);
~BasketNode(){ delete[] k_; }
	
/** <p>The vector of assets \f$S=(S_p,\dots,S_n)\f$ at the node.
 *  Natural indices j=p,...,n.
 *  <p>This is a view of this vector in a static workspace (speed, no memory allocation).
 *  If the vector has to be preserved it must be copied before commands execute
 *  which overwrite the workspace.
 */
const RealArray1D& Svect(int p);

/** The state Z_j=k[j]*a of the Brownian driver at this node.
 *  Here a=sqrt(dt) is the ticksize of a standard Brownian motion Z over a time step 
 *  of size dt.
 */
int* getIntegerTicks(){ return k_; }


/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os);


private:

	static RealArray1D S_;           // workspace for the assets H_j
	static RealArray1D V_;           // workspace for volatility parts V_j of the log(S_j)
	                                 // book, 3.11 and 8.1.
	BasketLatticeData* lattice;      // information about the Lattice the node lives in.
	int* k_;                         // Z_j=k[j]*a, state of the Brownian driver
	                                 // a=sqrt(dt) the tick size of a standard Brownian 
	                                 // motion over an interval of length dt.


}; // end BasketNode

	

	
MTGL_END_NAMESPACE(Martingale)

#endif
 
