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
#include "Array.h"
#include <list>

MTGL_BEGIN_NAMESPACE(Martingale)

using std::list;

// we are using
class LiborFactorLoading;             // LiborFactorLoading.h


/*! \file Node.h
 *  Nodes in a stochastic lattice in general and in lattices for the 
 *  driftless Libor market model {@link DriftlessLMM} and baskets of assets 
 *  with constant volatility in particular. See book, chapter 8.
 *
 * <p><a name="lmm-node">Nodes for a general LMM lattice:</a>
 *  A node in a lattice for the Libor market model lives at a time in some
 *  accrual interval \f$(T_{t-1},T_t]\f$. At this time the surviving Libors are
 *  \f$X_t,X_{t+1},\dots,X_{n-1}\f$, where n is the number of accrual intervals.
 *  The node however stores the vector of accrual factors
 *  \f[H=(H_t,...,H_n)\f]
 *  as the most suitable way to store the state of the underlying Libor process since
 *  all other quantities can be recovered from these with minimal computational effort.
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
 * (i,j,k) with
 * \f[Z_1=ia,\quad Z_2=ja,\quad Z_3=ka.\f]
 * Two factor nodes are completely similar but the implementation is kept separate
 * for greater clarity.
 */
 
 

/**********************************************************************************
 *
 *            PROBABILITY WEIGHTED EDGES
 *
 *********************************************************************************/ 
 

class Node;

/** Edge in a node of a general stochastic lattice. Provides pointer to the node
 *  the edge connects to and transition probability along the edge.
 */
struct Edge {
	
	 /** Node we transition to along this edge */
	 Node* node; 
	 /** Transition probability along this edge */
	 Real probability;        

}; // end Edge	
	


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
	list<Edge> edges;  
		                    
		
public:
	
// ACCESSORS
	
	/** Discrete time t at which the node lives. */
	int getTimeStep() const { return s; }
	Real getPi() const { return pi; }
	void setPi(Real x){ pi=x; }
	
	/** List of edges */
	list<Edge>& getEdges() { return edges; }
	
	
// CONSTRUCTOR
	
	/** @param s_ number of time steps to reach the node from time zero.
	 */
	Node(int s_) : s(s_), edges() {  }
			 		
    ~Node(){ getEdges().~list<Edge>(); }
		

	
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


/** <p>Lightweight nodes in an {@link LmmLattice}.
 *  Stores only the two state variables from which the accrual factors
 *  \f$H_j\f$ are computed. Here we store as little as possible in each node.
 *
 * <p>The only service provided is the computation of the vector H
 * of accrual factors at this node. This vector is the most efficient way
 * to store the state of the Libor process at this node.
 *
 * <p>The other functionals of the Libor process such as 
 * swaprates, swaption payoffs, caplet payoffs, annuities, bonds,...
 * are computed from the vector H by global functions which take H as a 
 * parameter. Reason: avoid the repeated computation of H.
 */
class LiteLmmNode : public Node {

public:	
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a of the driving Brownian motion Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
LiteLmmNode(int s, const IntArray1D& k, LmmLatticeData* latticeData);

~LiteLmmNode(){ delete[] k_; }


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

/** Prints what type of node it is.*/
static std::ostream& printType(std::ostream& os);

/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os);


private:

	static RealArray1D H_;           // workspace for accrual factors H_j
	static RealArray1D V_;           // workspace for volatility parts V_j of the log(U_j)
	                                 // book, 8.1

	LmmLatticeData* lattice;         // information about the LmmLattice the node lives in.
	int* k_;                         // Z_j=k[j]*a, state of the Brownian driver
	                                 // a=sqrt(dt) the tick size of a standard Brownian 
	                                 // motion over an interval of length dt.
	
/** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
 */
int get_t() const;


}; // end LiteLmmNode
	
	

/** Heavyweight nodes in an {@link LmmLattice}.
 *  Stores the entire vector \f$H=(H_j,\dots,H_n)\f$ of accrual factors
 *  still alive at the time at which the node lives. The values are computed
 *  in the constructor of the node.
 *
 * <p>The only service provided is the computation of the vector H
 * of accrual factors at this node. This vector is the most efficient way
 * to store the state of the Libor process at this node.
 *
 * <p>The other functionals of the Libor process such as 
 * swaprates, swaption payoffs, caplet payoffs, annuities, bonds,...
 * are computed from the vector H by global functions which take H as a 
 * parameter. Reason: conformity with the case of LiteLmmNodes.
 */
class HeavyLmmNode : public Node {

public:
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a, a=sqrt(dt), of the Brownian driver Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
HeavyLmmNode(int s, const IntArray1D& k, LmmLatticeData* latticeData);	
~HeavyLmmNode(){ delete[] k_; }


/** The state Z_j=k[j]*a of the Brownian driver at this node.
 *  Here a=sqrt(dt) is the ticksize of a standard Brownian motion Z over a time step 
 *  of size dt.
 */
int* getIntegerTicks(){ return k_; }
	
/** The vector \f$H=(H_p,\dots,H_n)\f$ of accrual factors at the node.
 *  Natural indices j=p,...,n.
 */
const RealArray1D& Hvect(int p){ return H_; }

/** Prints what type of node it is.*/
static std::ostream& printType(std::ostream& os);

/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os) const;


private:

	LmmLatticeData* lattice;         // information about the LmmLattice the node lives in.
	int* k_;                           // Z_j=k[j]*a, a=sqrt(dt) the tick size of a standard Brownian 
	                                   // motion over an interval of length dt.
	RealArray1D H_;                    // the vector H=(H_j,...,H_n) of accrual factors still alive
	                                   // at the time at which the node lives.
    static RealArray1D V_;             // workspace for the volatility parts V_j of the log(U_j),
                                       // book 8.1

/** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
 */
int get_t() const;


}; // end HeavyLmmNode
		
	
	

/**********************************************************************************
 *
 *            NODE IN GENERAL BASKET
 *
 *********************************************************************************/

// forward declaration
class BasketLatticeData;


/** <p>Lightweight nodes in a {@link BasketLattice}.
 *  Stores only the r state variables from which the assets
 *  \f$S_j\f$ are computed. Here we store as little as possible in each node.
 *
 * <p>The only service provided is the computation of the vector S
 * of accrual factors at this node.
 */
class LiteBasketNode : public Node {

public:	
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a of the driving Brownian motion Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
LiteBasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData);
~LiteBasketNode(){ delete[] k_; }
	
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

/** Prints what type of node it is.*/
static std::ostream& printType(std::ostream& os);

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


}; // end LiteBasketNode


	
	

/** Heavyweight node in a {@link BasketLattice}.
 *  Stores the entire asset vector \f$S=(S_1,\dots,S_n)\f$ at the node. 
 *  The values are computed in the constructor of the node.
 */
class HeavyBasketNode : public Node {

public:
/** 
 *  @param s number of time steps to reach the node from time zero.
 *  @param k state Z_j=k[j]*a, a=sqrt(dt), of the Brownian driver Z.
 *  @param info object encapsulating information about the lattice the node lives in.
 */
HeavyBasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData);	
~HeavyBasketNode(){ delete[] k_; }
	
	
/** The vector of assets \f$S=(H_p,\dots,H_n)\f$ at the node.
 *  Natural indices j=p,...,n.
 */
const RealArray1D& Svect(int p){ return S_; }

/** The state Z_j=k[j]*a of the Brownian driver at this node.
 *  Here a=sqrt(dt) is the ticksize of a standard Brownian motion Z over a time step 
 *  of size dt.
 */
int* getIntegerTicks(){ return k_; }

/** Prints what type of node it is.*/
static std::ostream& printType(std::ostream& os);

/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os) const;


private:

    static RealArray1D V_;             // workspace for the volatility parts V_j of the log(S_j),
                                       // book 3.11
	BasketLatticeData* lattice;        // information about the Lattice the node lives in.
	int* k_;                           // Z_j=k[j]*a, a=sqrt(dt) the tick size of a standard Brownian 
	                                   // motion over an interval of length dt.
	RealArray1D S_;                    // the vector S=(S_0,...,S_{n-1}) of assets the node.

/** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
 */
int get_t() const;


}; // end HeavyBasketNode
		
	

	
MTGL_END_NAMESPACE(Martingale)

#endif
 
