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
#include "Matrices.h"
#include <list>

MTGL_BEGIN_NAMESPACE(Martingale)


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
	std::list<Edge> edges;  
		                    
		
public:
	
// ACCESSORS
	
	/** Discrete time t at which the node lives. */
	int getTimeStep() const { return s; }
	Real getPi() const { return pi; }
	void setPi(Real x){ pi=x; }
	
	/** List of edges */
	std::list<Edge>& getEdges() { return edges; }
	
	
// CONSTRUCTOR
	
	/** @param s_ number of time steps to reach the node from time zero.
	 */
	Node(int s_) : s(s_), edges() {  }
				 		
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


/** Single <a href="lmm-node">node</a> in an {@link LmmLattice} containing the vector H
 *  of accrual factors \f$H_j\f$ which  survive to the time at which the node lives.
 *  These are used as the most efficient means to store the state of the underlying 
 *  Libor process. Contains methods to compute Libors, swap rates, bonds, etc from these. 
 *  The mechanics of evolving the \f$H_j\f$ is left to derived classes.
 *  For more details see the file reference for the file Node.h (click on
 *  "File List").
 *  
 */
class LmmNode : public Node {
	
protected:
	
	/** Factor loading of the underlying lMM. */
	LiborFactorLoading* factorLoading;
	
	/** Dimension of underlying LMM. */
	int n;
	
	/** Number of time steps in each Libor accrual interval. */
	int nSteps;
	
	/** The accrual factors H_j(T_t), j=t,...,n at this node,
	 *  natural indexation, index base t.
	 */
	RealArray1D H; 

	
public:
	
// ACCESSORS

    /** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
	 */
	int get_t() const
    {
		if(s%nSteps==0) return s/nSteps;
		return s/nSteps+1;
	}		
	
	RealArray1D& getH(){ return H; }
	
	
// CONSTRUCTOR
	
	/** @param fl factor loading of underlying Libor process.
	 *  @param s number of time steps to reach the node from time zero.
	 *  @param steps number of time steps in each Libor accrual interval.
	 */
	LmmNode(LiborFactorLoading* fl, int s, int steps=1);
				 
	

	
// DIAGNOSTIC

	
	/** Diagnostic. Prints the time t, vector H and field pi.
	 */
	void printState() const;
	
	/** What type of LMM node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);

	
	
// LIBORS, SWAPRATES, ANNUITY (PBV)
	
	/** Libor \f$X_j=\delta_jL_j\f$ at this node.
	 */
	Real X(int j) const;

	
	/** Forward price of the annuity (PBV)
	 *  \f$H_{p,q}=B_{p,q}/B_n=\sum\nolimits_{j=p}^{q-1}\delta_jH_{j+1}\f$
	 *  at this node.
	 */
	Real Hpq(int p, int q) const;
	
	
	/** Forward swaprate \f$S_{p,q}\f$ at this node, swap interval
	 *  \f$[T_p,T_q]\f$, where \f$t\leq p<q\f$. Here t is the discrete time
	 *  at which the node lives.
	 */	  
	Real swapRate(int p, int q) const { return (H[p]-H[q])/Hpq(p,q); }
		

// ZERO COUPON BONDS
		
    /** The zero coupon bond \f$B_i\f$ expiring at time \f$T_i\f$ at this node. 
	 */
	Real B_i(int i) const
	{ 
		int t=get_t();
		return H[i]/H[t]; 
	}
	

}; // end LmmNode
	


/**********************************************************************************
 *
 *            NODE IN 2 FACTOR LMM LATTICE
 *
 *********************************************************************************/




/** <p><a href="lmm-node-3f">Node</a> in a 2 factor lattice for the Libor market 
 *  model {@link LmmLattice2F}.
 *  For more details see the file reference for the file Node.h (click on
 *  "File List").
 */
class LmmNode2F : public LmmNode {
	
	int i,                  // Z1=ia
		j;                  // Z2=ja,
    // a=sqrt(dt) the tick size of a standard Brownian motion over an interval of length dt.
	
public:
		
// ACCESSORS
	
	int get_i() const { return i; }
	int get_j() const { return j; }
	
	
// CONSTRUCTOR
	
	/** @param s number of time steps to reach this node from time zero.
	 *  @param nSteps number of time steps in each Libor acrual interval.
	 *  @param k state Z_1=ka
	 *  @param l state Z_2=la
	 */
	LmmNode2F(LiborFactorLoading* fl, int s, int nSteps, int k, int l) : LmmNode(fl,s,nSteps),
    i(k), j(l)
	{	}
	
	/** What type of LMM node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);
	

}; // end LmmNode2F
	



/**********************************************************************************
 *
 *            NODE IN 3 FACTOR LMM LATTICE
 *
 *********************************************************************************/




/** <p><a href="lmm-node-3f">Node</a> in a 3 factor lattice for the Libor market model 
 *  {@link LmmLattice3F}.
 *  For more details see the file reference for the file Node.h (click on
 *  "File List").
 */
class LmmNode3F : public LmmNode {
	

	int i,                  // state Z1=ia
		j,                  // state Z2=ja
	    k;                  // state Z3=ka,  
    // a=sqrt(dt) the tick size of a standard Brownian motion over an interval of length dt.

	
public:
		
// ACCESSORS
	
	int get_i() const { return i; }
	int get_j() const { return j; }
	int get_k() const { return k; }
	
// CONSTRUCTOR
	
	/** @param s number of time steps to reach this node from time zero.
	 *  @param nSteps number of time steps in each Libor acrual interval.
	 *  @param p state Z1=pa.
	 *  @param q state Z2=qa.
	 *  @param r state Z3=ra.
	 */
	LmmNode3F(LiborFactorLoading* fl, int s, int nSteps, int p, int q, int r) : LmmNode(fl,s,nSteps),
    i(p), j(q), k(r)
	{	}
				 
	/** What type of LMM node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);
	

}; // end LmmNode3F
	

	

/**********************************************************************************
 *
 *            NODE IN GENERAL BASKET
 *
 *********************************************************************************/


/** Single <a href="lmm-node">node</a> in a basket of assets with constant volatilities.
 *  Contains the vector S of asset prices.
 *  For more details see the file reference for the file BasketLattice.h.
 *  
 */
class BasketNode : public Node {
	
protected:
	
	/** Number of assets */
	int n;
	
	/** Vector of asset prices, index base zero.
	 */
	RealArray1D S; 

	
public:
	
// ACCESSORS
	
	/** Vector of asset prices. */
	RealArray1D& getS() { return S; }
	
	
// CONSTRUCTOR
	
	/** @param n number of assets.
	 *  @param s number of time steps to reach the node from time zero.
	 */
	BasketNode(int _n, int _s): Node(s), n(_n), S(_n) { }
				 
	

	
// DIAGNOSTIC

	
	/** Diagnostic. Prints the time t, vector H and field pi.
	 */
	void printState() const;
		
	/** What type of LMM node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);
	 	

}; // end BasketNode
	



/**********************************************************************************
 *
 *            NODE IN 2 FACTOR BASKET 
 *
 *********************************************************************************/



/** Node in a 2 factor lattice for a basket of constant volatility assets.
 *  For more details see the file reference for the file BasketLattice.h.
 */
class BasketNode2F : public BasketNode {
	
	int i,                  // Z1=ia
		j;                  // Z2=ja,
    // a=sqrt(dt) the tick size of a standard Brownian motion over an interval of length dt.
	
public:
		
// ACCESSORS
	
	int get_i(){ return i; }
	int get_j(){ return j; }
	
	
// CONSTRUCTOR
	
	/** @param n number of assets.
	 *  @param s number of time steps to reach this node from time zero.  
	 *  @param k state Z_1=ka
	 *  @param l state Z_2=la
	 */
	BasketNode2F(int n, int s, int k, int l) : BasketNode(n,s),
    i(k), j(l)
	{	}
	
	/** What type of node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);
	

}; // end LmmNode2F
	



/**********************************************************************************
 *
 *            NODE IN 3 FACTOR BASKET LATTICE
 *
 *********************************************************************************/



/** Node in a 3 factor lattice for a basket of constant volatility assets.
 *  For more details see the file reference for the file BasketLattice.h.
 */
class BasketNode3F : public BasketNode {
	

	int i,                  // state Z1=ia
		j,                  // state Z2=ja
	    k;                  // state Z3=ka,  
    // a=sqrt(dt) the tick size of a standard Brownian motion over an interval of length dt.

	
public:
		
// ACCESSORS
	
	int get_i(){ return i; }
	int get_j(){ return j; }
	int get_k(){ return k; }
	
// CONSTRUCTOR
	
	/** @param n number of assets.
	 *  @param s number of time steps to reach this node from time zero.
	 *  @param p state Z1=pa.
	 *  @param q state Z2=qa.
	 *  @param r state Z3=ra.
	 */
	BasketNode3F(int n, int s, int p, int q, int r) : BasketNode(n,s),
    i(p), j(q), k(r)
	{	}
	
	/** What type of node it is (2 or 3 factors).*/
	static std::ostream& printType(std::ostream& os);
				 
	
}; // end BasketNode3F
	




MTGL_END_NAMESPACE(Martingale)

#endif
 
