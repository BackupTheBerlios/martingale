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

#include "LiborFactorLoading.h"
#include "Utils.h"
#include <math.h>
#include <list>


MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file Node.h
 *  Nodes in a stochastic lattice in general and in lattices 
 *  for the Libor market model in particular.
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
	int getTimeStep(){ return s; }
	Real getPi(){ return pi; }
	void setPi(Real x){ pi=x; }
	
	/** List of edges */
	std::list<Edge>& getEdges(){ return edges; }
	
	
// CONSTRUCTOR
	
	/** @param s_ number of time steps to reach the node from time zero.
	 */
	Node(int s_) : s(s_), edges() {  }
				 
		
    ~Node()
	{
		std::list<Edge>::const_iterator theEdge;         // pointer to Edge
		for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge)  delete &(*theEdge); 
	}
		

	
// DIAGNOSTIC
	
	/** Diagnostic. Prints the transition probabilities and corresponding Q_ij(t) values.
	 */
	void printTransitionProbabilities();
	 
	

}; // end Node
	


/**********************************************************************************
 *
 *            NODE IN GENERAL LMM LATTICE
 *
 *********************************************************************************/


/** Single node in an {@link LmmLattice} containing the vector \f$H=(H_t(T_t),...,H_n(T_t))\f$ 
 *  and methods to compute Libors, swaprates, bonds, etc from these. The mechanics of evolving the
 *  \f$H_j\f$ is left to derived classes.
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
	vector<Real> H; 

	
public:
	
// ACCESSORS

    /** Returns t such that the node lives in the accrual interval \f$(T_{t-1},T_t]\f$.
	 */
	int get_t()
    {
		if(s%nSteps==0) return s/nSteps;
		return s/nSteps+1;
	}		
	
	vector<Real>& getH(){ return H; }
	
	
// CONSTRUCTOR
	
	/** @param fl factor loading of underlying Libor process.
	 *  @param s number of time steps to reach the node from time zero.
	 *  @param steps number of time steps in each Libor accrual interval.
	 */
	LmmNode(LiborFactorLoading* fl, int s, int steps=1);
				 
	

	
// DIAGNOSTIC

	
	/** Diagnostic. Prints the time t, vector H and field pi.
	 */
	void printState()
	{
         cout << "\n\n Diagnostic: node at time step s = " << s
		      << "\nVector H: " << H
		      << "\nPrice pi = " << pi;
	}
	 
	
	
// LIBORS, SWAPRATES, ANNUITY (PBV)
	
	/** Libor \f$X_j=\delta_jL_j\f$ at this node.
	 */
	Real X(int j)
    {
		int t=get_t();             // node in (T_{t-1},T_t]
		if((t==n)&&(j==n-1)) return (H[j]-1.0);
		if((t<n)&&(j>=t)&&(j<n)) return (H[j]-H[j+1])/H[j+1];
		// else
		cout << "\n\nLmmNode:X(j): j = " << j << " not in ["<<t<<","<<n-1<<"]"
		     << "\nTime step s = " << s
		     << "Terminating.";
		exit(0);		
	}

	
	/** Forward price of the annuity (PBV)
	 *  \f$H_{p,q}=B_{p,q}/B_n=\sum\nolimits_{j=p}^{q-1}\delta_jH_{j+1}\f$
	 *  at this node.
	 */
	Real Hpq(int p, int q);
	
	
	/** Forward swaprate \f$S_{p,q}\f$ at this node, swap interval
	 *  \f$[T_p,T_q]\f$, where \f$t\leq p<q\f$. Here t is the discrete time
	 *  at which the node lives.
	 */	  
	Real swapRate(int p, int q){ return (H[p]-H[q])/Hpq(p,q); }
		

// ZERO COUPON BONDS
		
    /** The zero coupon bond \f$B_i\f$ expiring at time \f$T_i\f$ at this node. 
	 */
	Real B_i(int i)
	{ 
		int t=get_t();
		return H[i]/H[t]; 
	}
	

}; // end Node
	


/**********************************************************************************
 *
 *            NODE IN 2 FACTOR LMM LATTICE
 *
 *********************************************************************************/




/** Node in a 2 factor {@link LmmLattice2F}.
 *  Contains the vector \f$V(T_t)=(V_t(T_t),...,V_{n-1}(T_t))\f$
 *  as well as the variables \f$V1=V_{n-1}, V2=V_{n-2}-V_{n-1}\f$ which the lattice 
 *  evolves.
 */
struct LmmNode2F : public LmmNode {
	
	
	Real V1, V2;            // the variables which are evolved
	int i,                  // V1=V1(0)+ia_1
		j;                  // V2=V2(0)+ja_2,

	
public:
		
// ACCESSORS
	
	Real getV1(){ return V1; }
	Real getV2(){ return V2; }
	void setV1(Real v){ V1=v; }
	void setV2(Real v){ V2=v; }
	int get_i(){ return i; }
	int get_j(){ return j; }
	
	
// CONSTRUCTOR
	
	/** @param s number of time steps to reach this node from time zero.
	 *  @param nSteps number of time steps in each Libor acrual interval.
	 *  @param k state V_1=V_1(0)+ka
	 *  @param l state V_2=V_2(0)+la
	 */
	LmmNode2F(LiborFactorLoading* fl, int s, int nSteps, int k, int l) : LmmNode(fl,s,nSteps),
    i(k), j(l)
	{	}
				 
	
	 
	 
// COMPUTATION OF THE STATE	 

	
	
	/** Diagnostic. Checking the vector V computed at a node at time t.
	 */
	void checkState(int t, vector<Real>& V, vector<Real>& Z, Matrix<Real>& RQt);	
	

}; // end LmmNode2F
	


	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
