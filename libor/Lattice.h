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

#ifndef martingale_lattice_h    
#define martingale_lattice_h

#include "TypedefsMacros.h"
#include "Node.h"
#include <list>
#include <cstdlib> // exit()
#include <cmath>

MTGL_BEGIN_NAMESPACE(Martingale)

using std::list;


/*! \file Lattice.h
 *  General lattice as a vector of pointers to the list of nodes at
 *  each time t.
 */



/**********************************************************************************
 *
 *            GENERAL LATTICE
 *
 *********************************************************************************/


/** <p> A lattice is a compact representation of a large number of paths of a stochastic
 *  process with the property that the path set is invariant under path splitting at each
 *  time step. This makes the computation of a conditional expectation at time t
 *  \f$\pi_t=E_t(h)\f$
 *  at each node at time t a very simple procedure: applying the Double Expectation Theorem
 *  \f[E_t(h)=E_t[E_{t+1}(h)]\f]
 *  and assuming that the values of the random variable h are known by time t=s we can 
 *  compute \f$\pi_t\f$ by backward recursion starting from time t=s as follows:
 *  \f[\pi_s=h,\quad\; \pi_t=E_t[\pi_{t+1}],\quad t<s.\f]
 *  At each node N living at time t the conditional expectation \f$\pi_t=E_t[\pi_{t+1}]\f$
 *  is simply the probability weighted average
 * \f[\pi_t(N)=\sum\nolimits_{e}p(e)*\pi_{t+1}(N(e))\f]
 * where the sum extends over all edges e originating from the node N and p(e) and N(e)
 * denote the transition probability along the edge e and the node at time t+1 which the 
 * edge e connects to. At time t the values of \f$\pi_{t+1}\f$ have already been computed
 * at all nodes living at time t+1 and in particular at all the nodes N(e).
 *
 * <p>Here h could be a conditional expectation \f$h=E_u(k)\f$ of some
 * random variable k at a later time u. Thus the case of iterated conditional expectations
 * is handled with equal ease and this is the main advantage of lattices over Monte Carlo
 * path simulation.
 *
 * <p>The Lattice consists of nodes of type Node. The class Node must 
 * have a member function getEdges() returning a std::list<Edge> of the edges in the node.
 * Since the type of node varies with the lattice we make it a template parameter.

 *
 * @param Node the type of nodes in this lattice.
 */
template<typename Node>
class Lattice {
	
protected:
	
	/** Number of time steps in the lattice. */
	int T;   
				
	// nodeList[t] is the list of pointers to nodes at time t
	vector< list<Node*>* > nodeList; 

	
public:
	

// ACCESSORS

    int getTimeSteps() const { return T; }
	
	/** The list of nodes at time t.
	 */
	list<Node*>& getNodeList(int t) { return *(nodeList[t]); }
	
	/** The root node of the lattice.
	 */
	Node* getRoot(){ return nodeList[0]->front(); }
	
	
// CONSTRUCTOR
	
	/** @param s number of time steps in the lattice.
	 */
	Lattice(int s) : 
	T(s), nodeList(T+1)
	{  
		// the lists of nodes at each time t
		for(int t=0;t<=T;t++)
		nodeList[t]=new list<Node*>(); // empty list
	}

	
	~Lattice()
    {
		// the lists of nodes at each time t
		for(int t=0;t<=T;t++) getNodeList(t).~list<Node*>(); 
	}
		   		
	
	/** Goes through all the nodes in the lattice and checks if the transition 
	 *  probabilities are in [0,1] and sum to 1. Note however that this class does 
	 *  not allocate any nodes or edges. This is left to the derived classes.
	 */
	void selfTest()
    {	
		cout << endl << endl << "Checking nodes until time T-1 = " << T-1 << endl;
		// loop over times t
		for(int t=0;t<T;t++){
			
			std::cout << "\nNodes at time t = " << t << ";";
			
			// probabilities at each node and their sum
			Real p,psum=0.0;   
			bool prob_failure=false;
						
			list<Node*>* nodes=nodeList[t];

			// iterate over list elements, iterator is pointer to pointer to node
			list<Node*>::const_iterator theNode;
			for(theNode=nodes->begin(); theNode!=nodes->end(); ++theNode){
		
				Node* node=*theNode;
				list<Edge>& edges=node->getEdges();
		        list<Edge>::const_iterator theEdge;           // pointer to edge
				psum=0.0;
	        	for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge){
						
					p=theEdge->probability;
					if((p<0)||(p>1)) prob_failure=true;
					psum+=p;
				}
		        
				if(abs(psum-1.0)>0.0000001) prob_failure=true;
				if(prob_failure){
									
		         	cout << " *** t = " << t 
			             << "; transition probabilities defective, sum = " << psum;
					node->printTransitionProbabilities();
				    exit(0);
				 } // endif
				
			} // end for theNode
							
            cout << " OK";
			
		} // end for t
        	
	} // end selfTest()

				

}; // end Lattice



/** <a name="lattice-builder"><b>Lattice builders.</b></a>
 *  Standalone functions to build lattices with m time steps for which the state 
 *  variable is a standard Brownian motion Z in r dimensions (the number of factors). 
 *  Assumes that the nodes store the state Z_j=k[j]*a, a=sqrt(dt), j=0,1,...,r-1, 
 *  in the integer vector k and have construtors of the form
 *  <center><code>
 *  Node(int s, int* k, LatticeData* latticeData),
 *  </code></center>
 *  where s is the number of time steps to reach the node from time zero,
 *  k is the integer vector storing the state of Z and latticeData is  data object
 *  handed to the nodes from the lattice.
 *  Here dt is the size of the time step and a=sqrt(dt) the tick size of a standard 
 *  Brownian motion over a time step of size dt. The type Node must define a function 
 *  <code>int* getIntegerTicks()</code> returning 
 *  the state vector k.
 */
namespace LatticeBuilder {

/** Builds a <a href=#lattice-builder>two factor lattice</a>. 
 * @param T number of time steps in the lattice.
 * @param nodeList nodeList[t] is a pointer to the list of nodes at time t.
 * @param latticeData data object handed to nodes from lattice.
 */
template<class Node, class LatticeData>
void buildTwoFactorLattice
(int T, vector< std::list<Node*>* >& nodeList, LatticeData* latticeData)
{
	int nodes=0;                                            // counter
	list<Node*>& nodes_0=*(nodeList[0]);               // list of nodes at time t

	// array of integer tick multiples Z_j=k[j]*a, j=0,1, at each node
	IntArray1D k(2);

	// the initial node n(s,k0,k1), k0=k1=0.
	Node* nextNode=new Node(0,k,latticeData);
			
	// enter into node list at time s=0
	nodes_0.push_back(nextNode); 
	nodes++;

	// for s=0,...,T-1 build nodes at time s+1 from nodes at time s
	for(int s=0;s<T;s++) {

		list<Node*>& listCurrentNodes=*(nodeList[s]);    // list of nodes at time t
		list<Node*>& listNewNodes=*(nodeList[s+1]);      // list of nodes at time t+1

		// registry for possible nodes at time s+1, n(s+1,i,j), 
		// -(s+1) <= i=k[0],j=k[1] <= s+1
		Array2D<Node*> newNodes(2*s+3,2*s+3,-s-1,-s-1);
			
		// run through the list of nodes at time t and connect to the new nodes
		// iterator is pointer to pointer to node
		list<Node*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			Node* currentNode=*theNode;
			// state of this node
			int* l=currentNode->getIntegerTicks();
			// list of edges of this node
			list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions p,q=-1,+1 (up/down tick)
				for(int p=-1;p<2;p+=2)
				for(int q=-1;q<2;q+=2) 
				{	
					// integer state we transition to
					k[0]=l[0]+p; k[1]=l[1]+q;
					// connect to this node in the registry
					nextNode=newNodes(k[0],k[1]);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new Node(s,k,latticeData); 
						  newNodes(k[0],k[1])=nextNode;                       // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.25;
					edgeList.push_back(*edge);
					
			   } // end for p,q
	    } // end for theNode

				
		cout << "\nTime step = " << s+1 << ";  total nodes = " << nodes;
		
	} // end for s
} // end buildLattice

	

/** Builds a <a href=#lattice-builder>three factor lattice</a>. 
 * @param T number of time steps in the lattice.
 * @param nodeList nodeList[t] is a pointer to the list of nodes at time t.
 * @param latticeData data object handed to nodes from lattice.
 */
template<class Node, class LatticeData>
void buildThreeFactorLattice
(int T, vector< std::list<Node*>* >& nodeList, LatticeData* latticeData)
{

	int nodes=0;                                            // counter
	list<Node*>& nodes_0=*(nodeList[0]);               // list of nodes at time t

	// array of integer tick multiples Z_j=k[j]*a, j=0,1, at each node
	IntArray1D k(3);
	// the initial node n(s,k0,k1), k0=k1=0.
	Node* nextNode=new Node(0,k,latticeData);
			
	// enter into node list at time s=0
	nodes_0.push_back(nextNode); 
	nodes++;
	
	// for s=0,...,T-1 build nodes at time s+1 from nodes at time s
	for(int s=0;s<T;s++) {
	
		list<Node*>& listCurrentNodes=*(nodeList[s]);    // list of nodes at time t
		list<Node*>& listNewNodes=*(nodeList[s+1]);      // list of nodes at time t+1

		// registry for possible nodes at time s+1, n(s+1,i,j), 
		// -(s+1) <= i=k[0],j=k[1] <= s+1
		Array3D<Node*> newNodes(2*s+3,2*s+3,2*s+3,-s-1,-s-1,-s-1);
			
		// run through the list of nodes at time t and connect to the new nodes
		// iterator is pointer to pointer to node
		list<Node*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			Node* currentNode=*theNode;
			// state of this node
			int* l=currentNode->getIntegerTicks();
			// list of edges of this node
			list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions p,q=-1,+1 (up/down tick)
				for(int p=-1;p<2;p+=2)
				for(int q=-1;q<2;q+=2) 
				for(int r=-1;r<2;r+=2) 
				{	
					// integer state we transition to
					k[0]=l[0]+p; k[1]=l[1]+q; k[2]=l[2]+r;
					// connect to this node in the registry
					nextNode=newNodes(k[0],k[1],k[2]);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new Node(s,k,latticeData); 
						  newNodes(k[0],k[1],k[2])=nextNode;                  // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.125;
					edgeList.push_back(*edge);
					
			   } // end for p,q
	    } // end for theNode

				
		cout << "\nTime step = " << s+1 << ";  total nodes = " << nodes;
		
	} // end for s
} // end buildLattice


}; // end LatticeBuilder



	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
