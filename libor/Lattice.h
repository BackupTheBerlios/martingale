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

#include "Node.h"
#include "Matrices.h"
#include <list>

MTGL_BEGIN_NAMESPACE(Martingale)


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
 * \f[E_t(h)=E_t[E_{t+1}(h)]\f]
 * we obtain at each node nd at time t
 * \f[\pi_t(nd)=\sum\nolimits_{e}p(e)*\pi_{t+1}(nd(e))\f]
 * where the sum extends over all edges e originating from the node nd and p(e) and nd(e)
 * denote the transition probability along the edge e and the node at time t+1 which the 
 * edge e connects to. Here h could be a conditional expectation \f$h=E_s(k)\f$ of some
 * random variable k at a later time s. Thus the case of iterated conditional expectations
 * is handled with equal ease and this is the main advantage of lattices over Monte Carlo
 * path simulation.
 *
 * <p>The Lattice consists of nodes of type Node. The class Node must 
 * have a member function getEdges() returning a std::list<Edge> of the edges in the node.
 * Since the type of node varies with the lattice we make it a template parameter.
 * The concrete Node subclasses which will derive from Node can then be defined as local
 * classes in the various lattices that use them.
 *
 * @param Node the type of nodes in this lattice.
 */
template<typename Node_t>
class Lattice {
	
protected:
	
	/** Number of time steps in the lattice. */
	int m;
	/** Number of nodes. */
	int nodes;         
				
	// nodeList[t] is the list of pointers to nodes at time t
	vector< std::list<Node_t*>* > nodeList; 

	
public:
	

// ACCESSORS

    int getTimeSteps() const { return m; }
	
	/** The list of nodes at time t.
	 */
	std::list<Node_t*>& getNodeList(int t) { return *(nodeList[t]); }
	
	/** The root node of the lattice.
	 */
	Node_t* getRoot(){ return nodeList[0]->front(); }
	
	
// CONSTRUCTOR
	
	/** @param steps number of time steps in the lattice.
	 */
	Lattice(int steps) : 
	m(steps), nodes(0), nodeList(m+1)
	{  
		// the lists of nodes at each time t
		for(int t=0;t<=m;t++)
		nodeList[t]=new std::list<Node_t*>(); // empty list
	}
	
	~Lattice()
    {
		// the lists of nodes at each time t
		for(int t=0;t<=m;t++){
			
		   std::list<Node_t*>& nodes=getNodeList(t);
		   std::list<Node_t*>::const_iterator theNode;
		   for(theNode=nodes.begin(); theNode!=nodes.end(); ++theNode) delete &(*theNode);
		   }
	}
			 
		   
		
	
	/** Goes through all the nodes in the lattice and checks if the transition 
	 *  probabilities are in [0,1] and sum to 1. Note however that this class does 
	 *  not allocate any nodes or edges. This is left to the derived classes.
	 */
	void selfTest()
    {	
		cout << endl << endl
		     << "Checking nodes, total number: " << nodes << endl;
		// loop over times t
		for(int t=0;t<m;t++){
			
			cout << "\nNodes at time t = " << t << ";";
			
			// probabilities at each node and their sum
			Real p,psum=0.0;   
			bool prob_failure=false;
						
			std::list<Node_t*>* nodes=nodeList[t];
			// iterate over list elements, iterator is pointer to pointer to node
			std::list<Node_t*>::const_iterator theNode;
			for(theNode=nodes->begin(); theNode!=nodes->end(); ++theNode){
				
				Node_t* node=*theNode;
				std::list<Edge>& edges=node->getEdges();
		        std::list<Edge>::const_iterator theEdge;           // pointer to edge
	        	for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge){
					
			    	psum=0.0;
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
				
			} // end for it
							
            cout << " OK";
			
		} // end for t
		
	} // end selfTest()

				

}; // end Lattice




	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 