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

#include "LmmLattice.h"
#include "LiborFactorLoading.h"
#include "Node.h"
#include "Array.h"
#include "Utils.h"
#include <cmath>


MTGL_BEGIN_NAMESPACE(Martingale)






template<typename LmmNode>
std::ostream& operator << (std::ostream& os, const LmmLattice<LmmNode>& ltt)
{
	return ltt.printSelf(os);
}


				

/**********************************************************************************
 *
 *            CONSTANT VOLATILITY TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/




void 
ConstVolLmmLattice2F::
buildLattice(int m)
{
	std::cout << "\n\nBuilding lattice: " << *this
	          << "\nTime steps: " << m << endl << endl;
	          		  
	std::list<LmmNode2F*>& nodes_0=*(nodeList[0]);               // list of nodes at time t
		
	// the initial node, i=j=0.
	LmmNode2F* nextNode=new LmmNode2F(factorLoading,0,nSteps,0,0);
	for(int j=0;j<=n;j++) nextNode->getH()[j]=H0[j];
			
	// enter into node list at time t=0
	nodes_0.push_back(nextNode); 
	nodes++;
		
	// for t=0,...,m-1 build nodes at time t+1 from nodes at time t
	for(int s=0;s<m;s++) {
	
		std::list<LmmNode2F*>& listCurrentNodes=*(nodeList[s]);    // list of nodes at time t
		std::list<LmmNode2F*>& listNewNodes=*(nodeList[s+1]);      // list of nodes at time t+1

		// registry for possible nodes at time t+1, n(t+1,i,j), 
		// -(t+1) <= i,j <= t+1
		Array2D<LmmNode2F*> newNodes(2*s+3,2*s+3,-s-1,-s-1);
			
		// run through the list of nodes at time t and connect to the new nodes
		// iterator is pointer to pointer to node
		std::list<LmmNode2F*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			LmmNode2F* currentNode=*theNode;
			// state of this node
			int i=currentNode->get_i(), j=currentNode->get_j();
			// list of edges of this node
			std::list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions k,l=-1,1 (up/down tick)
				for(int k=-1;k<2;k+=2)
				for(int l=-1;l<2;l+=2) 
				{	
					// connect to node n(s+1,i+k,j+l) in registry
					nextNode=newNodes(i+k,j+l);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new LmmNode2F(factorLoading,s+1,nSteps,i+k,j+l); 
						  newNodes(i+k,j+l)=nextNode;                         // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  setStateVariables(nextNode);
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.25;
					edgeList.push_back(*edge);
					
			   } // end for k,l
	    } // end for theNode

				
		cout << "\nTotal nodes = " << nodes;
		
	} // end for s
} // end buildLattice
	

	

void 
ConstVolLmmLattice2F::
test(int n) const
{
	Timer watch; watch.start();
	LiborFactorLoading*
	fl=LiborFactorLoading::sample(n,VolSurface::CONST,Correlations::CS);
	ConstVolLmmLattice2F lattice(fl,n-3);
	lattice.selfTest();
	watch.stop();
	watch.report("3 factor LMM lattice construction and self test");
}



/**********************************************************************************
 *
 *            CONSTANT VOLATILITY THREE FACTOR LMM LATTICE
 *
 *********************************************************************************/
	
				



void 
ConstVolLmmLattice3F::
buildLattice(int m)
{
	std::cout << "\n\nBuilding lattice: " << *this
	          << "\nTime steps: " << m << endl << endl;
		  
	std::list<LmmNode3F*>& nodes_0=*(nodeList[0]);               // list of nodes at time t
		
	// the initial node, i=j=k=0.
	LmmNode3F* nextNode=new LmmNode3F(factorLoading,0,nSteps,0,0,0);
	for(int j=0;j<=n;j++) nextNode->getH()[j]=H0[j];
			
	// enter into node list at time t=0
	nodes_0.push_back(nextNode); 
	nodes++;
		
	// for t=0,...,m-1 build nodes at time t+1 from nodes at time t
	for(int s=0;s<m;s++) {
	
		std::list<LmmNode3F*>& listCurrentNodes=*(nodeList[s]);    // list of nodes at time t
		std::list<LmmNode3F*>& listNewNodes=*(nodeList[s+1]);      // list of nodes at time t+1

		// registry for possible nodes at time t+1, n(t+1,i,j), 
		// -(t+1) <= i,j <= t+1
		Array3D<LmmNode3F*> newNodes(2*s+3,2*s+3,2*s+3,-s-1,-s-1,-s-1);
			
		// run through the list of nodes at time t and connect to the new nodes
		// iterator is pointer to pointer to node
		std::list<LmmNode3F*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			LmmNode3F* currentNode=*theNode;
			// state of this node
			int i=currentNode->get_i(), 
			    j=currentNode->get_j(),
			    k=currentNode->get_k();
			// list of edges of this node
			std::list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions i,j,k += +/-1 (up/down tick of Z1,Z2,Z3)
				for(int p=-1;p<2;p+=2)
				for(int q=-1;q<2;q+=2)
				for(int r=-1;r<2;r+=2)       // p,q,r = +/- 1
				{	
					// connect to node n(s+1,i+p,j+q,k+r) in registry
					nextNode=newNodes(i+p,j+q,k+r);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new LmmNode3F(factorLoading,s+1,nSteps,i+p,j+q,k+r); 
						  newNodes(i+p,j+q,k+r)=nextNode;                         // register the node
						  listNewNodes.push_back(nextNode);                       // append to list
						  setStateVariables(nextNode);
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.125;            // 1/8
					edgeList.push_back(*edge);
					
			   } // end for p,q,r
		} // end for theNode

		cout << "\nTotal nodes = " << nodes;
	} // end for t
} // end buildLattice
	
	



void 
ConstVolLmmLattice3F::
test(int n) const
{
	Timer watch; watch.start();
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(n,VolSurface::CONST,Correlations::CS);
	ConstVolLmmLattice3F lattice(fl,n-3);
	lattice.selfTest();
	watch.stop();
	watch.report("3 factor LMM lattice construction and self test");
}




MTGL_END_NAMESPACE(Martingale)




