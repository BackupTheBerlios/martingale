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

#include "BasketLattice.h"
#include "Lattice.h"
#include "Matrices.h"
#include "Utils.h"
#include <string>
#include <iostream>
#include <math.h>
using namespace Martingale;




/**********************************************************************************
 *
 *            TWO FACTOR LATTICE
 *
 *********************************************************************************/
	
BasketLattice2F::
BasketLattice2F
(int _n, int _T, Real _dt, RealVector& _S0, RealVector& _sg, UTRRealMatrix& _rho) :
BasketLattice<BasketNode2F>(_n,_T,_dt,_S0,_sg,_rho),
R(_rho.rankReducedRoot(2))
{  	
	// scale the rows of R back to norm one 
	// this preserves rho_{jj}=1 and volatilities
	for(int i=0;i<n-1;i++){
			
		Real f=0.0;      // norm of row_i(R)
		for(int j=0;j<2;j++) f+=R(i,j)*R(i,j);
		f=sqrt(f);
		R.scaleRow(i,1.0/f);
	}
		
	buildLattice();       // build lattice with T time steps
	testFactorization();
		
} // end constructor
	
				
BasketLattice2F* 
BasketLattice2F::
sample(int n, int T)
{
	// initial asset prices
	RealVector S0(n);
	for(int j=0;j<n;j++) S0[j]=100.0;
			
	// volatilities
	RealVector sg(n);
	for(int j=0;j<n;j++) sg[j]=0.3;
			
	// correlation of returns
	UTRRealMatrix rho(n);
    for(int j=0;j<n;j++)
	for(int k=j;k<n;k++) rho(j,k)=exp(0.2*(j-k));
			
	// time step
	Real dt=0.1;
		
	return new BasketLattice2F(n,T,dt,S0,sg,rho);
}
	
	
// Covariance matrix rank 2 factorization error
	

void 
BasketLattice2F::
testFactorization() const
{
	cout << "\n\nRelative errors of the approximate rank 2 factorization"
	     << "\nrho=RR' (rank(R)=2) of the correlation matrix rho"
	     << "\n(trace norm): " << endl << endl;
  
	rho.testFactorization(2);

}   // end testFactorization



void BasketLattice2F::
buildLattice()
{
	cout << "\n\n\nBuilding lattice until time step s = " << T << endl << endl;
	
	// list of nodes at time s=0
	std::list<BasketNode2F*>& nodes_0=getNodeList(0);              
		
	// the initial node, i=j=0.
	BasketNode2F* nextNode=new BasketNode2F(n,0,0,0);
	for(int j=0;j<n;j++) nextNode->getS()[j]=S0[j];
			
	// enter into node list at time s=0
	nodes_0.push_back(nextNode); 
	nodes++;
		
	// for s=0,...,m-1 build nodes at time s+1 from nodes at time s
	for(int s=0;s<T;s++) {
	
		std::list<BasketNode2F*>& listCurrentNodes=getNodeList(s);       // list of nodes at time s
		std::list<BasketNode2F*>& listNewNodes=getNodeList(s+1);         // list of nodes at time s+1

		// registry for possible nodes at time s+1, n(s+1,i,j), 
		// -(s+1) <= i,j <= s+1
		Array2D<BasketNode2F*> newNodes(2*s+3,2*s+3,-s-1,-s-1);
			
		// run through the list of nodes at time s and connect to the new nodes
		// iterator is pointer to pointer to node
		std::list<BasketNode2F*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			BasketNode2F* currentNode=*theNode;
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
						
						  nextNode=new BasketNode2F(n,s+1,i+k,j+l); 
						  newNodes(i+k,j+l)=nextNode;                         // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  setAssetPrices(nextNode);
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.25;
					edgeList.push_back(*edge);
					
			   } // end for k,l
		} // end for theNode

        cout << "\nTime step s = " << s << ", total nodes = " << nodes;
		
	} // end for s
} // end buildLattice
	
	
		
void 
BasketLattice2F::
setAssetPrices(BasketNode2F* node)
{
     int s=node->getTimeStep();   // time step at which node lives
			         
     Real Z1=node->get_i()*a, 
          Z2=node->get_j()*a;     // Z2=ja
		 
	 // compute asset prices
	 RealArray1D& S=node->getS();
		 
	 // volatility parts V_j of Y_j=log(S_j)
	 for(int j=0;j<n;j++) S[j]=sg[j]*(R(j,0)*Z1+R(j,1)*Z2); 
		 
     // V_j -> Y_j: add on initial value and drift
	 for(int j=0;j<n;j++) S[j]+=log(S0[j])+s*driftunit[j];  
				 
	 // move from Y_j=log(S_j) to S_j
	 for(int j=0;j<n;j++) S[j]=exp(S[j]);    
			 
} // setAssetPrices

	
	
	
// TEST
	
void 
BasketLattice2F::
test(int n, int T)
{
	Timer watch; watch.start();
	BasketLattice2F* lattice=sample(n,T);
	lattice->selfTest();
	watch.stop();
	watch.report("Libor tree construction and self test");
}



/**********************************************************************************
 *
 *            THREE FACTOR LMM LATTICE
 *
 *********************************************************************************/
	
	
BasketLattice3F::
BasketLattice3F
(int _n, int _T, Real _dt, RealVector& _S0, RealVector& _sg, UTRRealMatrix& _rho):
BasketLattice<BasketNode3F>(_n,_T,_dt,_S0,_sg,_rho),
R(_rho.rankReducedRoot(3))
{  	
	// scale the rows of R back to norm one 
	// this preserves rho_{jj}=1 and volatilities
	for(int i=0;i<n-1;i++){
			
		Real f=0.0;      // norm of row_i(R)
		for(int j=0;j<3;j++) f+=R(i,j)*R(i,j);
		f=sqrt(f);
		R.scaleRow(i,1.0/f);
	}
		
	buildLattice();       // build lattice with T time steps
	testFactorization();
		
} // end constructor
	
				
BasketLattice3F* 
BasketLattice3F::
sample(int n, int T)
{
	// initial asset prices
	RealVector S0(n);
	for(int j=0;j<n;j++) S0[j]=100.0;
			
	// volatilities
	RealVector sg(n);
	for(int j=0;j<n;j++) sg[j]=0.3;
			
	// correlation of returns
	UTRRealMatrix rho(n);
    for(int j=0;j<n;j++)
	for(int k=j;k<n;k++) rho(j,k)=exp(0.2*(j-k));
			
	// time step
	Real dt=0.1;
		
	return new BasketLattice3F(n,T,dt,S0,sg,rho);
}
	
	
// Covariance matrix rank 2 factorization error
	
void 
BasketLattice3F::
testFactorization() const
{
	cout << "\n\nRelative errors of the approximate rank 3 factorization"
	     << "\nrho=RR' (rank(R)=3) of the correlation matrix rho"
	     << "\n(trace norm): " << endl << endl;
  
	rho.testFactorization(3);

}   // end testFactorization



void BasketLattice3F::
buildLattice()
{
	cout << "\n\n\nBuilding lattice until time step s = " << T << endl << endl;
	
	// list of nodes at time s=0
	std::list<BasketNode3F*>& nodes_0=getNodeList(0);            
		
	// the initial node, i=j=0.
	BasketNode3F* nextNode=new BasketNode3F(n,0,0,0,0);
	for(int j=0;j<n;j++) nextNode->getS()[j]=S0[j];
			
	// enter into node list at time s=0
	nodes_0.push_back(nextNode); 
	nodes++;
		
	// for s=0,...,m-1 build nodes at time s+1 from nodes at time s
	for(int s=0;s<T;s++) {
	
		std::list<BasketNode3F*>& listCurrentNodes=getNodeList(s);       // list of nodes at time s
		std::list<BasketNode3F*>& listNewNodes=getNodeList(s+1);         // list of nodes at time s+1

		// registry for possible nodes at time s+1, n(s+1,i,j,k), 
		// -(s+1) <= i,j,k <= s+1
		Array3D<BasketNode3F*> newNodes(2*s+3,2*s+3,2*s+3,-s-1,-s-1,-s-1);
			
		// run through the list of nodes at time t and connect to the new nodes
		// iterator is pointer to pointer to node
		std::list<BasketNode3F*>::const_iterator theNode;
		for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
		{
			BasketNode3F* currentNode=*theNode;
			// state of this node
			int i=currentNode->get_i(), 
			    j=currentNode->get_j(),
			    k=currentNode->get_k();
			// list of edges of this node
			std::list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions k,l=-1,1 (up/down tick)
				for(int p=-1;p<2;p+=2)
				for(int q=-1;q<2;q+=2)
				for(int r=-1;r<2;r+=2) 
				{	
					// connect to node n(s+1,i+k,j+l) in registry
					nextNode=newNodes(i+p,j+q,k+r);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new BasketNode3F(n,s+1,i+p,j+q,k+r); 
						  newNodes(i+p,j+q,k+r)=nextNode;                     // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  setAssetPrices(nextNode);
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=0.125;
					edgeList.push_back(*edge);
					
			   } // end for p,q,r
		} // end for theNode

		cout << "\nTime step s = " << s << ", total nodes = " << nodes;
		
	} // end for s
} // end buildLattice
	
	

void 
BasketLattice3F::
setAssetPrices(BasketNode3F* node)
{
     int s=node->getTimeStep();   // time step at which node lives
			         
     Real Z1=node->get_i()*a, 
	      Z2=node->get_j()*a,     // Z2=ja
		  Z3=node->get_k()*a;
		 
     // compute asset prices
	 RealArray1D& S=node->getS();
		 
	 // volatility parts V_j of Y_j=log(S_j)
	 for(int j=0;j<n;j++) S[j]=sg[j]*(R(j,0)*Z1+R(j,1)*Z2+R(j,2)*Z3); 
		 
     // V_j -> Y_j: add on initial value and drift
	 for(int j=0;j<n;j++) S[j]+=log(S0[j])+s*driftunit[j];  
				 
	 // move from Y_j=log(S_j) to S_j
	 for(int j=0;j<n;j++) S[j]=exp(S[j]);    
			 
} // setAssetPrices

	
	
	
// TEST
	
void 
BasketLattice3F::
test(int n, int T)
{
	Timer watch; watch.start();
	BasketLattice3F* lattice=sample(n,T);
	lattice->selfTest();
	watch.stop();
	watch.report("Libor tree construction and self test");
}

	
	
	
	
	
	
	
	
	
	
	
	


