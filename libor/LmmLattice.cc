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
#include "Matrices.h"
#include "Utils.h"
#include <cmath>
using namespace Martingale;





				

/**********************************************************************************
 *
 *            CONSTANT VOLATILITY TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/


ConstVolLmmLattice2F::
ConstVolLmmLattice2F
(LiborFactorLoading* fl, int t, int steps=1) : 
LmmLattice<LmmNode2F>(fl,t*steps,steps),
nSteps(steps),
delta(fl->getDeltas()[0]),
dt(delta/steps),
a(sqrt(dt)), 
sg(n-1,1),
R(fl->getRho().rankReducedRoot(2))
{  
	// check if Libor accrual periods are constant
	const RealArray1D& deltas=fl->getDeltas();
	for(int j=0;j<n;j++) if(deltas[j]!=delta) {
			
	   std::cout << "\n\nConstVolLmmLattice2F(): Libor accrual periods not constant."
	             << "\nTerminating.";
	   exit(1);
	}
		
	// check if volatilities are constant
	if(fl->getVolSurfaceType()!=VolSurface::CONST) {
			
	   std::cout << "\n\nConstVolLmmLattice2F(): volatilities not constant."
	             << "\nTerminating.";
	   exit(1);
    }
		
	// set constant vols
	for(int j=1;j<n;j++) sg[j]=factorLoading->sigma(j,0.0);
			
	// scale the rows of R back to norm one 
	// this preserves rho_{jj}=1 and volatilities
	for(int i=1;i<n;i++){
			
		Real f=0.0;      // norm of row_i(R)
		for(int j=0;j<2;j++) f+=R(i,j)*R(i,j);
		f=sqrt(f);
		R.scaleRow(i,1.0/f);
	}
		
	buildLattice(m);       // Lattice::m = number of time steps t*nSteps 
	testFactorization();
		
} // end constructor
	
				

void 
ConstVolLmmLattice2F::
testFactorization() const
{
	std::cout << "\n\nRelative errors of the approximate rank 2 factorization"
	          << "\nrho=RR' (rank(R)=2) of the correlation matrix rho"
	          << "\n(trace norm): " << endl << endl;
  
	factorLoading->getRho().testFactorization(2);

}   // end factorAnalysis



void 
ConstVolLmmLattice2F::
buildLattice(int m)
{
	printSelf(std::cout);
	std::cout << "\nTime steps: " << m << endl << endl;
	          		  
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
setStateVariables(LmmNode2F* node)
{
    int t=node->get_t(),         // node lives in (T_{t-1},T_t]
	    s=node->getTimeStep();   // time step at which node lives
	RealVector V(n-t,t);       // the volatility parts V_j of log(U_j)
			                          // j=t,...,n-1
    Real Z1=node->get_i()*a, 
         Z2=node->get_j()*a;     // Z2=ja
		 
	for(int j=t;j<n;j++) V[j]=sg[j]*(R(j,0)*Z1+R(j,1)*Z2); 
		 
	// let tau(s) denote he continuous time reached with time step s.
	// add the initial value log(U_j(0) and drift mu_j(s)=mu_j(0,tau(s)) 
	// to obtain the log(U_j(tau(s)))
	for(int j=t;j<n;j++) V[j]+=log(U0[j])+mu(s,j);
				 
	// move from log(U_j) to U_j
	for(int j=t;j<n;j++) V[j]=std::exp(V[j]);    // now V=U
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
	RealArray1D& H=node->getH();
	H[n]=1.0; 
    for(int j=n-1;j>=t;j--){ 
			 
		 H[j]=V[j]+H[j+1]; 
		 //if(H[j]>3.0) { printState(); exit(0); }
	}
} // setStateVariables

	

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
	
ConstVolLmmLattice3F::
ConstVolLmmLattice3F
(LiborFactorLoading* fl, int t, int steps=1) : 
LmmLattice<LmmNode3F>(fl,t*steps,steps),
nSteps(steps),
delta(fl->getDeltas()[0]),
dt(delta/steps),
a(sqrt(dt)), 
sg(n-1,1),
R(fl->getRho().rankReducedRoot(3))
{  
	// check if Libor accrual periods are constant
	const RealArray1D& deltas=fl->getDeltas();
	for(int j=0;j<n;j++) if(deltas[j]!=delta) {
			
	   std::cout << "\n\nConstVolLmmLattice2F(): Libor accrual periods not constant."
	             << "\nTerminating.";
	   exit(0);
    }
		
	// check if volatilities are constant
	if(fl->getVolSurfaceType()!=VolSurface::CONST) {
			
	   std::cout << "\n\nConstVolLmmLattice2F(): volatilities not constant."
	             << "\nTerminating.";
	   exit(1);
    }

		
	// set constant vols
	for(int j=1;j<n;j++) sg[j]=factorLoading->sigma(j,0.0);
			
	// scale the rows of R back to norm one 
	// this preserves rho_{jj}=1 and volatilities
	for(int i=1;i<n;i++){
			
		Real f=0.0;      // norm of row_i(R)
		for(int j=0;j<3;j++) f+=R(i,j)*R(i,j);
		f=sqrt(f);
		R.scaleRow(i,1.0/f);
	}
		
	buildLattice(m);       // Lattice::m = number of time steps t*nSteps 
	testFactorization();
		
} // end constructor
	
				
	
// Covariance matrix rank 2 factorization error
	

void 
ConstVolLmmLattice3F::
testFactorization() const
{
	cout << "\n\nRelative errors of the approximate rank 3 factorization"
	     << "\nrho=RR' (rank(R)=3) of the correlation matrix rho"
	     << "\n(trace norm): " << endl << endl;
  
	factorLoading->getRho().testFactorization(3);

}   // end factorAnalysis



void 
ConstVolLmmLattice3F::
buildLattice(int m)
{
	printSelf(std::cout);
	std::cout << "\nTime steps: " << m << endl << endl;
		  
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
setStateVariables(LmmNode3F* node)
{
    int t=node->get_t(),         // node lives in (T_{t-1},T_t]
	     s=node->getTimeStep();   // time step at which node lives
	 RealVector V(n-t,t);       // the volatility parts V_j of log(U_j)
			                          // j=t,...,n-1
     Real Z1=node->get_i()*a,    
	      Z2=node->get_j()*a,
	      Z3=node->get_k()*a;      // Z3=ka
		 
	 for(int j=t;j<n;j++) V[j]=sg[j]*(R(j,0)*Z1+R(j,1)*Z2+R(j,2)*Z3); 
		 
	 // let tau(s) denote he continuous time reached with time step s.
	 // add the initial value log(U_j(0) and drift mu_j(s)=mu_j(0,tau(s)) 
	 // to obtain the log(U_j(tau(s)))
	 for(int j=t;j<n;j++) V[j]+=log(U0[j])+mu(s,j);
				 
	 // move from log(U_j) to U_j
	 for(int j=t;j<n;j++) V[j]=std::exp(V[j]);    // now V=U
				 
	 // write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
	 RealArray1D& H=node->getH();
	 H[n]=1.0; 
     for(int j=n-1;j>=t;j--){ 
			 
		 H[j]=V[j]+H[j+1]; 
		 //if(H[j]>3.0) { printState(); exit(0); }
	}
} // setStateVariables



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


