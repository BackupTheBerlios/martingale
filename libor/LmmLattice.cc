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
using namespace Martingale;


		

/**********************************************************************************
 *
 *            GENERAL TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
	
    LmmLattice2F::
	LmmLattice2F(LiborFactorLoading* fl, int s) : LmmLattice<LmmNode2F>(fl,s),
	a(0), 
	k1(n-1), k2(n-1),
	K1(n-1), K2(n-1),
	D11(n-2), D22(n-2), D12(n-2), 
	Q11(n-2), Q22(n-2), Q12(n-2), 
	R(n-1),
	RQ(n-1)
    {  
		// the covariance integrals, note the decorrelations
		// V1 <--> V_{n-1}, V2 <--> V_{n-2}-V_{n-1}
		Real* Tc=factorLoading->getTenorStructure();           // Tc[t]=T_t
		for(int t=0;t<m;t++){
			
			Real d11, d12, d22;
			d11=factorLoading->integral_sgi_sgj_rhoij(n-1,n-1,Tc[t],Tc[t+1]);
			d22=factorLoading->integral_sgi_sgj_rhoij(n-2,n-2,Tc[t],Tc[t+1]);
			d12=factorLoading->integral_sgi_sgj_rhoij(n-2,n-1,Tc[t],Tc[t+1]);
			D11[t]=d11;
			D12[t]=d12-d11;
			D22[t]=d22+d11-2*d12;
			
		} // end for s

		setRank2CovariationMatrixRoots(m);
        setTickSize(m);
		buildLattice(m);
		testFactorization();
		
	} // end constructor
	
				
	
// Covariance matrix rank 2 factorization error
	

void LmmLattice2F::
testFactorization()
{
	cout << "\n\nRelative errors of the approximate rank 2 factorizations of all"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank 2,"
	     << "\n(trace norm): " << endl << endl;
	
	Real* Tc=factorLoading->getTenorStructure();
	for(int t=1;t<n-2;t++){
		
		Real T_t=Tc[t];
		UTRMatrix<Real> Ct(n-t,t);
		for(int i=t;i<n;i++)
		for(int j=i;j<n;j++) Ct(i,j)=factorLoading->integral_sgi_sgj_rhoij(i,j,0.0,T_t);
        
		Ct.testFactorization(2);
	}
}   // end factorAnalysis


	
	
	
// TEST
	
	void LmmLattice2F::
    test(int n)
    {
		Timer watch; watch.start();
		LiborFactorLoading* fl=JR_FactorLoading::sample(n);
		LmmLattice2F lattice(fl,n-3);
		lattice.selfTest();
		watch.stop();
		watch.report("Libor tree construction and self test");
	}
		
	
// methods used to build the lattice
	
	UTRMatrix<Real>& LmmLattice2F::
	decorrelatedCovariationMatrix(int t)
    {
	    Real* Tc=factorLoading->getTenorStructure();
		Real T_t=Tc[t];  
		UTRMatrix<Real>& C=*(new UTRMatrix<Real>(n-t,t));
		for(int i=t;i<n;i++){
			
			// the covariances Ci[j]=E[V_i(T_t)V_j(T_t)],
			// j=i,...,n-1
			vector<Real> Ci(n-i,i);
		    for(int j=i;j<n;j++) Ci[j]=factorLoading->integral_sgi_sgj_rhoij(i,j,0.0,T_t);
			
			// set covariances, decorrelated in the last two components
			for(int j=i;j<n;j++)
			if(j!=n-2) C(i,j)=Ci[j];
			else                                    // j=n-2, i!=n-2
			if(i!=n-2) C(i,j)=Ci[n-2]-Ci[n-1];
			else{                                   // j=n-2, i=n-2
				
				// C_{n-1,n-1}
				Real Cnn=factorLoading->integral_sgi_sgj_rhoij(n-1,n-1,0.0,T_t);
				C(i,j)=Ci[n-2]-2*Ci[n-1]+Cnn;
			} // end for j
		} // end for i
		
		return C;
	} // end decorrelatedCovariationMatrix
	
			
				

	void LmmLattice2F::
	setRank2CovariationMatrixRoots(int m)
    {
	    for(int t=1;t<=m;t++){
			
			// this is indexed as i,j=t,...,n-1
			UTRMatrix<Real>& Ct=decorrelatedCovariationMatrix(t);
			Matrix<Real>& Rt=Ct.rankReducedRoot(2); // row index i=t,...,n-1
			// save
			R[t]=&Rt;
			
			// RQ[t], 2 by 2 inverse of the last two rows of R[t]
			// rows of R[t] indexed from base t
			Matrix<Real>* RQt_ptr=new Matrix<Real>(2);
			Matrix<Real>& RQt=*RQt_ptr;
			
			// the last two rows of Rt
			Real A=Rt(n-2,0), B=Rt(n-2,1),
			     C=Rt(n-1,0), D=Rt(n-1,1);
			
			// inversion by hand, determinant
			Real det=A*D-B*C;
			// inverse
			RQt(0,0)=D/det; RQt(0,1)=-B/det;
			RQt(1,0)=-C/det; RQt(1,1)=A/det;
			// save
			RQ[t]=RQt_ptr;

			
		}
	} // end setRank2CovariationMatrixRoots
	
	

	void LmmLattice2F::
	setTickSize(int m)
    {
		// search for a
		a=0.05;
		bool go_on=true;
		
		while(go_on &&(a>0.0)){
			
			go_on=false;    // reset to true if condition violated
			a-=0.001;

			for(int t=0;t<m;t++){
				     
		    	k1[t]=1+(int)(sqrt(D11[t])/a);
		        k2[t]=1+(int)(sqrt(D22[t])/a);
				
			    Real A1=k1[t]*a,
				     A2=k2[t]*a,
				     p1,p2;
				
				Q11[t]=D11[t]/(A1*A1),
				Q22[t]=D22[t]/(A2*A2),
				Q12[t]=D12[t]/(A1*A2),
				p1=Q11[t]+Q22[t]+Q12[t]-1,
				p2=Q11[t]+Q22[t]-Q12[t]-1;
				
			    go_on=((p1<0)||(p2<0));
			    if(go_on) break;
			} // end for t
		} // end while
				
		
		cout << "\n k1[t]: " << k1
		     << "\n k2[t]: " << k2;
		
		// bounds K_j(t) for the number of states
		K1[0]=0;
		K2[0]=0;
		
		unsigned int nodes=(2*K1[0]+1)*(2*K2[0]+1);
		for(int t=0;t<m;t++){ 
			
			K1[t+1]=K1[t]+k1[t]; 
			K2[t+1]=K2[t]+k2[t]; 
			nodes+=(2*K1[t+1]+1)*(2*K2[t+1]+1); 
		}
		
		cout << "\n a= " << a;
		cout << "\nMaximum possible number of nodes = "  << nodes;
			
	} // end setTickSize
	
	
	
	 void LmmLattice2F::
	 setStateVariables(LmmNode2F* node)
	 {
         int t=node->getTime();
		 vector<Real> V(n-t,t);   // the volatility parts V_j of log(U_j)
			                          // j=t,...,n-1
		 V[n-1]=node->getV1();
		 V[n-2]=node->getV2();
			 
		 // matrices needed to derive V_j, j<n-2, from V_{n-2}, V_{n-1}
		 Matrix<Real>& Rt=getR(t);
		 Matrix<Real>& RQt=getRQ(t);
			 
		 // the "factors" Z_1, Z_2 for V_{n-2}, V_{n-1}
		 vector<Real> Z(2);   
		 Z[0]=V[n-2]; Z[1]=V[n-1]; Z*=RQt;
		 
		 // use these to compute the remaining V_j
		 // Rt=*R[t] has row index base t.
		 for(int j=t;j<n;j++) V[j]=Rt(j,0)*Z[0]+Rt(j,1)*Z[1]; 
		 node->checkState(t,V,Z,RQt); 
		 
		 // eliminate the differencing V_{n-2}-V_{n-1} in coordinate n-2
		 V[n-2]+=V[n-1];

		 // add the initial value log(U_j(0) and drift mu_j(t)=mu_j(0,T_t) 
		 // to obtain the log(U_j(T_t))
		 for(int j=t;j<n;j++) V[j]+=log(U0[j])+mu(j,t);
				 
		 // move from log(U_j) to U_j
		 for(int j=t;j<n;j++) V[j]=exp(V[j]);    // now V=U
				 
		 // write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
		 vector<Real>& H=node->getH();
		 H[n]=1.0; 
	     for(int j=n-1;j>=t;j--){ 
			 
			 H[j]=V[j]+H[j+1]; 
			 //if(H[j]>3.0) { printState(); exit(0); }
		}
			 
	} // setStateVariables


		

	void LmmLattice2F::
	buildLattice(int m)
    {
		cout << "\n\n\nBuilding lattice until time t = " << m << endl << endl;
		  
		std::list<LmmNode2F*>& nodes_0=*(nodeList[0]);               // list of nodes at time t
		
		// the initial node
		LmmNode2F* nextNode=new LmmNode2F(factorLoading,0,0,0);
		nextNode->setV1(0.0);
		nextNode->setV2(0.0);
		for(int j=0;j<=n;j++) nextNode->getH()[j]=H0[j];
			
	    // enter into node list at time t=0
		nodes_0.push_back(nextNode); 
		nodes++;
		
		// loop over times s=0,...,t-1
		for(int t=0;t<m;t++) {
	
			Real A1=k1[t]*a, A2=k2[t]*a;                             // tick sizes at time t
			std::list<LmmNode2F*>& listCurrentNodes=*(nodeList[t]);    // list of nodes at time t+1
			std::list<LmmNode2F*>& listNewNodes=*(nodeList[t+1]);      // list of nodes at time t+1

			// registry for possible nodes at time t+1, n(t+1,i,j)
			// i=-K_1(t+1),...,K_1(t+1), j=-K_2(t+1),...,K_2(t+1)
			Matrix<LmmNode2F*> newNodes(2*K1[t+1]+1,2*K2[t+1]+1,-K1[t+1],-K2[t+1]);
			
			// run through the list of nodes at time t and connect to the new nodes
			// iterator is pointer to pointer to node
			std::list<LmmNode2F*>::const_iterator theNode;
			for(theNode=listCurrentNodes.begin(); theNode!=listCurrentNodes.end(); ++theNode)			
			{
				LmmNode2F* currentNode=*theNode;
				// state of this node
				Real V1=currentNode->getV1(), V2=currentNode->getV2();
				int i=currentNode->get_i(), j=currentNode->get_j();
				// list of edges of this node
				std::list<Edge>& edgeList=currentNode->getEdges();
				
				// loop over transitions k,l=-1,0,1
				for(int k=-1;k<2;k++)
				for(int l=-1;l<2;l++) if(!(k==0&&l==0)) 
				{	
					// connect to node n(s+1,i+k,j+l) in registry
					nextNode=newNodes(i+k*k1[t],j+l*k2[t]);
					// if it does not exist yet
					if(!nextNode){ 
						
						  nextNode=new LmmNode2F(factorLoading,t+1,i+k*k1[t],j+l*k2[t]); 
						  newNodes(i+k*k1[t],j+l*k2[t])=nextNode;             // register the node
						  listNewNodes.push_back(nextNode);                   // append to list
						  nextNode->setV1(V1+k*A1);
						  nextNode->setV2(V2+l*A2);
						  setStateVariables(nextNode);
						  nodes++; 
					}
					
					// make an edge to the new node
					Edge* edge=new Edge();
					edge->node=nextNode;
					edge->probability=probability(t,k,l);
					edgeList.push_back(*edge);
					
			   } // end for k,l
			} // end for theNode

				
			cout << "\nTotal nodes = " << nodes;
		
		} // end for t
		
	} // end buildLattice
	
	

	Real LmmLattice2F::
	probability(int t, int i, int j)
    {
		if(j==0) return 0.5*(1-Q22[t]);
		if(i==0) return 0.5*(1-Q11[t]);
		if(i+j==0) return 0.25*(Q11[t]+Q22[t]-Q12[t]-1);
		else return 0.25*(Q11[t]+Q22[t]+Q12[t]-1);
	}
	
		


