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

#ifndef martingale_libortree2f_h    
#define martingale_libortree2f_h

#include "LiborFactorLoading.h"
#include "Utils.h"
#include <math.h>
#include <list>

MTGL_BEGIN_NAMESPACE(Martingale)



/** <p>Two factor lattice of the driftless Libor market model.
 *  See book 6.8, 8.1.1 for the details and notation. The variables
 *  evolved in the tree are \f$V1=V_{n-1},\ V2=V_{n-2}-V_{n-1}\f$ and 
 *  \f$t=0,1,\dots,n-2\f$ indexes the time steps \f$T_t\to T_{t+1}\f$.
 *  Current implementation constructs the lattice only up to discrete
 *  time \f$t=n-2\f$.
 *
 * <p>The lattice evolves the state variables \f$U_j\f$ of the driftless 
 * Libor market model but the nodes store the accrual factors \f$H_j\f$
 * since these allow us to compute Libors, swaprates, etc with the minimum effort.
 * See book, 6.8. 
 * 
 * <p>The variables \f$V_j\f$ are the unbounded variation parts of the
 * logarithms \f$log(U_j)\f$. Only these need to be modelled after the standard
 * simplifications See book 3.11, 8.1.1.
 */
class LiborTree2F {
	
	int n,             // dimension (number of accrual intervals)
	    m,             // lattice allocated until time t=m (T_m)
	    nodes;         // number of nodes
	
	vector<Real> logU0;     // logU0[j]=log(U_j(0))
// <------------- TO DO ------------> 
// work with 2 independent tick sizes	
	Real a;            // minimum tick size a_1=a_2=a
	
	vector<int> k1, k2;      // k_j(t)    
	vector<int> K1, K2;      // K_j(t)
	
	vector<Real> D11, D22, D12;  // D_{ij(t)
	vector<Real> Q11, Q22, Q12;  // Q_{ij(t)
	
	// forward declarations
	struct Edge;
	struct Node;
	
	// nodeList[t] is the list of pointers to nodes at time t
	// need this up to time t=n-2 (lattice ends t=n-2).
	vector< std::list<Node*>* > nodeList; 

	// factor loading of the underlying driftless LMM
	LiborFactorLoading* factorLoading; 
	
	// deterministic drifts mu(j,t)=-0.5*integral_0^{T_t}sigma_j^2(s)ds, 
	// 1<=j<=t<=n-1, index base 1, natural indexation
	LTRMatrix<Real> mu;
	
	// the sequence of rank 2 roots R(t) of the covariation matrices C(t) 
	// of the vectors V(t)=(V_t(T_t),...,V_{n-1}(T_t)).
	vector< Matrix<Real>* > R;
	
    // RQ[t] is the 2 by 2 inverse of the last two rows of R[t].
	vector< Matrix<Real>* > RQ;
	
public:

// ACCESSORS

    /** Dimension of underlying LMM (number of accrual periods).
	 */
    int getDimension(){ return n; }
	
	/** The matrix of deterministic drifts
	 *  \f[\mu(j,t)=\mu_j(0,T_t)=-{1\over2}\int_0^{T_t}\sigma_j^2(s)ds,\q t\leq j,\f]
	 *  of the logarithms \f$Y_j=log(U_j)\f$.
	 */
    LTRMatrix<Real>& getDrifts(){ return mu; }
	
	/** The rank 2 approximate root R(t) of the covariation matrix C(t) 
	 *  of the vector V(t)=(V_t(T_t),...,V_{n-1}(T_t)). See book 8.1.1.
	 */
    Matrix<Real>& getR(int t){ return *R[t]; }
	
    /** The 2 by 2 inverse of the last two rows of
	 *  {@link #getRank2CovariationMatrixRoot(int t)}.
	 */
    Matrix<Real>& getRQ(int t){ return *RQ[t]; }
	
    /** The vector logU0[j]=log(U_j(0)), initial values.
	 */
    vector<Real>& getLogU0(){ return logU0; }
	
	
// CONSTRUCTOR
	
	/** The lattice can only be built up to time \f$t=n-2\f$ at most.
	 *  The last step must be handled differently since only one state variable
	 *  remains to make this step. Will fix this later.
	 *
	 *  @param fl factor loading of the underlying LMM.
	 *  @param s lattice is built up to discrete time \f$t=s\leq n-2\f$ inclusively.
	 */
	LiborTree2F(LiborFactorLoading* fl, int s) :
	n(fl->getDimension()), 
	m(s),
	nodes(0),
	logU0(n),
	a(0), 
	k1(n-1), k2(n-1),
	K1(n-1), K2(n-1),
	D11(n-2), D22(n-2), D12(n-2), 
	Q11(n-2), Q22(n-2), Q12(n-2), 
	nodeList(n-1),
	factorLoading(fl),
	mu(n-1,1),
	R(n-1),
	RQ(n-1)
    {  
		// set log(U_j(0)), j=0,1,...,n-1
        Real* x=factorLoading->getInitialXLibors();     // x[j]=X_j(0)
        for(int j=0;j<n;j++){ 
			
			// U_j(0)=X_j(0)(1+X_{j+1}(0))...(1+X_{n-1}(0))
			Real Uj0=x[j]; for(int k=j+1;k<n;k++) Uj0*=1+x[k]; 
			logU0[j]=log(Uj0);
		}

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
		
		// the lists of nodes at each time s
		for(int t=0;t<=m;t++)
		nodeList[t]=new std::list<Node*>(); // empty lists
		
		// write the deterministic drifts mu_j(t)=mu_j(0,T_t)
		for(int j=1;j<n;j++)
		for(int t=1;t<=j;t++) 
			mu(j,t)=-0.5*factorLoading->integral_sgi_sgj_rhoij(j,j,0.0,Tc[t]);

		setRank2CovariationMatrixRoots(m);
        setTickSize(m);
		buildLattice(m);
		
	} // end constructor
	
	
	/** Goes through all the nodes in the lattice and checks if the transition 
	 *  probabilities are in [0,1] and sum to 1.
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
						
			std::list<Node*>& nodes_t=*nodeList[t];
			// iterate over list elements
			std::list<Node*>::const_iterator it;
			for(it=nodes_t.begin(); it!=nodes_t.end(); ++it){
				
				Node* node=*it;
				Matrix<Edge*>& edges=node->edges;
				
				psum=0.0;
		    	for(int k=-1;k<2;k++)
			    for(int l=-1;l<2;l++)
				if((k!=0)||(l!=0))
				{	
					Edge* edge=edges(k,l);
					p=edge->probability;
					if((p<0)||(p>1)) prob_failure=true;
					psum+=p;
				}
		        
				if(abs(psum-1.0)>0.0000001) prob_failure=true;
				if(prob_failure){
									
		         	cout << " *** t = " << t 
			             << "; transition probabilities defective, sum = " << psum;
					node->printTransitionProbabilities();
					cout << "\n\nQ_11=" << Q11[t] << ", Q_22=" << Q22[t] << ", Q_12=" << Q12[t];
				    exit(0);
				 } // endif
				
			} // end for it
							
            cout << " OK";
			
		} // end for t
		
	} // end selfTest()
				
	
	
	
	
	
// TEST
	
	static void test(int n)
    {
		Timer watch; watch.start();
		LiborFactorLoading* fl=JR_FactorLoading::sample(n);
		LiborTree2F Tree(fl,n-3);
		Tree.selfTest();
		watch.stop();
		watch.report("Time");
	}
		
		

private:
	
	/** <p>Computes the vector of matrices R=(R(0),R(1),...,R(m)), where R(t) is 
	 *  the rank 2 approximate root of the covariation matrix C(t) of the vector
	 *  \f$V(t)=(V_t(T_t),...,V_{n-1}(T_t))\f$. This matrix is needed to compute
	 *  the vector V(t) from the variables V1, V2 at the nodes at time t.
	 *
	 *  <p>Also computes the sequence of matrices RQ[t] (the 2 by 2 inverse of the last 
	 *  two rows of R[t]). Assumes \f$m\leq n-2\f$.
	 */
	void setRank2CovariationMatrixRoots(int m)
    {
	    for(int t=1;t<=m;t++){
		
			// R[t]
			Real* Tc=factorLoading->getTenorStructure(); 
			Real T_t=Tc[t];       // T_t
			
			// this is indexed as i,j=t,...,n-1
			UTRMatrix<Real>& Ct=factorLoading->logLiborCovariationMatrix(t,n,0.0,T_t);
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
	
	
	/** Computes the tick size a>0 and sets the corresponding quantities Q_ij(t).
	 * @param m lattice is built up to time t=m inclusively.
	 */
	void setTickSize(int m)
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
	
	
	// build lattice up to discrete time t=m (inclusive)
	void buildLattice(int m)
    {
		cout << "\n\n\nBuilding lattice until time t = " << m << endl << endl;
		  
		std::list<Node*>& nodes_0=*(nodeList[0]);               // list of nodes at time t
		
		// the initial node
		Node* NextNode=new Node(this,0,0,0);
		NextNode->V1=0.0;
		NextNode->V2=0.0;
	    // enter into node list at time t=0
		nodes_0.push_back(NextNode); 
		nodes++;
		
		// loop over times s=0,...,t-1
		for(int t=0;t<m;t++) {
	
			Real A1=k1[t]*a, A2=k2[t]*a;                          // tick sizes at time t
			std::list<Node*>& listCurrentNodes=*(nodeList[t]);    // list of nodes at time t+1
			std::list<Node*>& listNewNodes=*(nodeList[t+1]);      // list of nodes at time t+1

			// registry for possible nodes at time t+1, n(t+1,i,j)
			// i=-K_1(t+1),...,K_1(t+1), j=-K_2(t+1),...,K_2(t+1)
			Matrix<Node*> NewNodes(2*K1[t+1]+1,2*K2[t+1]+1,-K1[t+1],-K2[t+1]);
			
			// run through the list of nodes at time t and connect to the new nodes
			std::list<Node*>::const_iterator nds;
			for(nds=listCurrentNodes.begin(); nds!=listCurrentNodes.end(); ++nds)			
			{
				Node* CurrentNode=*nds;
				// state of this node
				Real V1=CurrentNode->V1, V2=CurrentNode->V2;
				int i=CurrentNode->i, j=CurrentNode->j;
				
				// loop over transitions k,l=-1,0,1
				for(int k=-1;k<2;k++)
				for(int l=-1;l<2;l++) if(!(k==0&&l==0)) 
				{	
					// connect to node n(s+1,i+k,j+l) in registry
					NextNode=NewNodes(i+k*k1[t],j+l*k2[t]);
					// if it does not exist yet
					if(!NextNode){ 
						
						  NextNode=new Node(this,t+1,i+k*k1[t],j+l*k2[t]); 
						  NewNodes(i+k*k1[t],j+l*k2[t])=NextNode;             // register the node
						  listNewNodes.push_back(NextNode);                   // append to list
						  NextNode->V1=V1+k*A1;
						  NextNode->V2=V2+l*A2;
						  NextNode->setStateVariables();
						  nodes++; 
					}
					
					Edge* edge=CurrentNode->edges(k,l);
					edge->node=NextNode;
					edge->probability=probability(t,k,l);

				} // end for k,l
			} // end for nds

				
			cout << "\nTotal nodes = " << nodes;
		
		} // end for t
		
	} // end buildLattice
	
	
	// the transition probability p_{ij}
	Real probability(int t, int i, int j)
    {
		if(j==0) return 0.5*(1-Q22[t]);
		if(i==0) return 0.5*(1-Q11[t]);
		if(i+j==0) return 0.25*(Q11[t]+Q22[t]-Q12[t]-1);
		else return 0.25*(Q11[t]+Q22[t]+Q12[t]-1);
	}
	
		
						
// EDGES AND NODES


	
	struct Edge {
	
	    Node* node;              // node we transition to along this edge	
	    Real probability;        // transition probability

	}; // end edge	
	

	/** Single node containing the vector \f$V(T_t)=(V_t(T_t),...,V_{n-1}(T_t))\f$
	 *  as well as the variables \f$V1=V_{n-1}, V2=V_{n-2}-V_{n-1}\f$ which the lattice 
	 *  evolves.
	 */
    struct Node {
	
    	LiborTree2F* theTree;   // the tree the node lives in
		
		Real V1, V2;            // the variables which are evolved
		int n,                  // n dimension of the underlying Libor process.     
		    t,                  // discrete time at which the node lives
		    i,                  // V2=V2(0)+ja
		    j;                  // V1=V1(0)+ia, 
		
		vector<Real> H;         // the accrual factors H_j(T_t), j=t,...,n at this node
		                        // natural indexation, index base t
	    Real  F;                // value of some path functional at this node
	
	    Matrix<Edge*> edges;    // edges with transition probabilities p_ij, i,j=-1,0,1.
		                        // natural indexation, index base -1 
		

		/** @param tree the Libor tree the node lives in.
		 *  @param s discrete time t (continuous time T_t) at which the node lives.
		 *  @param k state V_1=V_1(0)+ka
		 *  @param l state V_2=V_2(0)+la
		 */
		Node(LiborTree2F* tree, int s, int k, int l) : 
		theTree(tree),
		n(tree->getDimension()), t(s), i(k), j(l), 
		H(n-t+1,t), edges(3,3,-1,-1) 
	    {
		     // loop over transitions k,l=-1,0,1
			 for(int k=-1;k<2;k++)
			 for(int l=-1;l<2;l++) edges(l,k)=new Edge();
		}
				 
		
		~Node()
	    {
			for(int k=-1;k<2;k++)
			for(int l=-1;l<2;l++) delete edges(k,l); 
				
	 	}
		
		
		/** Prints the transition probabilities and corresponding Q_ij(t) values.
		 *  Diagnostic.
		 */
		void printTransitionProbabilities()
	    {
			 for(int k=-1;k<2;k++)
			 for(int l=-1;l<2;l++)
			 if((k!=0)||(l!=0))
		     {	
				  Edge* edge=edges(k,l);
			      Real p=edge->probability;
			      cout << "\np(" << k << "," << l << ") = " << p;
		     }		
		 } // end printTransitionProbabilities

		 
		 /** Computes the vector V from the variables V1,V2 which are evolved in the 
		  *  lattice. Needs the vector logU0[j]=log(U_j(0)) and the matrices Rt=*R[t], 
		  *  RQt=*RQ[t] from the LiborTree2F. Unfortunately we have to hand these in as 
		  *  parameters since the class Node is unaware
		  *
		  */
		 void setStateVariables()
	     {
             vector<Real> V(n-t,t);   // the volatility parts V_j of log(U_j)
			                          // j=t,...,n-1
			 V[n-1]=V1;
			 V[n-2]=V1+V2;
			 
			 // matrices needed to derive V_j, j<n-2, from V_{n-2}, V_{n-1}
			 Matrix<Real>& Rt=theTree->getR(t);
			 Matrix<Real>& RQt=theTree->getRQ(t);
			 
			 // the "factors" Z_1, Z_2 for V_{n-2}, V_{n-1}
			 vector<Real> Z(2);   
			 Z[0]=V[n-2]; Z[1]=V[n-1]; Z*=RQt;
		 
			 // use these to compute the remaining V_j
			 // Rt=*R[t] has row index base t.
			 for(int j=t;j<n-2;j++) V[t]=Rt(t,0)*Z[0]+Rt(t,1)*Z[1];
				 
			 // add the initial value log(U_j(0) and drift mu_j(t)=mu_j(0,T_t) 
			 // to obtain the log(U_j(T_t))
			 vector<Real>& logU0=theTree->getLogU0();
			 LTRMatrix<Real>& mu=theTree->getDrifts();
			 for(int j=t;j<n;j++) V[j]+=logU0[j]+mu(j,t);
				 
			 // move from log(U_j) to U_j
			 for(int j=t;j<n;j++) V[j]=exp(V[j]);    // now V=U
				 
			 // write the H_j=U_j+H_{j+1}, H_n=1, j=t,...,n
			 H[n]=1.0; Real f=1.0;
	         for(int j=n-1;j>=t;j--){ f+=V[j]; H[j]=f; }
			 
		 }

		
		

	
    }; // end Node
	
	
}; // end LiborTree2F



	
	
	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
