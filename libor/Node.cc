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

#include "Node.h"
#include "LiborFactorLoading.h"
#include "LmmLattice.h"
#include "BasketLattice.h"
#include "Matrix.h"
#include <iostream>
#include <cstdlib>                        // exit()



MTGL_BEGIN_NAMESPACE(Martingale)



/**********************************************************************************
 *
 *            GENERAL NODES 
 *
 *********************************************************************************/
 
Node::
~Node()
{ 
	vector<Edge*>::iterator theEdge=edges.begin();
	while(theEdge!=edges.end()){ delete *theEdge; ++theEdge; }
	// the direct member edges is deallocated automatically
} 	


void 
Node::
printTransitionProbabilities() const
{
	cout << "\n\n\nTransition probabilities\n: ";
	vector<Edge*>::const_iterator theEdge=edges.begin();
	while(theEdge!=edges.end()) {
	
		Edge* edge = *theEdge;
		cout << edge->probability << ", ";
		++theEdge;
	}		
} // end printTransitionProbabilities



/**********************************************************************************
 *
 *            BASE TYPES FOR LMM NODES
 *
 *********************************************************************************/

// LITE-LMM-NODES

// out of class initialization necessary
RealArray1D LmmNode_LiteBase::H_(LMM_MAX_DIM); 
RealArray1D LmmNode_LiteBase::V_(LMM_MAX_DIM);


LmmNode_LiteBase::
LmmNode_LiteBase(int s, const IntArray1D& k, LmmLatticeData* latticeData) : Node(s), 
lattice(latticeData), k_(new int[latticeData->r]) 
{
	// set the state
	int r=latticeData->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}

	

const RealArray1D& 
LmmNode_LiteBase::
Hvect(int p)
{
	int n=lattice->n,                           // dimension of Libor process
	    r=lattice->r;                           // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;                // tick size a over one time step
	const RealArray1D& log_U0=lattice->log_U0;  // U_j(0), initial values
	const RealArray1D& sg=lattice->sg;          // constant vols sigma_j, see book, 6.11.1
	const RealArray1D& mu=lattice->driftUnit;   // drifts over one time step
	const RealMatrix& R=lattice->R;             // r factor pseudo square root of correlation matrix

	
	if(s==0){          // needs special treatment since R starts with row index 1
		
		H_[n]=1.0;
        for(int j=n-1;j>=0;j--) H_[j]=std::exp(log_U0[j])+H_[j+1]; 
		return H_;
	}
			
    // the volatility parts V_j of log(U_j), j=t,...,n-1
    V_.setDimension(n-p);
	V_.setIndexBase(p);
	for(int j=p;j<n;j++){
		
	    V_[j]=0.0;
		for(int u=0;u<r;u++) V_[j]+=R(j,u)*k_[u];
	    V_[j]*=a*sg[j];
	}                                                       // now V_=V
		 
    // add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=p;j<n;j++) {
		
        V_[j]+=log_U0[j]+s*mu[j];          // now V=log(U)
		V_[j]=std::exp(V_[j]);             // now V=U
	}
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1
	// into the static workspace H_:
	H_.setDimension(n-p+1);
	H_.setIndexBase(p);
	H_[n]=1.0; 
    for(int j=n-1;j>=p;j--) H_[j]=V_[j]+H_[j+1]; 

    return H_;
	
} // Hvect



int 
LmmNode_LiteBase:: 
get_t() const
{
	int nSteps=lattice->nSteps;
	if(s%nSteps==0) return s/nSteps;
	return s/nSteps+1;
}



std::ostream& 
LmmNode_LiteBase:: 
printType(std::ostream& os) 
{
	return
	os << "Lightweight LmmNode";
}



std::ostream& 
LmmNode_LiteBase:: 
printSelf(std::ostream& os)
{
	int t=get_t();
	return
	os << "\n\nLite LmmNode, factors: " << lattice->r
	   << "\nAccrual interval: (T_"<<t-1<<",T_"<<t<<"], time step s = " << s
	   << "\nVector of accrual factors:" << Hvect(t);
}




// HEAVY-LMM-NODES

// out of class initialization necessary
RealArray1D LmmNode_HeavyBase::V_(LMM_MAX_DIM);
	
LmmNode_HeavyBase::
LmmNode_HeavyBase(int s, const IntArray1D& k, LmmLatticeData* latticeData) : Node(s),
lattice(latticeData),                        // get_t() needs this for nSteps
k_(new int[latticeData->r]), 
H_((latticeData->n)-get_t()+1,get_t()) 
{ 
	int n=lattice->n,                           // dimension of Libor process
	    r=lattice->r;                           // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;                // tick size a over one time step
	const RealArray1D& log_U0=lattice->log_U0;  // log(U_j(0)), initial values
	const RealArray1D& sg=lattice->sg;          // constant vols sigma_j, see book, 6.11.1
	const RealArray1D& mu=lattice->driftUnit;   // drifts over one time step
	const RealMatrix& R=lattice->R;             // r factor pseudo square root of correlation matrix
	
	// set the state of the Brownian driver
	for(int i=0;i<r;i++) k_[i]=k[i];

	// compute the accrual factors H_j
	
	int t=get_t();             // node lives in (T_{t-1},T_t].
	if(t==0){                  // needs special treatment since R starts with row index 1
		
		H_[n]=1.0;
        for(int j=n-1;j>=t;j--) H_[j]=std::exp(log_U0[j])+H_[j+1]; 
			
	} else {                   // computation from the dynamics
	   
		// the volatility parts V_j of log(U_j), j=t,...,n-1.  
	    for(int j=t;j<n;j++){
		
	    	V_[j]=0.0;
		    for(int u=0;u<r;u++) V_[j]+=R(j,0)*k_[u];
		    V_[j]*=sg[j]*a;
	    }                                        // now V_=V
	 
	    // add the initial value log(U_j(0) and drift mu_j(s)
	    for(int j=t;j<n;j++){
		
		    V_[j]+=log_U0[j]+s*mu[j];             // now V=log(U)      
		    V_[j]=std::exp(V_[j]);                // now V=U
	    }
			 
	    // write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
	    H_[n]=1.0; 
        for(int j=n-1;j>=t;j--) H_[j]=V_[j]+H_[j+1];
	} // end else

} // end constructor



int 
LmmNode_HeavyBase:: 
get_t() const
{
	int nSteps=lattice->nSteps;
	if(s%nSteps==0) return s/nSteps;
	return s/nSteps+1;
}



std::ostream& 
LmmNode_HeavyBase:: 
printType(std::ostream& os) 
{
	return
	os << "Heavyweight LmmNode";
}



std::ostream& 
LmmNode_HeavyBase:: 
printSelf(std::ostream& os) const
{
	int t=get_t();
	return
	os << "\n\nHeavy LmmNode, factors: " << lattice->r
	   << "\nAccrual interval: (T_"<<t-1<<",T_"<<t<<"], time step s = " << s
	   << "\nVector of accrual factors:" << H_;
}




/**********************************************************************************
 *
 *            ASSET BASKET NODES 
 *
 *********************************************************************************/


// LITE-BASKET-NODES

// out of class initialization necessary
RealArray1D LiteBasketNode::S_(BASKET_MAX_DIM); 
RealArray1D LiteBasketNode::V_(BASKET_MAX_DIM);


LiteBasketNode::
LiteBasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData) : Node(s),
lattice(latticeData), k_(new int[latticeData->r])
{  
	// set the state
	int r=lattice->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}

	

const RealArray1D& 
LiteBasketNode::
Svect(int p)
{
	int n=lattice->n,                           // dimension of Libor process
	    r=lattice->r;                           // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;                // tick size a over one time step
	const RealArray1D& log_S0=lattice->log_S0;  // initial asset price logs
	const RealArray1D& sg=lattice->sg;          // volatilities, see ConstVolFactorLoading.
	const RealArray1D& mu=lattice->driftUnit;   // drifts over one time step
	const RealMatrix& R=lattice->R;             // r factor pseudo square root of correlation matrix

	
	// the volatility parts V_j of log(S_j), j=p,...,n-1
	S_.setDimension(n-p);
	S_.setIndexBase(p);
	for(int j=p;j<n;j++){
		
		S_[j]=0.0;
		for(int u=0;u<r;u++) S_[j]+=R(j,u)*k_[u];
		S_[j]*=a*sg[j];
	}                                                       // now S_=V
		 
    // add the initial value log(S_j(0) and drift mu_j(s)
	for(int j=p;j<n;j++){
		
		S_[j]+=log_S0[j]+s*mu[j];           // now S_=log(S)
        S_[j]=std::exp(S_[j]);              // now S_=S
	}

    return S_;
	
} // Svect



std::ostream& 
LiteBasketNode:: 
printType(std::ostream& os) 
{
	return
	os << "Lightweight BasketNode";
}



std::ostream& 
LiteBasketNode:: 
printSelf(std::ostream& os) 
{
	return
	os << "\n\nLightweight BasketNode, factors: " << lattice->r
	   << "\nVector of assets factors:" << Svect(0);
}




// HEAVY-BASKET-NODES

// out of class initialization necessary
RealArray1D HeavyBasketNode::V_(BASKET_MAX_DIM);
	
HeavyBasketNode::
HeavyBasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData) : Node(s), 
lattice(latticeData), k_(new int[latticeData->r]), S_(latticeData->n)
{  
	int n=lattice->n,                            // dimension of Libor process
	    r=lattice->r;                            // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;                 // tick size a over one time step
	const RealArray1D& log_S0=lattice->log_S0;   // initial asset price logs
	const RealArray1D& sg=lattice->sg;           // volatilities, see ConstVolFactorLoading.
	const RealArray1D& mu=lattice->driftUnit;    // drifts over one time step
	const RealMatrix& R=lattice->R;              // r factor pseudo square root of correlation matrix

		
	// set the state of the Brownian driver
	for(int i=0;i<r;i++) k_[i]=k[i];

	// compute the assets S_j
    
	// volatility parts V_j of log(S_j)
	for(int j=0;j<n;j++){
		
		S_[j]=0.0;
		for(int u=0;u<r;u++) S_[j]+=R(j,0)*k_[u];
		S_[j]*=sg[j]*a;
	}                                                          // now S_=V
		 
	// add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=0;j<n;j++){
		
		S_[j]+=log_S0[j]+s*mu[j];              // now S_=log(S)      
		S_[j]=std::exp(S_[j]);                 // now S_=S
	}

} // end constructor



std::ostream& 
HeavyBasketNode:: 
printType(std::ostream& os) 
{
	return
	os << "Heavyweight BasketNode";
}



std::ostream& 
HeavyBasketNode:: 
printSelf(std::ostream& os) const
{
	return
	os << "\n\nHeavyweight BaskeNode, factors: " << lattice->r
	   << "\nAsset vector:" << S_;
}



MTGL_END_NAMESPACE(Martingale)	

