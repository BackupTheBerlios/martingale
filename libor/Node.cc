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
#include "Lmmlattice.h"
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
	std::list<Edge>::const_iterator theEdge;         // pointer to Edge
	for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge)  delete &(*theEdge); 
}


void 
Node::
printTransitionProbabilities() const
{
	cout << "\n\n\nTransition probabilities\n: ";
	std::list<Edge>::const_iterator theEdge;
	for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge) 
	{	
		Real p=theEdge->probability;
		cout << p << ", ";
	}		
} // end printTransitionProbabilities



/**********************************************************************************
 *
 *            BASE TYPES FOR LMM NODES
 *
 *********************************************************************************/

// LITE-LMM-NODES

// out of class initialization necessary
LiteLmmNode::H_(LMM_MAX_DIM); 
LiteLmmNode::V_(LMM_MAX_DIM);


LiteLmmNode::
LiteLmmNode(int s, const IntArray1D& k, int LmmLatticeData* latticeData) : 
s_(s), lattice(latticeData), k_(new int[info->r]) 
{  
	// set the state
	int r=lattice->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}

	
LiteLmmNode::
const RealArray1D& Hvect(int p)
{
	int n=lattice->n,                          // dimension of Libor process
	    r=lattice->r;                          // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;               // tick size a over one time step
	const RealArray1D& U0=lattice->U0;         // U_j(0), initial values
	const RealArray1D& sg=lattice->sg;         // constant vols sigma_j, see book, 6.11.1
	const RealArray1D& mu=lattice->driftUnit;  // drifts over one time step
	const RealMatrix& R=lattice->R;            // r factor pseudo square root of correlation matrix

	
	// the volatility parts V_j of log(U_j), j=t,...,n-1
	V_.setDimension(n-p);
	V_.setIndexBase(p);
	for(int j=p;j<n;j++){
		
		V_[j]=0.0;
		for(int u=0;u<r;u++) V_[j]+=R(j,u)*k[u];
		V_[j]*=a*sg[j];
	}
		 
    // add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=p;j<n;j++) V_[j]+=log(U0[j])+s*mu[j];         // now V=log(U)
				 
	// move from log(U_j) to U_j
	for(int j=p;j<n;j++) V_[j]=std::exp(V_[j]);              // now V=U
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1
	// into the static workspace H_:
	H_.setDimension(n-p+1);
	H_.setIndexBase(p);
	H_[n]=1.0; 
    for(int j=n-1;j>=p;j--) H_[j]=V_[j]+H_[j+1]; 

    return H_;
	
} // Hvect


LiteLmmNode:: 
int get_t() const
{
	int nSteps=lattice->nSteps;
	if(s%nSteps==0) return s/nSteps;
	return s/nSteps+1;
}


LiteLmmNode:: 
std::ostream& printType(std::ostream& os)
{
	return
	os << "Lightweight LmmNode";
}


LiteLmmNode:: 
std::ostream& printSelf(std::ostream& os)
{
	int t=get_t();
	return
	os << "\n\nLite LmmNode, factors: " << lattice->r
	   << "\nAccrual interval: (T_"<<t-1<<",T_"<<t<<"], time step s = " << s
	   << "\nVector of accrual factors:" << Hvect(t);
}




// HEAVY-LMM-NODES

// out of class initialization necessary
HeavyLmmNode::V_(LMM_MAX_DIM);
	
HeavyLmmNode::
HeavyLmmNode(int s, const IntArray1D& k, LmmLatticeData* latticeData) : 
s_(s), k_(new int[info->r]), lattice(latticeData), H(n-get_t()+1,get_t())
{  
	int n=lattice->n,                          // dimension of Libor process
	    r=lattice->r;                          // number of factors
	Real dt = lattice->timestep,         
	     a  = lattice->ticksize;               // tick size a over one time step
	const RealArray1D& U0=lattice->U0;         // U_j(0), initial values
	const RealArray1D& sg=lattice->sg;         // constant vols sigma_j, see book, 6.11.1
	const RealArray1D& mu=lattice->driftUnit;  // drifts over one time step
	const RealMatrix& R=lattice->R;            // r factor pseudo square root of correlation matrix

	
	// set the state of the Brownian driver
	for(int i=0;i<r;i++) k_[i]=k[i];

	// compute the accrual factors H_j
	
	int t=get_t();             // node lives in (T_{t-1},T_t].
	// the volatility parts V_j of log(U_j), j=t,...,n-1.
	RealVector V(n-t,t);       
	for(int j=t;j<n;j++){
		
		V[j]=0.0;
		for(int u=0;u<r;u++) V[j]+=(R(j,0)*k_[u];
		V[j]*=sg[j]*a;
	}
		 
	// add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=t;j<n;j++) V[j]+=log(U0[j])+s*mu[j];           // now V=log(U)      
				 
	// move from log(U_j) to U_j
	for(int j=t;j<n;j++) V[j]=std::exp(V[j]);                 // now V=U
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
	H[n]=1.0; 
    for(int j=n-1;j>=t;j--) H[j]=V[j]+H[j+1]; 

} // end constructor


HeavyLmmNode:: 
int get_t() const
{
	int nSteps=lattice->nSteps;
	if(s%nSteps==0) return s/nSteps;
	return s/nSteps+1;
}


HeavyLmmNode:: 
std::ostream& printType(std::ostream& os)
{
	return
	os << "Heavyweight LmmNode";
}


HeavyLmmNode:: 
std::ostream& printSelf(std::ostream& os)
{
	int t=get_t();
	return
	os << "\n\nHeavy LmmNode, factors: " << lattice->r
	   << "\nAccrual interval: (T_"<<t-1<<",T_"<<t<<"], time step s = " << s
	   << "\nVector of accrual factors:" << H;
}




/**********************************************************************************
 *
 *            ASSET BASKET NODES 
 *
 *********************************************************************************/


// LITE-BASKET-NODES

// out of class initialization necessary
LiteBaskeNode::S_(BASKET_MAX_DIM); 
LiteBasketNode::V_(BASKET_MAX_DIM);


LiteBaskeNode::
LiteBaskeNode(int s, const IntArray1D& k, int BasketLatticeData* latticeData) : 
s_(s), lattice(latticeData), k_(new int[info->r]) 
{  
	// set the state
	int r=lattice->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}

	
LiteBaskeNode::
const RealArray1D& Svect(int p)
{
	int n=lattice->n,                          // dimension of Libor process
	    r=lattice->r;                          // number of factors
	Real dt = lattice->timestep;         
	     a  = lattice->ticksize;               // tick size a over one time step
	const RealArray1D& S0=lattice->S0;         // initial asset prices
	const RealArray1D& sg=lattice->sg;         // volatilities, see ConstVolFactorLoading.
	const RealArray1D& mu=lattice->driftUnit;  // drifts over one time step
	const RealMatrix& R=lattice->R;            // r factor pseudo square root of correlation matrix

	
	// the volatility parts V_j of log(S_j), j=p,...,n-1
	V_.setDimension(n-p);
	V_.setIndexBase(p);
	for(int j=p;j<n;j++){
		
		V_[j]=0.0;
		for(int u=0;u<r;u++) V_[j]+=R(j,u)*k[u];
		V_[j]*=a*sg[j];
	}
		 
    // add the initial value log(S_j(0) and drift mu_j(s)
	for(int j=p;j<n;j++) V_[j]+=log(S0[j])+s*mu[j];          // now V=log(S)
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1
	// into the static workspace H_:
	S_.setDimension(n-p);
	S_.setIndexBase(p);
    for(int j=n-p;j<p;j++) S_[j]=std::exp(V_[j]); 

    return S_;
	
} // Svect


LiteBasketNode:: 
std::ostream& printType(std::ostream& os)
{
	return
	os << "Lightweight BasketNode";
}


LiteBaskeNode:: 
std::ostream& printSelf(std::ostream& os)
{
	return
	os << "\n\nLightweight BasketNode, factors: " << lattice->r
	   << "\nVector of assets factors:" << Svect(0);
}




// HEAVY-BASKET-NODES

// out of class initialization necessary
HeavyBasketNode::V_(BASKET_MAX_DIM);
	

	
HeavyBasketNode::
HeavyBasketNode(int s, const IntArray1D& k, BasketLatticeInfo* info) : 
s_(s), k_(new int[info->r]), lattice(info), S(n)
{  
	int n=lattice->n,                   // dimension of Libor process
	    r=lattice->r;                   // number of factors
	Real dt = lattice->timestep;         
	     a  = lattice->ticksize;               // tick size a over one time step
	const RealArray1D& S0=lattice->S0;         // initial asset prices
	const RealArray1D& sg=lattice->sg;         // volatilities, see ConstVolFactorLoading.
	const RealArray1D& mu=lattice->driftUnit;  // drifts over one time step
	const RealMatrix& R=lattice->R;            // r factor pseudo square root of correlation matrix

		
	// set the state of the Brownian driver
	for(int i=0;i<r;i++) k_[i]=k[i];

	// compute the assets S_j
    
	// volatility parts V_j of log(S_j)
	for(int j=0;j<n;j++){
		
		V[j]=0.0;
		for(int u=0;u<r;u++) V[j]+=(R(j,0)*k_[u];
		V[j]*=sg[j]*a;
	}
		 
	// add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=t;j<n;j++) V[j]+=log(U0[j])+s*mu[j];           // now V=log(U)      
				 
	// move from log(U_j) to U_j
	for(int j=t;j<n;j++) V[j]=std::exp(V[j]);                 // now V=U
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1, 
	H[n]=1.0; 
    for(int j=n-1;j>=t;j--) H[j]=V[j]+H[j+1]; 

} // end constructor


HeavyBasketNode:: 
std::ostream& printType(std::ostream& os)
{
	return
	os << "Heavyweight BasketNode";
}


HeavyBasketNode:: 
std::ostream& printSelf(std::ostream& os)
{
	int t=get_t();
	return
	os << "\n\nHeavyweight BaskeNode, factors: " << lattice->r
	   << "\nAsset vector:" << S;
}



MTGL_END_NAMESPACE(Martingale)	

