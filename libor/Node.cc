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
#include "LiborFunctional.h"
#include <iostream>
#include <cstdlib>                        // exit()
#include <cmath>


using std::cout;
using std::exp;
using std::ostream;


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
 *             LMM NODES
 *
 *********************************************************************************/



// out of class initialization necessary
RealArray1D LmmNode::H_(LMM_MAX_DIM); 
RealArray1D LmmNode::V_(LMM_MAX_DIM);


LmmNode::
LmmNode(int s, const IntArray1D& k, LmmLatticeData* latticeData) : Node(s), 
lattice(latticeData), k_(new int[latticeData->r]) 
{
	// set the state
	int r=latticeData->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}



int 
LmmNode:: 
get_t() const
{
	int nSteps=lattice->nSteps;
	if(s%nSteps==0) return s/nSteps;
	return s/nSteps+1;
}


const RealArray1D& 
LmmNode::
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
	
	H_.setDimension(n-p+1);
	H_.setIndexBase(p);	
	if(s==0){          // needs special treatment since R starts with row index 1
		
		H_[n]=1.0;
        for(int j=n-1;j>=p;j--) H_[j]=std::exp(log_U0[j])+H_[j+1]; 
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
		V_[j]=std::exp(V_[j]);                  // now V=U
	}
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1
	// into the static workspace H_:
	H_.setDimension(n-p+1);
	H_.setIndexBase(p);
	H_[n]=1.0; 
    for(int j=n-1;j>=p;j--) H_[j]=V_[j]+H_[j+1]; 

    return H_;
	
} // Hvect



Real 
LmmNode:: 
H_pq(int p, int q)
{
	const RealArray1D& H=Hvect(p);
	Real delta = lattice->delta;
	return LiborFunctional::H_pq(p,q,H,delta);
}
	
	

Real 
LmmNode:: 
swapRate(int p, int q)
{
	const RealArray1D& H=Hvect(p);
	Real delta = lattice->delta;
	return LiborFunctional::swapRate(p,q,H,delta);
}
	
	

Real 
LmmNode:: 
forwardSwaptionPayoff(int p, int q, Real kappa)
{
	const RealArray1D& H=Hvect(p);
	Real delta = lattice->delta;
	return LiborFunctional::forwardSwaptionPayoff(p,q,kappa,H,delta);
}


Real 
LmmNode:: 
forwardCapletPayoff(Real kappa)
{
	int i=get_t();
	const RealArray1D& H=Hvect(i);
	Real delta = lattice->delta;
	return LiborFunctional::forwardCapletPayoff(i,kappa,H,delta);
}


Real 
LmmNode:: 
forwardBondCallPayoff(BondCall* bc)
{
	int t=get_t();
	const RealArray1D& H=Hvect(t);
    return LiborFunctional::forwardBondCallPayoff(bc,H);
}


ostream& 
LmmNode:: 
printSelf(ostream& os)
{
	int t=get_t();
	return
	os << "\n\nLmmNode, factors: " << lattice->r
	   << "\nAccrual interval: (T_"<<t-1<<",T_"<<t<<"], time step s = " << s
	   << "\nVector of accrual factors:" << Hvect(t);
}






/**********************************************************************************
 *
 *            ASSET BASKET NODES 
 *
 *********************************************************************************/


// LITE-BASKET-NODES

// out of class initialization necessary
RealArray1D BasketNode::S_(BASKET_MAX_DIM); 
RealArray1D BasketNode::V_(BASKET_MAX_DIM);


BasketNode::
BasketNode(int s, const IntArray1D& k, BasketLatticeData* latticeData) : Node(s),
lattice(latticeData), k_(new int[latticeData->r])
{  
	// set the state
	int r=lattice->r;  // number of factors
	for(int i=0;i<r;i++) k_[i]=k[i];
}

	

const RealArray1D& 
BasketNode::
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
        S_[j]=std::exp(S_[j]);                   // now S_=S
	}

    return S_;
	
} // Svect



ostream& 
BasketNode:: 
printSelf(ostream& os) 
{
	return
	os << "\n\nBasketNode, factors: " << lattice->r
	   << "\nVector of assets factors:" << Svect(0);
}




MTGL_END_NAMESPACE(Martingale)	

