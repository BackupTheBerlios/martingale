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
#include "Node.h"
#include "Lattice.h"
#include "Utils.h"
#include "Matrix.h"
#include "LiborFactorLoading.h" 
#include "LiborFunctional.h"
#include "Option.h"
#include <cmath>


MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::log;


// LMM LATTICE DATA



// GENERAL LMM LATTICE

LmmLattice::
LmmLattice(int q, LiborFactorLoading* fl, int t, int steps, bool verbose) : 
// m=t*steps is the number of time steps
Lattice<StandardBrownianNode>(t*steps),
factorLoading(fl),
n(fl->getDimension()), 
r(q),
nSteps(steps),
delta(fl->getDeltas()[0]),
dt(delta/steps),
a(sqrt(dt)),
sg(n-1,1),
log_U0(n),
mu(n-1,1),
R(fl->getRho().rankReducedRoot(r))
{ 
	// check the number of factors
	if((r!=2)&&(r!=3)){
		
	   cout << "\n\nLmmLattice constructor: number of factors must be 2 or 3 but is " << r
	             << "\nTerminating.";
	   exit(1);
    }		
		
	// check if Libor accrual periods are constant
	const RealArray1D& deltas=fl->getDeltas();
	for(int j=0;j<n;j++) if(deltas[j]!=delta) {
			
	   cout << "\n\nLmmLatticeconstructor: Libor accrual periods not constant."
	             << "\nTerminating.";
	   exit(1);
    }
	
	// check if volatilities are constant
	if(fl->getType()->volType!=VolSurface::CONST) {
			
	   cout << "\n\nLmmLattice-constructor: volatility surface not constant."
	             << "\nTerminating.";
	   exit(1);
    }	
		
	// set log(U_j(0)), j=0,1,...,n-1
    const RealArray1D& x=factorLoading->getInitialXLibors();     // x[j]=X_j(0)
    for(int j=0;j<n;j++){ 
			
	    // U_j(0)=X_j(0)(1+X_{j+1}(0))...(1+X_{n-1}(0))
		Real Uj0=x[j]; for(int k=j+1;k<n;k++) Uj0*=1+x[k]; 
		log_U0[j]=log(Uj0);
	}

    // set constant vols
	for(int j=1;j<n;j++) sg[j]=factorLoading->sigma(j,0.0);
	
	// drifts of Y_j=log(U_j) over a single time step
    for(int i=1;i<n;i++) mu[i]=-sg[i]*sg[i]*dt/2; 
	
	// m=t*steps is the number of time steps
	buildLattice(t*steps,verbose);
	
} // end constructor



LmmLattice* 
LmmLattice::
sample(int r, int n, int p, int nSteps, bool verbose=false) 
{
	LiborFactorLoading*
	fl=LiborFactorLoading::sample(n,VolSurface::CONST,Correlations::CS);
	LmmLattice* lattice = new LmmLattice(r,fl,p,nSteps,verbose);
	return lattice;
}



void 
LmmLattice::
rescaleVols()
{
	for(int i=1;i<n;i++){
	
		Real f=0.0;      // norm of row_i(R)
		for(int j=0;j<r;j++) f+=R(i,j)*R(i,j); 
		f=sqrt(f);

		R.scaleRow(i,1.0/f);
	}
}



void 
LmmLattice::
testFactorization() const
{
	cout << "\n\nRelative errors of the approximate rank " << r << " factorization"
	     << "\nrho=RR' (rank(R)=" << r << ") of the correlation matrix rho"
	     << "\n(trace norm): " << endl << endl;
  
	factorLoading->getRho().testFactorization(r);

}   // end factorAnalysis


ostream& 
LmmLattice::
printSelf(ostream& os) const
{
	return
	os << r << "factor LMM lattice: "
	   << "\nNumber of time steps in each accrual interval: " << nSteps 
	   << endl << *factorLoading;
}


const RealArray1D& 
LmmLattice::
Hvect
(int p, StandardBrownianNode* node, int s)
{
	H_.setDimension(n-p+1);
	H_.setIndexBase(p);	
	H_[n]=1.0;
	// needs special treatment since R starts with row index 1
	if(s==0){          
		
        for(int j=n-1;j>=p;j--) H_[j]=std::exp(log_U0[j])+H_[j+1]; 
		return H_;
	}
			
    // the volatility parts V_j of log(U_j), j=t,...,n-1
	int* k=node->getIntegerTicks();                         // Z_j=k[j]*a
    V_.setDimension(n-p);
	V_.setIndexBase(p);
	for(int j=p;j<n;j++){
		
	    V_[j]=0.0;
		for(int u=0;u<r;u++) V_[j]+=R(j,u)*k[u];
	    V_[j]*=a*sg[j];
	}                                                       // now V_=V
		 
    // add the initial value log(U_j(0) and drift mu_j(s)
	for(int j=p;j<n;j++) {
		
        V_[j]+=log_U0[j]+s*mu[j];          // now V=log(U)
		V_[j]=std::exp(V_[j]);                  // now V=U
	}
				 
	// write the  H_n=1, H_j=U_j+H_{j+1}, j=t,...,n-1
	// into the static workspace H_:
    for(int j=n-1;j>=p;j--) H_[j]=V_[j]+H_[j+1]; 

    return H_;
	
} // Hvect



Real 
LmmLattice::
L(int j, StandardBrownianNode* node, int s)
{
	const RealArray1D& H=Hvect(j,node,s);
	return (H[j]/H[j+1]-1)/delta;
}
	
	
Real 
LmmLattice::
H_pq(int p, int q, StandardBrownianNode* node, int s)
{
	const RealArray1D& H=Hvect(p,node,s);
	return LiborFunctional::H_pq(p,q,H,delta);
}


Real 
LmmLattice::
swapRate(int p, int q, StandardBrownianNode* node, int s)
{
	const RealArray1D& H=Hvect(p,node,s);
	return LiborFunctional::swapRate(p,q,H,delta);
}
		

Real 
LmmLattice::
forwardSwaptionPayoff(int p, int q, Real kappa, StandardBrownianNode* node, int s)
{
	const RealArray1D& H=Hvect(p,node,s);
	return LiborFunctional::forwardSwaptionPayoff(p,q,kappa,H,delta);
}


Real 
LmmLattice::
forwardCapletPayoff(int i, Real kappa, StandardBrownianNode* node, int s)
{
	const RealArray1D& H=Hvect(i,node,s);
	return LiborFunctional::forwardCapletPayoff(i,kappa,H,delta);
}



Real 
LmmLattice::
forwardBondCallPayoff(BondCall* bc, StandardBrownianNode* node, int s)
{
	int t = bc->getPeriodsToExpiry();     // option expires at T_t 
	const RealArray1D& H=Hvect(t,node,s);
    return LiborFunctional::forwardBondCallPayoff(bc,H);
}



Real 
LmmLattice::
tau(int s)
{
    const RealArray1D& T_=factorLoading->getTenorStructure();
  	int t=s/nSteps;
   	Real delta_t=T_[t+1]-T_[t];
		
    return T_[t]+(delta_t*(s%nSteps))/nSteps;
}


void 
LmmLattice::
test(int r, int n)
{
    Timer watch; watch.start();
	LmmLattice* lattice = sample(r,n,n,1,true);
	lattice->rescaleVols();
	lattice->selfTest();
	delete lattice;
	watch.stop();
	watch.report("LmmLattice test");
}


void 
LmmLattice::
buildLattice(int m, bool verbose)
{
	if(verbose)
	cout << "\n\nBuilding lattice: " << *this
	          << "\n\nTime steps: " << m << endl << endl;
	switch(r){
		
		case 3 :  LatticeBuilder::buildThreeFactorLattice(this,m,verbose); 
			      break;
		default : LatticeBuilder::buildTwoFactorLattice(this,m,verbose); 
	}

} // end buildLattice



// out of class intialization
RealArray1D LmmLattice::H_(LMM_MAX_DIM); 
RealArray1D LmmLattice::V_(LMM_MAX_DIM);

	


MTGL_END_NAMESPACE(Martingale)

