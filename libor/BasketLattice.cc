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
#include "Matrix.h"
#include "Utils.h"
#include "Node.h"
#include "FactorLoading.h"
#include <cmath>


MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::log;
using std::exp;

 
	

// GENERAL BASKET LATTICE

void 
BasketLattice::
testFactorization() const 
{
	fl->getCorrelationMatrix().testFactorization
	(r_,"Relative error of the rank r factorization of the correlation matrix (trace norm)");
}
	
	
BasketLattice::
BasketLattice
(int r, ConstantFactorLoading* fl, int T, Real ds, RealArray1D S0) : 
Lattice<StandardBrownianNode>(T),
n_(fl->getDimension()), 
T_(T), 
r_(r),
dt(ds), 
sg(fl->getVols()),
log_S0(n_), 
mu(n_),
R(fl->getCorrelationMatrix().rankReducedRoot(r))
{   
	for(int i=0;i<n_;i++) log_S0[i]=log(S0[i]);
	for(int i=0;i<n_;i++) mu[i]=-sg[i]*sg[i]*dt/2;
}


BasketLattice* 
BasketLattice::
sample(int r, int n, int T)
{
	// initial asset prices
	RealArray1D S0(n);
	for(int j=0;j<n;j++) S0[j]=100.0;
			
	// volatilities
	RealArray1D sg(n);
	for(int j=0;j<n;j++) sg[j]=0.3;
			
	// correlation of returns
	UTRRealMatrix rho(n);
    for(int j=0;j<n;j++)
	for(int k=j;k<n;k++) rho(j,k)=exp(0.1*(j-k));
	
	// the factorloading
	ConstantFactorLoading* fl = new ConstantFactorLoading(n,sg,rho);
			
	// time step
	Real dt=0.1;
		
	return new BasketLattice(r,fl,T,dt,S0);
}



ostream& 
BasketLattice::
printSelf(ostream& os) const
{
	return
	os << "Basket lattice: " << r_ << "factors.";
}




void 
BasketLattice::
test(int r, int n, int T)
{
	BasketLattice* lattice = sample(r,n,T);
	lattice->selfTest();
	delete lattice;
}
		

void 
BasketLattice::
buildLattice(int m, bool verbose)
{
	std::cout << "\n\nBuilding lattice: " << *this
	          << "\nTime steps: " << m << endl << endl;
	switch(r_){
		
	    case 3 :  LatticeBuilder::buildThreeFactorLattice(this,m,verbose); 
			      break;
		default : LatticeBuilder::buildTwoFactorLattice(this,m,verbose); 
	}      	

} // end buildLattice
	


	

	
MTGL_END_NAMESPACE(Martingale)

 
