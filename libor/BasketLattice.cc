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
#include <iostream>
#include <cmath>


using std::ostream;
using std::cout;
using std::log;
using std::exp;


MTGL_BEGIN_NAMESPACE(Martingale)

 
// BASKET LATTICE DATA


BasketLatticeData::
BasketLatticeData
(int steps,
 Real dt,
 const RealArray1D& vols,
 const RealArray1D& Y0,
 const RealMatrix& Q
) :
n(Q.rows()), T(steps), r(Q.cols()), 
timestep(dt), ticksize(sqrt(dt)),
sg(vols), log_S0(Y0), driftUnit(n), R(Q)
{  
	// drifts per time step
	for(int i=0;i<n;i++) driftUnit[i]=-sg[i]*sg[i]*dt/2;
}
	

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
(int r, ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0) : 
Lattice<BasketNode>(T),
n_(fl->getDimension()), T_(T), 
dt_(dt), Y0_(n_), sg_(fl->getVols()),
R_(fl->correlationMatrix().rankReducedRoot(r)),
latticeData(0)
{   
	for(int i=0;i<n_;i++) Y0_[i]=log(S0[i]);
	latticeData = new BasketLatticeData(T,dt,sg_,Y0_,R_);
}


ostream& 
BasketLattice::
printSelf(ostream& os) const
{
	int r = latticeData->r;
	return
	os << "Basket lattice: " << r << "factors.";
}



// TWO FACTOR BASKET LATTICE

BasketLattice2F::
BasketLattice2F
(ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0) : 
BasketLattice(2,fl,T,dt,S0) 
{    
	buildLattice(T);
	testFactorization();
}


BasketLattice2F* 
BasketLattice2F::
sample(int n,int T)
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
		
	return new BasketLattice2F(fl,T,dt,S0);
}


void 
BasketLattice2F::
test(int n, int T)
{
	BasketLattice2F* lattice = sample(n,T);
	lattice->selfTest();
	delete lattice;
}
		

void 
BasketLattice2F::
buildLattice(int m)
{
	std::cout << "\n\nBuilding lattice: " << *this
	          << "\nTime steps: " << m << endl << endl;
	LatticeBuilder::
	buildTwoFactorLattice<BasketNode,BasketLatticeData>(m,nodeList,latticeData);         	

} // end buildLattice
	


// THREE FACTOR LATTICE


BasketLattice3F::
BasketLattice3F
(ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0) : 
BasketLattice(3,fl,T,dt,S0) 
{    
	buildLattice(T);
	testFactorization();
}


BasketLattice3F* 
BasketLattice3F::
sample(int n, int T)
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
	for(int k=j;k<n;k++) rho(j,k)=exp(0.2*(j-k));
	
	// the factorloading
	ConstantFactorLoading* fl = new ConstantFactorLoading(n,sg,rho);
			
	// time step
	Real dt=0.1;
		
	return new BasketLattice3F(fl,T,dt,S0);
}



void 
BasketLattice3F::
test(int n, int T)
{
	BasketLattice3F* lattice = sample(n,T);
	lattice->selfTest();
	delete lattice;
}



void 
BasketLattice3F::
buildLattice(int m)
{
	std::cout << "\n\nBuilding lattice: " << *this
	          << "\nTime steps: " << m << endl << endl;
	LatticeBuilder::
	buildThreeFactorLattice<BasketNode,BasketLatticeData>(m,nodeList,latticeData);         	

} // end buildLattice
	

	
MTGL_END_NAMESPACE(Martingale)

 
