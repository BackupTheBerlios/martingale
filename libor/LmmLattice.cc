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
#include <cmath>


MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::log;


// LMM LATTICE DATA


LmmLatticeData::
LmmLatticeData
(int steps,
 Real dt,
 const RealArray1D& vols,
 const RealArray1D& Y_0, 
 const RealMatrix& Q
) :
n(1+Q.rows()), nSteps(steps), r(Q.cols()), 
timestep(dt), ticksize(sqrt(dt)), delta(steps*dt),
sg(vols), log_U0(Y_0), driftUnit(Q.rows(),1), R(Q)
{  
	for(int i=1;i<n;i++) driftUnit[i]=-sg[i]*sg[i]*dt/2; 
}
	


// GENERAL LMM LATTICE

LmmLattice::
LmmLattice(int q, LiborFactorLoading* fl, int t, int steps=1) : 
Lattice<LmmNode>(t*steps),
n(fl->getDimension()), 
r(q),
nSteps(steps),
delta(fl->getDeltas()[0]),
dt(delta/steps),
a(sqrt(dt)),
Y0(n),
factorLoading(fl),
sg(n-1,1),
R(fl->getRho().rankReducedRoot(r)),
latticeData(0)
{ 
	// check the number of factors
	if((r!=2)&&(r!=3)){
		
	   cout << "\n\nLmmLatticeconstructor: number of factors must be 2 or 3 but is " << r
	             << "\nTerminating.";
	   exit(1);
    }		
		
	// check if Libor accrual periods are constant
	const RealArray1D& deltas=fl->getDeltas();
	Real delta=deltas[0];
	for(int j=0;j<n;j++) if(deltas[j]!=delta) {
			
	   cout << "\n\nLmmLatticeconstructor: Libor accrual periods not constant."
	             << "\nTerminating.";
	   exit(1);
    }
	
	// check if volatilities are constant
	if(fl->getType().volType!=VolSurface::CONST) {
			
	   cout << "\n\nLmmLattice-constructor: volatility surface not constant."
	             << "\nTerminating.";
	   exit(1);
    }	
		
	// set log(U_j(0)), j=0,1,...,n-1
    const RealArray1D& x=factorLoading->getInitialXLibors();     // x[j]=X_j(0)
    for(int j=0;j<n;j++){ 
			
	    // U_j(0)=X_j(0)(1+X_{j+1}(0))...(1+X_{n-1}(0))
		Real Uj0=x[j]; for(int k=j+1;k<n;k++) Uj0*=1+x[k]; 
		Y0[j]=log(Uj0);
	}

    // set constant vols
	for(int j=1;j<n;j++) sg[j]=factorLoading->sigma(j,0.0);
			
	latticeData = new LmmLatticeData(steps,dt,sg,Y0,R);
	
} // end constructor



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
	os << "LMM lattice:"
	   << "\nNumber of time steps in each accrual interval: " << nSteps 
	   << endl << *factorLoading;
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
	


// TWO FACTOR LMM LATTICE


LmmLattice2F::
LmmLattice2F
(LiborFactorLoading* fl, int t, int steps=1,bool verbose=false) :
LmmLattice(2,fl,t,steps) 
{   
	buildLattice(t*steps,verbose);     // m=t*steps is the number of time steps
	if(verbose) testFactorization(); 
}



LmmLattice2F* 
LmmLattice2F::
sample(int n, int p, int nSteps, bool verbose=false) 
{
	LiborFactorLoading*
	fl=LiborFactorLoading::sample(n,VolSurface::CONST,Correlations::CS);
	LmmLattice2F* lattice = new LmmLattice2F(fl,p,nSteps,verbose);
	return lattice;
}


void 
LmmLattice2F::
test(int n)
{
    Timer watch; watch.start();
	LmmLattice2F* lattice = sample(n,n,1,true);
	lattice->rescaleVols();
	lattice->selfTest();
	delete lattice;
	watch.stop();
	watch.report("Two factor LmmLattice test");
}
		

void 
LmmLattice2F::
buildLattice(int m, bool verbose)
{
	if(verbose)
	cout << "\n\nBuilding lattice: " << *this
	          << "\n\nTime steps: " << m << endl << endl;
	LatticeBuilder::
	buildTwoFactorLattice<LmmNode,LmmLatticeData>(m,nodeList,latticeData,verbose);         		  

} // end buildLattice



// THREE FACTOR LMM LATTICE


LmmLattice3F::
LmmLattice3F
(LiborFactorLoading* fl, int t, int steps=1, bool verbose=false) :	
LmmLattice(3,fl,t,steps) 
{   
	buildLattice(t*steps,verbose);     // m=t*steps is the number of time steps
	if(verbose) testFactorization();
}


LmmLattice3F* 
LmmLattice3F::
sample(int n, int p, int nSteps, bool verbose=false) 
{
	LiborFactorLoading*
	fl=LiborFactorLoading::sample(n,VolSurface::CONST,Correlations::CS);
	LmmLattice3F* lattice = new LmmLattice3F(fl,p,nSteps,verbose);
	return lattice;
}


void 
LmmLattice3F::
test(int n)
{
	Timer watch; watch.start();
	LmmLattice3F* lattice = sample(n,n,1,true);
	lattice->rescaleVols();
	lattice->selfTest();
	delete lattice;
	watch.stop();
	watch.report("Three factor LmmLattice test");
}
		
		
void 
LmmLattice3F::
buildLattice(int m, bool verbose)
{
	if(verbose)
	cout << "\n\nBuilding lattice: " << *this
	          << "\n\nTime steps: " << m << endl << endl;
	LatticeBuilder::
	buildThreeFactorLattice<LmmNode,LmmLatticeData>(m,nodeList,latticeData,verbose);         	

} // end buildLattice



MTGL_END_NAMESPACE(Martingale)

