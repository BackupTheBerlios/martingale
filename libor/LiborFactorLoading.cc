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



#include "LiborFactorLoading.h"
#include "VolatilityAndCorrelation.h"
#include "Random.h"
#include "Utils.h"
#include <string>
#include <cmath>
#include <iostream>

using std::ostream;
using std::string;

MTGL_BEGIN_NAMESPACE(Martingale)



/*******************************************************************************
 *
 *                     LiborFactorLoadingType
 * 
 ******************************************************************************/

LiborFactorLoadingType::
LiborFactorLoadingType(int volSurface, int correlations) :
volType(volSurface), corrType(correlations) 
{    }                                                                         


string 
LiborFactorLoadingType::
volSurfaceType() const { return VolSurface::volSurfaceType(volType); }


string 
LiborFactorLoadingType::
correlationType() const { return Correlations::correlationType(corrType); }


ostream& 
LiborFactorLoadingType::
printSelf(ostream& os) const
{
	return
	os << "\nLibor factor loading: VolSurface: " << volSurfaceType() 
	   << ", Correlations: " << correlationType();
}


/*******************************************************************************
 *
 *                     LiborFactorLoading
 * 
 ******************************************************************************/


const UTRRealMatrix& 
LiborFactorLoading::
getRho() const { return corr->getCorrelationMatrix(); }

 
LiborFactorLoading::
LiborFactorLoading
(const RealArray1D& L0, const RealArray1D& deltas, const RealArray1D& _k,
VolSurface* vols, Correlations* corrs) :
flType(vols->getType(),corrs->getType()),
n(L0.getDimension()), 
delta(deltas), T(n+1), l(L0), x(n), k(_k),
vol(vols), corr(corrs), rho(*corrs)
{   
    // tenor structure, initial XLibors
	T[0]=0;
	for(int i=0;i<n;i++){ 
			
		T[i+1]=T[i]+delta[i];
		x[i]=delta[i]*l[i];
	}
		
	if(n!=corrs->getDimension()){
	    std::cout << "\n\nLiborfactorLoading(): Libors-Correlations dimensions don't match."
		          << "Libors n = " << n << ", correlations n = " << corrs->getDimension()
	              << "\nTerminating.";
	    exit(1);
	}

} // end constructor
	
	
	
LiborFactorLoading* 
LiborFactorLoading::
sample(int n, int volType=VolSurface::JR, int corrType=Correlations::CS)
{
	RealArray1D deltas(n);
	RealArray1D c(n);        // k_j
	RealArray1D L0(n);
    for(int i=0;i<n;i++)
	{ deltas[i]=0.25; c[i]=0.2+0.1*Random::U01(); L0[i]=0.04; }				
	VolSurface* vols = VolSurface::sample(volType);
	Correlations* corrs = Correlations::sample(n,corrType);
			
	return new LiborFactorLoading(L0,deltas,c,vols,corrs);
}
	
	
	
ostream& 
LiborFactorLoading::
printSelf(ostream& os) const
{
	RealVector L0(n,l.getData());      // initial Libors		
	RealVector vols(n-1,1);            // annualized volatilities
	for(int i=1;i<n;i++) vols[i]=annualVol(i);
	
	return
	os << "\n\nLiborFactorLoading: dimension = " << n 
	   << "\nInitial Libors: " << L0
	   << "\nAnnualized volatilities: " << vols
	   << *vol << *corr;
}
	
	
Real 
LiborFactorLoading::	
sigma(int i, Real t) const 
{ 
	return k[i]*(vol->sigma(t,T[i])); 
}
	
	
Real 
LiborFactorLoading::
annualVol(int i) const
{
	if(i==0) return 0.0;
	Real   T_i=T[i],
	       volSqr=integral_sgi_sgj_rhoij(i,i,0.0,T_i);
	return sqrt(volSqr/T_i);
}
					

// LOG-COVARIATION MATRICES 
	
	
Real 
LiborFactorLoading::
integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const
{  
   return k[i]*k[j]*rho(i,j)*(vol->integral_sgsg(s,t,T[i],T[j])); 
}

   	
	
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrix(int p,int q, Real s, Real t) const
{
    int size=q-p;       // matrix size
	UTRRealMatrix& logCVM=*(new UTRRealMatrix(size,p));

    for(int i=p;i<q;i++)
    for(int j=i;j<q;j++)
       logCVM(i,j)=integral_sgi_sgj_rhoij(i,j,s,t);

    return logCVM;
}// end logCovariationMatrix
   
   
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrix(int t) const
{
   return logLiborCovariationMatrix(t+1,n,T[t],T[t+1]);
}
   
   
const UTRRealMatrix& 
LiborFactorLoading::
logLiborCovariationMatrixRoot(int t) const
{
   return logLiborCovariationMatrix(t).utrRoot();
}


const RealMatrix& 
LiborFactorLoading::
reducedRankLogLiborCovariationMatrixRoot(int t, int r) const
{
   return logLiborCovariationMatrix(t).rankReducedRoot(r);
}

   
// TEST PROGRAMS 
		 
void 
LiborFactorLoading::
selfTest() const
{
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
    std::cout << "\nLIBOR FACTORLOADING SELFTEST:" << endl << *this;
	//printSelf(std::cout);
	
	cout << "\nTesting the root L of the matrix C=logLiborCovariationMatrix(t):" << endl;
	for(int t=0;t<n-1;t++){
		
		const UTRRealMatrix& C=logLiborCovariationMatrix(t);
		const UTRRealMatrix& L=logLiborCovariationMatrixRoot(t);
		Matrix<Real,UpperTriangular<Real> >& LLt=L.aat();
		LLt.setRowIndexBase(t+1); LLt.setColIndexBase(t+1);   // make equal to C
		cout << "\nt = " << t << ": ";
		C.testEquals(LLt,precision,epsilon,"C=L*L'");
	}

} // end test


void 
LiborFactorLoading::
factorizationTest(int r) const
{
	cerr << "\n\nApproximate factorization of all single time step log-Libor"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank " << r
	     << "\n\nRelative errors (trace norm): " << endl;
	
	for(int t=0;t<n-2;t++){
		
		const UTRRealMatrix& Ct=logLiborCovariationMatrix(t);
		Ct.testFactorization(r);
	}
}


// CALIBRATION

// use the first coordinates of x to set the Libor volatility scaling factors 
// k[j] (including k[0] although its useless) then the next coordinates to set 
// VolSurface::a,b,c,d and the follwing coordinates to set 
// Correlations::alpha,beta,r_oo
void 
LiborFactorLoading::
setParameters(const RealArray1D& u)
{
	vol->setParameters(u[0],u[1],u[2],u[3]);
	corr->setParameters(u[4],u[5],u[6]);
}
	

// Global Insertion

std::ostream& operator << 
(std::ostream& os, const LiborFactorLoading& fl){ return fl.printSelf(os); }



MTGL_END_NAMESPACE(Martingale)





