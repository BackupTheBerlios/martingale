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

#include "TestProbability.h"
#include "Random.h"
#include "RandomVariables.h"
#include "StochasticProcesses.h"
#include "QuasiMonteCarlo.h"


MTGL_BEGIN_NAMESPACE(Martingale)


void 
testSTNCovarianceMatrix(int d, int N)
{
    StandardNormalVector Z(d);
	Timer watch;
    watch.start();
	UTRMatrix<Real> cvm=Z.covarianceMatrix(N,"Standard normal covariance matrix");
	
	// get the largest and smallest diagonal element and the largest and smallest 
	// off diagonal element
	Real mindiag=100.0, maxdiag=0.0,
	     minoffdiag=100000.0, maxoffdiag=-100000.0;
	for(int i=0;i<d;i++)
	for(int j=i;j<d;j++){
		
		Real cvm_ij=cvm(i,j);
		if(i==j){
			if(cvm_ij<mindiag)mindiag=cvm_ij;
			if(cvm_ij>maxdiag)maxdiag=cvm_ij;
		}
		
		if(i!=j){
			if(cvm_ij<minoffdiag)minoffdiag=cvm_ij;
			if(cvm_ij>maxoffdiag)maxoffdiag=cvm_ij;
		}
	} // end for i
			
	watch.stop();	
    std::cout << endl << "Standard normal covariance matrix, dimension: " << d
	          << endl << "Maximal diagonal element: " << maxdiag
	          << endl << "Minimal diagonal element: " << mindiag
	          << endl << "Maximal offdiagonal element: " << maxoffdiag
	          << endl << "Minimal diagonal element: " << minoffdiag;
    watch.report("Finished");
}



void 
testPathFunctional(int t, int T, int nPath)
{ VectorBrownianMotion::testPathFunctional(t,T,nPath); }
	
	     
	
MTGL_END_NAMESPACE(Martingale)
