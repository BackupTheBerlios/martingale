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

#include "TestOptimizers.h"
#include "Optimizer.h"
#include "Utils.h"

using std::cout;


MTGL_BEGIN_NAMESPACE(Martingale)


void 
testDownhillSimplex(int n, int steps)
{
     printStars();
	 cout << "\n\nTesting Downhill Simplex.";
	 RealArray1D x(n);
     for(int i=0;i<n;i++)x[i]=1.6+(i%2)*0.6;
         
     Real delta=0.3;
     bool verbose=true;
			 			 
     Optimizer* optimizer=new ConcreteDownhillSimplex
     (&(ObjectiveFunction::function_1),x,delta,steps,verbose);         
     optimizer->search();
		 
} // end test          



void 
testBFGS(int n, int nVals)
{
     printStars();
	 cout << "\n\nTesting BFGS.";
	 RealArray1D x(n), h(n);
     for(int i=0;i<n;i++){ x[i]=1.6+(i%2)*0.6; h[i]=0.1; }
         
     Real stepmax=0.3;
	 int nRestarts=5;
     bool verbose=true;
			 			 
     Optimizer* optimizer=new ConcreteBFGS
	 (&(ObjectiveFunction::function_1),x,nVals,stepmax,h,nRestarts,verbose);          
     optimizer->search();
		 
} // end test  
	 
	

void 
testSobolSearch(int n, int nVals)
{
     printStars();
	 cout << "\n\nTesting Sobol search.";
     RealArray1D x(n), delta(n);
     for(int i=0;i<n;i++){ x[i]=1.6+(i%2)*0.6; delta[i]=0.5; }
         
     bool verbose=true;
			 			 
     Optimizer* optimizer=new ConcreteSobolSearch
	 (&(ObjectiveFunction::function_1),x,nVals,delta,verbose);         
     optimizer->search();
		 
} // end test          


	     	
MTGL_END_NAMESPACE(Martingale)
