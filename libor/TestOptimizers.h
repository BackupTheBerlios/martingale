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

#ifndef martingale_testoptimizers_h    
#define martingale_testoptimizers_h

#include "Optimizer.h"


MTGL_BEGIN_NAMESPACE(Martingale)

/** A collection of free standing short test programs.
 */
 

/*******************************************************************************    
    
              TEST OF DownHillSimplex
	
*******************************************************************************/

	 /** Test program in dimension n. Very nasty test function 
	  *  \f[f(x)=1\bigg/\sum\nolimits_jexp(-x_j^2)\f]
      *  consists of three very narrow valleys intersecting at right angles. 
	  *  The optimizer must crawl along each valley minimizing each variable separately. 
	  *  Takes a huge number of steps. Minimum is assumed at the origin.
	  *
	  * @param n dimension.
	  * @param steps maximum number steps.
      */
	 void testDownhillSimplex(int n, int steps)
     {
         Array1D<Real> x(n);
         for(int i=0;i<n;i++)x[i]=1.6+i*0.6;
         
         Real delta=0.3;
         bool verbose=true;
			 			 
         Optimizer* optimizer=new ConcreteDownhillSimplex
		 (&(ObjectiveFunction::function_1),n,x,delta,steps,verbose);         
         optimizer->search();
		 
     } // end test          


/*******************************************************************************    
    
              TEST OF BFGS
	
*******************************************************************************/

	 /** <p>Test program in dimension n. Very nasty test function 
	  *  \f[f(x)=1\bigg/\sum\nolimits_jexp(-x_j^2)\f]
      *  consists of n very narrow valleys intersecting at right angles. 
	  *  The optimizer must crawl along each valley minimizing each variable separately. 
	  *  Takes a huge number of steps. Minimum is assumed at the origin. The search starts 
	  *  the point \f$x_j=1.6+j*0.6\f$.
	  *
	  * <p>For \f$n\geq 3\f$ this fails miserably. It seems that we have to rescale the 
	  * objective function so that the gradient has coordinates of roughly the same
	  * magnitude. 
	  *
	  * @param n dimension.
	  * @param nVals maximum number of function evaluations.
      */
	 void testBFGS(int n, int nVals)
     {
         RealArray1D x(n), h(n);
         for(int i=0;i<n;i++){ x[i]=1.6+i*0.6; h[i]=0.1; }
         
         Real stepmax=0.3;
		 int nRestarts=5;
         bool verbose=true;
			 			 
         Optimizer* optimizer=new ConcreteBFGS
		 (&(ObjectiveFunction::function_1),n,x,nVals,stepmax,h,nRestarts,verbose);          
         optimizer->search();
		 
     } // end test     	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 