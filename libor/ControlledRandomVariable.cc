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


#include "ControlledRandomVariable.h"
using namespace Martingale;



		 
Real 
ControlledRandomVariable::
ControlledVariable::
nextValue() 
{
               
	Real mean_y=XC->getControlVariateMean(),     // control variate mean E(Y)
		 beta=XC->getBeta();                     // beta coefficient
    RealVector v=XC->nextValue();
    Real  x=v[0], y=v[1], xc=x-beta*(y-mean_y);
   	return xc;
} //end nextValue

			

RandomVariable* 
ControlledRandomVariable::
ControlledVariable::
controlled(){ return new ControlledVariable(this); }
		



Real 
ControlledRandomVariable::
expectation(int N){ return controlled()->expectation(N); }
    


Real 
ControlledRandomVariable::
expectation(int N, string message){ return controlled()->expectation(N,message); }
    
    


void 
ControlledRandomVariable::
controlVariateMeanTest(int N)
{
     cout << "\nTesting control variate mean:\n";
      // analytic control variate mean
     cout << "analytic: " << getControlVariateMean() << endl;
      // Monte Carlo control variate mean
     cout << "Monte Carlo: " << RandomVector::expectation(N)[1] << endl;
         
} // end controlVariateMeanTest
     
     

    
Real 
ControlledRandomVariable::
betaCoefficient()
{
     // Recall that Cov(X,Y)=E(XY)-E(X)E(Y) and Var(X)=E(X^2)-E(X)^2.
	int N=nBeta;
	Real sum_X=0, sum_Y=0,
	     sum_XX=0, sum_XY=0;
        
    for(int n=0;n<N;n++){
			
       RealVector v=nextValue();
       Real x=v[0],y=v[1];
            
       sum_X+=x; sum_Y+=y;
       sum_XX+=x*x; sum_XY+=x*y;
    }
        
    return (N*sum_XY-sum_X*sum_Y)/(N*sum_XX-sum_X*sum_X);
        
} //end betaCoefficient
    
    
     

Real 
ControlledRandomVariable::
correlationWithControlVariate(int N){ return correlation(0,1,N); } 
        
 


 