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


#include "DirichletProblem.h"
#include "StochasticProcesses.h"
#include <iostream>
using namespace Martingale; 



DirichletProblem::	
DirichletProblem(EuclideanRegion* D, int T, Real dt=0.01) : 
dim(D->getDimension()),
G(D), 
X(new VectorBrownianMotion(dim,T,dt)),
tau(new FirstExitTime<RealVector,Real>(X,G))
{   }
	
	
Real 
DirichletProblem::
solution(const RealVector& x, int nPath=50000, bool reportHits=false)
{		
	X->setStart(x);
	int c=0;                                 // hits on the boundary
	Real sum=0;
	for(int j=0;j<nPath;j++){
			
		int t=X->newPathSegment(tau);        // first time outside G
		int s=t-1;
		const RealVector& u=X->currentPath(s);
		const RealVector& v=X->currentPath(t);
		if(!(G->isMember(v))) c++;            // hit on the boundary            
		// point where the boundary is hit
        const RealVector& z=G->boundaryIntersection(u,v); 
		sum+=boundaryFunction(z);
	}
		
	if(reportHits) 
	cout << "\nDirichletProblem::solution: " << 100.0*c/nPath 
		 << "% of paths hit the boundary.";
	return sum/nPath;
} // end solution
	



DirichletProblemExample::
DirichletProblemExample(int dim, int T, Real dt=0.01) : 
DirichletProblem(new Ball(dim),T,dt) 
{  }
	

Real 
DirichletProblemExample::
boundaryFunction(const RealVector& u)
{ Real f=u[0]; for(int i=1;i<dim;i++) f+=u[i]; return f; }
		
	
void 
DirichletProblemExample::	
runExample(int d, int T)
{
	cout << "\n\nSolution f(x) of the Dirichlet problem on the unit ball"
	     << "\nfor the boundary function h(x_1,x_2,...,x_d)=x_1+x_2+...+x_d"
	     << "\nat the point x_1=x_2=...=x_d=0.75."
	     << "\nDimension  d = " << d;
		
	DirichletProblemExample Example(d,T);
	RealVector x(d);
	for(int i=0;i<d;i++) x[i]=1.0/(2*sqrt(d));
			
	cerr << "\n\nAnalytic: f(x) = " << Example.boundaryFunction(x);
	cerr << "\nMonte Carlo computation...";
	Real fx=Example.solution(x,30000,true);
	cerr << "\nMonte Carlo: f(x) = " << fx 
	     << "\nQuasi Monte Carlo computation...";
	Example.X->switchToQMC();
	fx=Example.solution(x,30000,true);
	cerr << "\nQuasi Monte Carlo: f(x) = " << fx;
} 


