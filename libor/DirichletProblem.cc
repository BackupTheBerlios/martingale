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
#include "Utils.h"
#include <iostream>
#include <cmath>

using std::floor;

MTGL_BEGIN_NAMESPACE(Martingale)



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
	std::cout << "\nDirichletProblem::solution: " << 100.0*c/nPath 
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
runExample()
{
    int do_again=1;
    // main loop
    while(do_again==1){	
        
        printStars();
	cout << "\n\nSolution f(x) of the Dirichlet problem on the unit ball"
	     << "\nfor the boundary function h(x_1,x_2,...,x_d)=x_1+x_2+...+x_d"
	     << "\nat the point x_1=x_2=...=x_d=1/(2*sqrt(d))."
	     << "\n\nEnter dimension  d = ";
	
        int d; cin >> d;
        int T =(int) floor(10.0*(d+20)/d);	
	DirichletProblemExample Example(d,T);
	RealVector x(d);
	for(int i=0;i<d;i++) x[i]=1.0/(2*sqrt(d));
			
	cout << "\n\nAnalytic: f(x) = " << Example.boundaryFunction(x);
	cout << "\nMonte Carlo computation...";
	Real fx=Example.solution(x,100000,true);
	cout << "\nMonte Carlo: f(x) = " << fx 
	     << "\nQuasi Monte Carlo computation...";
	Example.X->switchToQMC();
	fx=Example.solution(x,30000,true);
	cout << "\nQuasi Monte Carlo: f(x) = " << fx;

	cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
	cin >> do_again;
    }
} 




MTGL_END_NAMESPACE(Martingale)
