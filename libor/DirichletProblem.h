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

#ifndef martingale_dirichletproblem_h
#define martingale_dirichletproblem_h

#include "TypedefsMacros.h"
#include "Matrix.h"                           // typedef problem in forward declarations

MTGL_BEGIN_NAMESPACE(Martingale) 


// we are using
class DirichletProblemExample;
class EuclideanRegion;
class Ball;
class VectorBrownianMotion;
class StoppingTime;
// class RealVector


/*! \file DirichletProblem.h
 *  Solving the Dirichlet problem for a EuclideanRegion G (in any dimension): 
 *  for a given function h defined on the boundary of G we seek a function f satisfying 
 *  \f$\nabla f=0\f$ on the interior, \f$f=h\f$ on the boundary of G.
 *  Solution computed at arbitrary points x in G using Brownian motion.
 *  For more details see {@link DirichletProblem}.
 */
 
 
/** <p>This class solves the Dirichlet problem for a EuclideanRegion G (in any dimension): 
 *  for a given function h defined on the boundary of G we seek a function f satisfying 
 *  \f[\nabla f=0,\hbox{ on the interior, }f=h\hbox{ on the boundary of }B.\f]
 *
 *  <p>Brownian motion can be used to compute f(x) as follows: let X be a Brownian 
 *  motion starting at the point x in G. If the function f satisfies \f$\nabla f=0\f$
 *  then \f$f(X_t)\f$ is a martingale and hence \f$E[f(X_{\tau})]=f(X_0)=f(x)\f$,
 *  for each stopping time \f$\tau\f$. In particular if \f$\tau\f$ is the first exit 
 *  time from G (equivalently the hitting time for the boundary of G), then 
 *  \f$f(X_{\tau})=h(X_{\tau})\f$ since f=h on the boundary of G. Consequently we can 
 *  compute f(x) as
 *  \f[f(x)=E[h(X_{\tau})]\f]
 *  where \f$\tau\f$ is the first exit time from G.
 */
class DirichletProblem {
	
protected:
	
	int dim;                    // dimension
	EuclideanRegion* G;         // region on which it is solved
	VectorBrownianMotion* X;
	StoppingTime* tau;          // the first exit time of X from G

public: 
	
	/** The boundary function h. */
	virtual Real boundaryFunction(const RealVector& u) = 0;
	
// CONSTRUCTOR
	
	/** The number T of time steps is fixed in advance to enable a
     *  Quasi Monte Carlo computation where this number determines the dimension
     *  of the Sobol generator driving the Brownian motion. This number has to be 
     *  large enough to ensure that the boundary is hit. The solution routine can
	 *  report what percentage of paths hit the boundary. If this is too small
	 *  increase T.
	 *
	 *  @param D the region on which the problem is solved.
	 *  @param T time steps alloted to hit the boundary.
	 *  @param dt size of time steps.
	 */
	DirichletProblem(EuclideanRegion* D, int T, Real dt=0.01);
	
	
	/** The solution \f$f(x), x\in G\f$, 
	 *  each value computed from nPath paths of X.
	 *
	 * @param nPath number of paths launched for the boundary.
	 * @param reportHits report what percentage of paths hit the boundary.
	 */
	Real solution(const RealVector& x, int nPath=50000, bool reportHits=false);
	
	
}; // end DirichletProblem



// AN EXAMPLE


/** Dirichlet problem on the unit ball with boundary function
 *  \f$h(x_1,x_2,\ldots,x_d)=x_1+x_2+\dots+x_d\f$. This function is harmonic 
 *  everywhere so the solution is the same function on the interior of the ball.
 */
class DirichletProblemExample : public DirichletProblem {
	
public:
	
	/** The number of time steps necessary to reliably hit the boundary
	 *  depnds on the size dt of the time steps and the dimension.
	 *  For the default dt=0.01, T=500 suffices in dimension d=2,
	 *  while T=10 suffices in dimension d=30.
	 */
	DirichletProblemExample(int dim, int T, Real dt=0.01);
	
	Real boundaryFunction(const RealVector& u);
		
	
	/** <p>Computes f(x) for the boundary function 
	 *  \f$h(x_1,x_2,\ldots,x_d)=x_1+x_2+\dots+x_d\f$ on the unit ball
	 *  in \f$R^d\f$. This function is harmonic (\f$\nabla h=0\f$) and so
	 *  the solution on the interior has the same form.
	 *  We compute the solution f(x) for \f$x_1=x_2=\dots=x_n=1/2\sqrt{d}\f$
	 *  and then compare the known analytic value with what we get.
	 *  Size of time steps 0.01. Number of paths: 30000.
	 *
	 * @param d dimension.
	 * @param T time steps alloted to hit boundary.
	 */
	static void runExample(int d, int T);
	
	
}; // end DirichletProblemExample
		


MTGL_END_NAMESPACE(Martingale) 


#endif

