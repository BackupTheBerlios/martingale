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

#ifndef martingale_lmmlattice_h    
#define martingale_lmmlattice_h

// template class LmmLattice must be fully defined in header
#include "TypedefsMacros.h"
#include "Lattice.h"              // base class
#include "Node.h"                 // base class
#include "Utils.h"
#include "Array.h"                // direct member
#include "Matrix.h"               // direct member
#include "LiborFactorLoading.h" 
#include "Node.h"                 // problem with typedefed forward declarations


class std::ostream;



MTGL_BEGIN_NAMESPACE(Martingale)


// forward declaration
class BondCall;


/*! \file LmmLattice.h
 *  Lattices for the driftless Libor market model {@link DriftlessLMM}
 *  See book, chapter 8. It is assumed that the
 *  volatilities \f$\sigma_j(t)\f$ of the forward transported Libors \f$U_j\f$
 *  (see book, 6.8) are <i>constant</i> and that all Libor accrual intervals 
 *  have equal length \f$\delta\f$. Recall that the volatilities of the Libors
 *  themselves are stochastic in this model. The lattice makes time steps of
 *  equal size \f$\delta/nSteps\f$ where nSteps denotes the number of time steps 
 *  in each Libor accrual interval.
 *
 *  <p>The lattice uses {@link StandardBrownianNodes} since all functionals can
 *  be computed from the state of the driving Brownian motion 
 *  \f[Z=(Z_1,Z_2,\dots,Z_r)\f]
 *  at each node. Here r is th number of factors. Only r=2,3 are possible since 
 *  number of nodes <a href="#num-nodes">explodes</a> with the number of factors 
 *  (and time steps).
 *
 * <p><a name="lmm-lattice-3f">Three factor LMM lattice:</a>
 *  Nodes in a 3 factor lattice for the Libor market model {@link LmmLattice3F}
 *  compute the \f$H_j\f$ from the volatility parts \f$V_j\f$ of the forward transported
 *  Libors \f$U_j\f$ approximated as
 *  \f[V_j(s)\simeq\sigma_j\left[R_{j1}Z_1(s)+R_{j2}Z_2(s)+R_{j3}Z_3(s)\right],\f]
 *  where the matrix R is the approximate root of rank 3 of the correlation matrix 
 *  \f$\rho\f$ of the underlying driftless LMM: \f$\rho\simeq RR'\f$. 
 *  <p><a name="rescale"></a>The rows of R have to be
 *  scaled back to unit norm to preserve the correlations \f$\rho_{jj}=1\f$.
 *  This diminishes the quality of the approximation \f$\rho\simeq RR'\f$ but
 *  preserves the volatilities \f$\sigma_j\f$ of the \f$V_j\f$. Otherwise volatility
 *  is lost. See book, 8.1.2, for details and notation. Wether or not rescaling
 *  is indicated has to be determined by experiment. For example when pricing at the money
 *  swaptions rescaling leads to a significant deterioration in accuracy.
 *
 * <p>The \f$Z_j(t)\f$ are independent standard Brownian motions which evolve from the state
 *  \f$Z_j(0)=0\f$ and then tick up or down in ticks of size \f$a=\sqrt{dt}\f$ where dt is the 
 *  size of the time step. The state at any node is then given by the triple of integers
 *  (k0,k1,k2) with
 *  \f[Z_0=k0*a,\quad Z_1=k1*a,\quad Z_2=k2*a.\f]
 *  Two factor lattices are completely similar.
 *
 * <p><b>Arbitrage</b>
 * Recall that the \f$V_j\f$ are the volatility parts (unbounded variation parts) of the 
 * logarithms \f$Y_j=log(U_j)\f$. The \f$Y_j\f$ are then recovered from the \f$U_j\f$ using
 * the <i>continuous time</i> drifts of the \f$Y_j\f$. This does not preserve the martingale 
 * property of the \f$U_j\f$ modelled in the lattice and consequently the lattice is not 
 * arbitrage free. Our view here is that the lattice is an approximation of the arbitrage
 * free continuous time dynamics.
 *
 * <p><a name="num-nodes"><b>Number of nodes.</b></a>
 * Each node has four edges in the case of a two factor lattice and eight edges in the case of
 * a three factor lattice. The total number of nodes allocated depends on the number of time steps
 * in the lattice:
 * <ul>
 *     <li>A two factor lattice has \f$(t+1)^2\f$ nodes at time step t and 
 *         \f$(t+1)(t+2)(2t+3)/6\f$ nodes up to time step t (inclusively).
 *     </li>
 *     <li>A three factor lattice has \f$(t+1)^3\f$ nodes at time step t and
 *         \f$[(t+1)(t+2)/2]^2\f$ nodes up to time step t (inclusively).
 *     </li>
 * </ul>
 * We want to have a large number of nodes but not so large as to exceed main memory. 
 * With 1GB main memory we can tolerate about 11 million nodes corresponding to about 310
 * time steps in a two factor lattice and 6 million 
 * nodes corresponding to about 70 time steps in a three factor model. The number of time 
 * steps in the lattice can be controlled by choosing the number of time steps in each
 * Libor accrual interval. 
 */




/**********************************************************************************
 *
 *            GENERAL LMM LATTICE
 *
 *********************************************************************************/
 
 
	
 

/** <p>Lattice of the driftless Libor market model 
 *  with constant volatility surface and constant Libor accrual periods.
 *  See book 6.8, 8.1.1 for the details and notation. 
 *  For more details see the file reference for the file LmmLattice.h.
 *  Does not build a lattice, only the subclasses {@link LmmLattice2F} 
 *  {@link LmmLattice3F} do that.
 */
class LmmLattice : public Lattice<StandardBrownianNode> {
	
protected:
	
	// factor loading of the underlying driftless LMM
	LiborFactorLoading* factorLoading; 
	
	int n,                  // dimension (number of accrual intervals)
	    r,                  // number of factors
	    nSteps;             // number of time steps in each accrual interval
	
	/** constant Libor accrual period length.*/
	Real delta;
	/** time step.*/
	Real dt;
	/** ticksize a=sqrt(dt) of a standard Brownian motion over a single time step.*/
	Real a;
	/** the volatility scaling factors c_j, book, 6.11.1.*/
	RealArray1D sg;
	/** the initial values \f$log(U_j(0))\f$, book 6.8.*/
	RealArray1D log_U0;
    /** drift -sigma_j^2*dt/2 of Y_j=log(U_j) over a single time step.*/
	RealArray1D mu;  
	/** rank r approximate pseudo square root of the log(U_j) covariance matrix.*/
	RealMatrix R;     

	
public:


/** The number of factors. */
int getNumberOfFactors(){ return r; }

/** The factor loading of the underlying LMM.*/
LiborFactorLoading* getFactorLoading(){ return factorLoading; } 
	
/** The number of time steps per Libor compounding period. */
int getSteps(){ return nSteps; }
	
/** The size of the time step.*/
Real getTimeStep(){ return delta/nSteps; }

/** Returns largest t such that \f$T_t\leq s*dt\f$.
 * @param s number of time step.
 */
int getTenor(int s){ return (int)(s*dt/delta); }

/** Time and state independent transition probability along edge i.
 *  Declared final (hence nonvirtual) for speed. */

Real transitionProbability(int i){ return 1.0/(1<<r); }

		
// CONSTRUCTOR
	
/** 
 *  @param q number of factors: must be 2 or 3.
 *  @param fl factor loading of the underlying LMM, must have {@link CONST_VolSurface}.
 *  @param t lattice is built until Libor reset time T_t.
 *  @param steps number of equal sized time steps in each Libor accrual interval.
 *  @param verbose messages during build.
 */
LmmLattice(int q, LiborFactorLoading* fl, int t, int steps=1, bool verbose=false);

virtual ~LmmLattice(){ }


// SAMPLE LATTICE AND TEST

/** A sample r=2,3 factor lattice in dimension n (number of Libor accrual periods)
 *  built up to time T_p with nSteps time steps per accrual period.
 *  @param verbose details on lattice during build.
 */
static LmmLattice* sample(int r, int n, int p, int nSteps, bool verbose=false);


/** <a href="#rescale">Rescales</a> the rows of the rank reduced pseudo 
 *  square root R of the log(U_j)-correlation matrix C to unity (book, 8.1.2) 
 *  This preserves the Libor volatilities but makes the approximation 
 *  \f$C\simeq RR'\f$ less accurate. For example in the case of swaptions 
 *  the rescaled matrix will overestimate prices while prices will be 
 *  underestimated without rescaling. Rescaling is irreversible.
 */
void rescaleVols();

/** Tests the accuracy of the rank r factorization of the correlation matrix.*/
void testFactorization() const;

std::ostream& printSelf(std::ostream& os) const;


/** Builds an r=2,3 factor lattice in in dimension n (number of Libor accrual periods)
 *  up tp time \f$T_{n-3}\f$ with 3 time steps per Libor accrual period) and runs the 
 *  selfTest().
 */
static void test(int r, int n);


// LIBOR PATH FUNCTIONALS

/** <p>The vector \f$H=(H_p,\dots,H_n)\f$ of accrual factors at the node.
 *  Natural indices j=p,...,n.
 *  <p>This is a view of this vector in a static workspace (speed, no memory allocation).
 *  If the vector has to be preserved it must be copied before commands execute
 *  which overwrite the workspace.
 *
 * @param s time step at which the node lives.
 */
// The parameters is only needed to recognize the case s=0.
// It then propagates through all functions using this one.
// Keep it anyway, it might be useful for more complicated option payoffs.
const RealArray1D& Hvect(int p, StandardBrownianNode* node, int s);

/** Libor \f$L_j\f$ at the node.
 * @param s time step at which the node lives.*/
Real L(int j, StandardBrownianNode* node, int s);


/** The forward price \f$H_{p,q}\f$ of the annuity \f$B_{p,q}\f$ over the 
 *  interval [T_p,T_q] at this node.
 * @param s time step at which the node lives.*/
Real H_pq(int p, int q, StandardBrownianNode* node, int s);


/** The swaprate \f$S_{p,q}\f$ for a swap on [T_p,T_q] at the node.
 * @param s time step at which the node lives.*/
Real swapRate(int p, int q, StandardBrownianNode* node, int s);
		

/** Payoff of a forward swaption with strike rate kappa exercising into 
 *  a swap on the interval [T_p,T_q] if exercised at this node. 
 *  Payoff is compounded forward to the horizon T_n. 
 */
Real forwardSwaptionPayoff(int p, int q, Real kappa, StandardBrownianNode* node, int s);


/** Forward compounded payoff of caplet(i) with strike rate kappa at this node. 
 *  Assumes that the node lives at the Libor reset point \f$T_i\f$ at which Libor 
 *  for this caplet is set.
 *  @param s time step at which the node lives.*/
Real forwardCapletPayoff(int i, Real kappa, StandardBrownianNode* node, int s);


/** Forward payoff of the BondCall bc if exercised at this node. 
 * @param s time step at which the node lives.*/
Real forwardBondCallPayoff(BondCall* bc, StandardBrownianNode* node, int s);

			
private:
	
// build lattice with m time steps.
void buildLattice(int m, bool verbose);

// continuous time reached after time step s=0,1,...
Real tau(int s);

	static RealArray1D H_;     // workspace for accrual factors H_j
	static RealArray1D V_;     // workspace for volatility parts V_j of the log(U_j)
	                           // book, 8.1	

		
}; // end LmmLattice





MTGL_END_NAMESPACE(Martingale)

#endif
 
