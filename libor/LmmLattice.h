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
#include "Lattice.h"                      
#include "Utils.h"
#include "Matrix.h"
#include "LiborFactorLoading.h" 
//#include <iostream>

class std::ostream;



MTGL_BEGIN_NAMESPACE(Martingale)


// forward declaration
class LmmNode;                  // Node.h



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
 * <p>Lattices depend on the number f of factors and the type of LmmNode (heavy or lite).
 *  Even though the lattice is fully recombining the number of nodes explodes with the 
 *  number of factors and time steps. Realistically only 3 factor latttices are possible
 *  on PC equipment.
 *
 * <p> Heavy nodes store the entire vector H of accrual factors in each node
 *  All other functionals are computed from the accrual factors. A lattice built with
 *  these nodes requires a lot of memory and builds itself slowly. Once built however 
 *  it prices all derivatives very fast. This version is to be used if several derivatives 
 *  are to be priced based on the same lattice. 
 *  
 * <p>Lite nodes store only the state variables from which the accrual factors 
 *  H_j are computed. A lattice built with these nodes uses much less memory and builds 
 *  more quickly. Such a lattice can use many more nodes. However pricing derivatives 
 *  is slower. This is to be used if many nodes are necessary or if only one derivative
 *  is to be priced.
 *
 * <p><a name="lmm-lattice-3f">Three factor LMM lattice:</a>
 *  Nodes in a 3 factor lattice for the Libor market model {@link LmmLattice3F}
 *  compute the \f$H_j\f$ from the volatility parts \f$V_j\f$ of the forward transported
 *  Libors \f$U_j\f$ approximated as
 *  \f[V_j(s)\simeq\sigma_j\left[R_{j1}Z_1(s)+R_{j2}Z_2(s)+R_{j3}Z_3(s)\right],\f]
 *  where the matrix R is the approximate root of rank 3 of the correlation matrix 
 *  \f$\rho\f$ of the underlying driftless LMM: \f$\rho\simeq RR'\f$. 
 *  <p><a name="rescale"></a>The rows of R
 *  are scaled back to unit norm to preserve the correlations \f$\rho_{jj}=1\f$.
 *  This diminishes the quality of the approximation \f$\rho\simeq RR'\f$ but
 *  preserves the volatilities \f$\sigma_j\f$ of the \f$V_j\f$. If we don't do this we 
 *  lose volatility. See book, 8.1.2, for details and notation. 
 *
 * <p>The \f$Z_j(t)\f$ are independent standard Brownian motions which evolve from the state
 *  \f$Z_j(0)=0\f$ and then tick up or down in ticks of size \f$a=\sqrt{dt}\f$ where dt is the 
 *  size of the time step. The state at any node is then given by the triple of integers
 *  (i,j,k) with
 *  \f[Z_1=ia,\quad Z_2=ja,\quad Z_3=ka.\f]
 *  Two factor lattices are completely similar but the implementation is kept separate for
 *  greater clarity.
 *
 * <p><b>Arbitrage</b>
 * Recall that the \f$V_j\f$ are the volatility parts (unbounded variation parts) of the 
 * logarithms \f$Y_j=log(U_j)\f$. The \f$Y_j\f$ are then recovered from the \f$U_j\f$ using
 * the <i>continuous time</i> drifts of the \f$Y_j\f$. This does not preserve the martingale 
 * property of the \f$U_j\f$ modelled in the lattice and consequently the lattice is not 
 * arbitrage free. Our view here is that the lattice is an approximation of the arbitrage
 * free continuous time dynamics.
 *
 * <p><b>Number of nodes.</b>
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
 * With 1GB main memory we can tolerate about 5.3 million lightweight nodes
 * corresponding to 250 time steps in a two factor lattice and 3.5 million lightweight 
 * nodes corresponding to 60 time steps in a three factor model. The number of time 
 * steps in the lattice can be controlled by choosing the number of time steps in each
 * Libor accrual interval. 
 */




/**********************************************************************************
 *
 *            GENERAL LMM LATTICE
 *
 *********************************************************************************/
 
 
/** LmmLattice data which must be conveyed to the nodes.
 *  All arrays are indexed with natural indices from the Libor process context.
 */
struct LmmLatticeData {
	
	int 
	    /** dimension of Libor process: */                          n,        
	    /** number of time steps in each Libor accrual interval:*/  nSteps,   
	    /** number of factors: */                                   r;
	
	/** time step.*/
	Real timestep;
	/** ticksize a=sqrt(dt) of a standard Brownian motion over a single time step.*/
	Real ticksize;
	/** constant Libor accrual period length.*/
	Real delta;
	/** the volatility scaling factors c_j, book, 6.11.1.*/
	const RealArray1D& sg;
	/** the initial values \f$log(U_j(0))\f$, book 6.8.*/
	const RealArray1D& log_U0;
    /** drift -sigma_j^2*dt/2 of Y_j=log(U_j) over a single time step.*/
	RealArray1D driftUnit;  
	/** rank r approximate pseudo square root of the log(U_j) covariance matrix.*/
	const RealMatrix& R;     
	

	/** @param steps number of time steps in each Libor accrual interval.
	 *  @param dt size of time step.
	 *  @param vols the constant volatilities \f$\sigma_j\f$, book, 6.11.1.
	 *  @param Y_0 the initial values \f$Y_j(0)=log(U_j(0))\f$, book 6.8.
	 *  @param drifts drifts of the log(U_j) over a single time step.
	 *  @param Q low factor pseudo square root of the log(U_j) covariation matrix.
	 */
	LmmLatticeData
	(int steps, Real dt, const RealArray1D& vols,
	 const RealArray1D& Y_0, const RealMatrix& Q);
	
}; // end LmmLatticeInfo

	
 

/** <p>Lattice of the driftless Libor market model with constant volatility 
 *  surface. This is a lattice with nodes of type 
 *  <center><code>
 *  LmmNode&lt;LmmNodeBase&gt;, 
 *  </code></center>
 *  where LmmNodeBase=LmmNode::BaseType.
 *  The base determines if the nodes are lightweight or heavyweight.
 *  See {@link Node}, {@link LmmNode_LiteBase} or {@link LmmNode_HeavyBase}. 
 *  See book 6.8, 8.1.1 for the details and notation. 
 *  For more details see the file reference for the file LmmLattice.h.
 */
class LmmLattice : public Lattice<LmmNode> {
	
protected:
	
	int n,                  // dimension (number of accrual intervals)
	    r,                  // number of factors
	    nSteps;             // number of time steps in each accrual interval
	
	Real delta;             // constant accrual period
	Real dt;                // size of time step
	Real a;                 // a=sqrt(dt) ticksize of Brownian motion over one time step
	
	RealArray1D Y0;         // Y0[j]=Y_j(0)=log(U_j(0))

	// factor loading of the underlying driftless LMM
	LiborFactorLoading* factorLoading; 
	
	RealArray1D sg;         // the constant volatilities sigma_j, book, 6.11.1.
	RealMatrix& R;          // rank r approximate pseudosquare root of the log(U_j) covariance matrix
	
	LmmLatticeData* latticeData;    // data object passed to nodes

	
public:

// ACCESSORS

    /** The factor loading of the underlying LMM.
	 */
    LiborFactorLoading* getFactorLoading(){ return factorLoading; } 

    /** Dimension of underlying LMM (number of accrual periods).
	 */
    const LmmLatticeData* getData(){ return latticeData; }
	
	
// CONSTRUCTOR
	
/** 
 *  @param q number of factors: must be 2 or 3.
 *  @param fl factor loading of the underlying LMM, must have {@link CONST_VolSurface}.
 *  @param t lattice is built until Libor reset time T_t.
 *  @param steps number of equal sized time steps in each Libor accrual interval.
 */
LmmLattice(int q, LiborFactorLoading* fl, int t, int steps=1);

~LmmLattice(){ delete &R, delete latticeData; }



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
	
private:


// continuous time reached after time step s=0,1,...
Real tau(int s);
	

		
}; // end LmmLattice



/**********************************************************************************
 *
 *            CONSTANT VOLATILITY TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
		

//<----------------To Do----------------->
// make this a class template depending on the Node type (lite/heavy).


/** <p><a href="lmm-lattice-3f">Two factor lattice</a> for driftless Libor market model  
 *  {@link DriftlessLMM} with constant volatility functions \f$\sigma_j(t)=\sigma_j\f$. 
 *
 * <p>It is assumed that all Libor accrual intervals have equal length \f$\delta\f$..
 *  The lattice makes nSteps time steps of equal length \f$dt=\delta/nSteps\f$ in 
 *  each accrual interval. Consequently discrete time t corresponds to continuous time
 *  \f[t*dt=T_{t/nSteps}+(t\;mod\;nSteps)*dt.\f]
 *  For more details see the file reference for the file LmmLattice.h.
 */
class LmmLattice2F : public LmmLattice {

public:

	
// CONSTRUCTOR
	
/** The {@link VolSurface} of the factor loading must be of type CONST.
 *  @param fl factor loading of the underlying LMM, must have {@link CONST_VolSurface}.
 *  @param t lattice is built until Libor reset time T_t.
 *  @param steps number of time steps in each Libor accrual interval.
 *  @param verbose report details on lattice being built.
 */
LmmLattice2F
(LiborFactorLoading* fl, int t, int steps=1,bool verbose=false);

/** The number of factors. */
int nFactors(){ return 2; }
	
	
// SAMPLE LATTICE

/** A sample three factor lattice in dimension n (number of Libor accrual periods)
 *  built up to time T_p with nSteps time steps per accrual period.
 *  @param verbose details on lattice during build.
 */
static LmmLattice2F* sample(int n, int p, int nSteps, bool verbose=false);

/** Builds a lattice in in dimension n (number of Libor accrual periods)
 *  build with n time steps (one per accrual period) and runs the selfTest().
 */
static void test(int n);
		

private:

// build lattice with m time steps
void buildLattice(int m, bool verbose);

	

}; // end LmmLattice2F



/**********************************************************************************
 *
 *            CONSTANT VOLATILITY THREE FACTOR LMM LATTICE
 *
 *********************************************************************************/


//<----------------To Do----------------->
// make this a class template depending on the Node type (lite/heavy).



/** <p><a href="lmm-lattice-3f">Three factor lattice</a> for driftless Libor market model  
 *  {@link DriftlessLMM} with constant volatility functions \f$\sigma_j(t)=\sigma_j\f$. 
 *
 * <p>It is assumed that all Libor accrual intervals have equal length \f$\delta\f$..
 *  The lattice makes nSteps time steps of equal length \f$dt=\delta/nSteps\f$ in 
 *  each accrual interval. Consequently discrete time t corresponds to continuous time
 *  \f[t*dt=T_{t/nSteps}+(t\;mod\;nSteps)*dt.\f]
 */
class LmmLattice3F : public LmmLattice {
	
public:

	
// CONSTRUCTOR
	
/** The {@link VolSurface} of the factor loading must be of type CONST.
 *  @param fl factor loading of the underlying LMM, must have {@link CONST_VolSurface}.
 *  @param t t lattice is built until Libor reset time T_t.
 *  @param steps number of time seps in each Libor accrual interval.
 *  @param verbose report details on lattice being built.
 */
LmmLattice3F
(LiborFactorLoading* fl, int t, int steps=1, bool verbose=false);


/** The number of factors. */
int nFactors(){ return 3; }
	
	
// SAMPLE LATTICE AND TEST

/** A sample three factor lattice in dimension n (number of Libor accrual periods)
 *  built up to time T_p with nSteps time steps per accrual period.
 *  @param verbose details on lattice during build.
 */
static LmmLattice3F* sample(int n, int p, int nSteps, bool verbose=false);

/** Builds a lattice in in dimension n (number of Libor accrual periods)
 *  build with n-3 time steps (one per accrual period) and runs the selfTest().
 */
static void test(int n);
		
		
private:
	
// build lattice with m time steps.
void buildLattice(int m, bool verbose);

	
}; // end LmmLattice3F





MTGL_END_NAMESPACE(Martingale)

#endif
 
