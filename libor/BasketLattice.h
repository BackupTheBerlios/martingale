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

#ifndef martingale_basketlattice_h    
#define martingale_basketlattice_h

#include "TypedefsMacros.h"
#include "Matrix.h"             // direct members must have fully specified type
#include "Lattice.h"            // base class template
#include "Node.h"               // base class template parameter


MTGL_BEGIN_NAMESPACE(Martingale)



/*! \file BasketLattice.h
 *  Lattices for a vector of assets \f$S_j,\ j=0,\dots,n-1\f$, following a driftless 
 *  dynamics
 *  \f[dS_j(t)=S_j(t)\sigma_ju_j\cdot dW(t)\f]
 *  with constant volatilities \f$\sigma_j\f$ and unit vectors \f$u_j\f$.
 *  Here W(t) is an n-dimensional Brownian motion where n is the number of assets
 *  \f$S_j\f$ in the basket. The driftless nature of the asset
 *  dynamics implies that the asset prices are martingales. This means that we are 
 *  working with asset prices relative to a numeraire evolved in the corresponding
 *  numeraire measure. For example we could be working with discounted prices in the
 *  so called risk neutral probability (numeraire is the bank account) or with forward 
 *  prices at some fixed date T in the forward martingale measure (numeraire is the
 *  zero coupon bond maturing at time T).
 *
 *  <p>The asset price dynamics can be solved explicitly.
 *  Let \f$Y_j=log(S_j)\f$ denote the return on the asset \f$S_j\f$. Then
 *  \f[dY_j(t)=-\sigma_j^2/2dt+u_j\cdot dW(t)\f]
 *  resulting in
 *  \f[Y_j(t)=Y_j(0)-\sigma_j^2t/2+u_j\cdot W(t)=Y_j(0)-\sigma_j^2t/2+V_j(t)\f]
 *  with \f$V_j(t)=u_j\cdot W(t)\f$. Taking exponentials we get
 *  \f[S_j(t)=S_j(0)exp[-\sigma_j^2t/2+V_j(t)].\f]
 *  The unit vectors \f$u_j\f$ determine the instantaneous
 *  correlations of the return increments \f$dY_j\f$ as \f$\rho_{ij}=u_i\cdot u_j\f$,
 *  that is 
 *  \f[\rho=UU',\f]
 *  where U is the matrix with rows \f$r_i(U)=u_i\f$. Typically the unit vectors 
 *  \f$u_i\f$ are not given but the correlation matrix \f$\rho\f$ is a model parameter.
 *  From this the matrix \f$U\f$ can be computed as a pseudo squareroot of \f$\rho\f$.
 *  See book, Appendix A.1. To keep the number of nodes manageable we must 
 *  reduce the number of factors, that is, the dimension of the Brownian motion W.
 *  If we want to run an r-factor model we compute U as an approximate r-factor
 *  pseudo square root of \f$\rho\f$. Then U is an \f$n\times r\f$ matrix and the 
 *  \f$V_j(t)\f$ are approximated as
 *  \f[V_j(t)\simeq\sigma_j\left[U_{j1}Z_1(t)+\dots+U_{jr}Z_r(t)\right],\f]
 *  where the \f$Z_j(t)\f$ are independent one dimensional Brownian motions
 *  (the components of W). Thus the lattice can be built to evolve the Brownian motion 
 *  \f$Z_j(t)\f$ (the factors) and this is the simplest of all lattices.
 *  In practice at most three factors can be handled in the lattice otherwise the
 *  number of nodes explodes. 
 *  
 *  <p>Note that the use of the continuous time drift \f$-\sigma_j^2t/2\f$ 
 *  does not preserve the martingale property of \f$S_j\f$ (absence of arbitrage).
 *  This is not significant since the lattice is an approximation to the continuous time
 *  dynamics in which the martingale property does hold. However the drifts can be fixed
 *  so that the \f$S_j\f$ in the lattice do have the martingale property if this is 
 *  deemed desirable. 
 *
 * <p><a name="lmm-lattice-3f"><B>Three factor basket lattice.</B></a>
 *  the lattice evolves the variables in time steps of equal length and each 
 *  {@link Node} stores the vector of assets \f$S_j\f$ and these are computed from the
 *  volatilities \f$V_j\f$ of the returns \f$Y_j\f$ which are now givn as
 *  \f[V_j(s)\simeq\sigma_j\left[U_{j1}Z_1(s)+U_{j2}Z_2(s)+U_{j3}Z_3(s)\right],\f]
 *  where the matrix U is the approximate root of rank 3 of the correlation matrix 
 *  \f$\rho\f$. The rows of U are scaled backe to unit norm to preserve the 
 *  correlations \f$\rho_{jj}=1\f$.
 *  This diminishes the quality of the approximation \f$\rho\simeq UU'\f$ but
 *  preserves the volatilities \f$\sigma_j\f$ of the \f$V_j\f$. If we don't do this we 
 *  lose volatility. See book, 8.1.1, 8.1.2, for details and notation. 
 *
 * <p>The \f$Z_j(t)\f$ are independent standard Brownian motions which evolve from the state
 *  \f$Z_j(0)=0\f$ and then tick up or down in ticks of size \f$a=\sqrt{dt}\f$ where dt is the 
 *  size of the time step. The state at any node is then given by the triple of integers
 *  (i,j,k) with
 *  \f[Z_1=ia,\quad Z_2=ja,\quad Z_3=ka.\f]
 *  Two factor lattices are completely similar but the implementation is kept separate for
 *  greater clarity.
 *
 * <p><b>Number assets and of nodes.</b>
 * The lattice can handle any number of assets and the number of nodes in the lattice depends 
 * depends on the number of factors and the number of time steps in the lattice.
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
 * With 1GB main memory we can tolerate about 5.3 million nodes
 * corresponding to 250 time steps in a two factor lattice and 3.5 million 
 * nodes corresponding to 60 time steps in a three factor model. 
 */

 
 

/**********************************************************************************
 *
 *         LATTICE FOR A BASKET OF ASSETS WITH CONSTANT VOLATILITIES
 *
 *********************************************************************************/
 
 
// forward declarations
class ConstantFactorLoading;
class std::ostream;
 
 
 
/** Information about the LmmLattice which must be conveyed to the nodes.
 */
struct BasketLatticeData {
	
	int 
	    /** number of assets: */     n,        
	    /** number of time steps:*/  T,   
	    /** number of factors: */    r;
	
	/** time step.*/
	Real timestep;
	/** ticksize a=sqrt(dt) of a standard Brownian motion over a single time step.*/
	Real ticksize;
	/** the volatilities \f$\sigma_j\f$.*/
	const RealArray1D& sg;
	/** the initial values \f$Y_j(0)=log(S_j(0))\f$.*/
	const RealArray1D& log_S0;
	/** drift -sigma_j^2*dt/2 of Y_j over a single time step.*/	
	RealArray1D driftUnit;  
	/** rank r approximate pseudo square root of the log(S_j) covariance matrix.*/
	const RealMatrix& R;     
	

	/** @param steps number of time steps.
	 *  @param dt size of time step.
	 *  @param vols the volatilities \f$\sigma_j\f$.
	 *  @param Y0 the initial values \f$Y_j(0)=log(S_j(0))\f$.
	 *  @param Q rank r approximate pseudo square root of the log(S_j) covariation matrix.
	 */
	BasketLatticeData
	(int steps, Real dt, const RealArray1D& vols,
	 const RealArray1D& Y0, const RealMatrix& Q);

	
}; // end BasketLatticeInfo


/** Lattice for a basket of constant volatility assets.
 *  No decision on the number of factors.
 *  For more details see the file reference BasketLattice.h.
 */
class BasketLattice : public Lattice<BasketNode> {
	
protected:
	
	ConstantFactorLoading* fl;   // the asst factor loading
	
	int n_,       // number of assets
	    T_,       // number of time steps to the horizon
	    r_;       // the number of factors
	    
	Real dt_;     // size of time step
	
	// all indices start at zero
	RealArray1D Y0_;         // initial asset price logs Y0[j]=log(S_j(0))
	RealArray1D sg_;         // volatilities sg[j]=sigma_j
	/** rank r approximate pseudo square root of the log(S_j) covariance matrix.*/
	const RealMatrix& R_; 
	
	// data object handed to nodes
	BasketLatticeData* latticeData;
	
	
public:

// ACCESSORS

    /** Number of assets
	 */
    int getDimension() const { return n_; }
	
	/** The size of the time step. */
	Real getTimeStep(){ return dt_; }
	
    /** The vector of initial asset price logs.
	 */
    RealArray1D& get_log_S0() { return Y0_; }
	
	/** The factor loading of the asset prices.*/
	ConstantFactorLoading* getFactorLoading(){ return fl; }

	
// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
/** Computes the relative errors in the trace norm of the approximate rank r 
 *  factorization rho=RR' of the correlation matrix rho. See book, Appendix A.1.
 */
void testFactorization() const;
	
	
// CONSTRUCTOR
	
/** Use index base zero throughout.
 *  @param r the number of factors.
 *  @param fl the factorloading of the assets.
 *  @param T number of time steps.
 *  @param dt size of time step.
 *  @param S0 initial asset prices.
 *  @param fl the factorloading of the assets.
 */
BasketLattice
(int r, ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0);


std::ostream& printSelf(std::ostream& os) const;

  
}; // end BasketLattice





/**********************************************************************************
 *
 *            TWO FACTOR BASKET LATTICE
 *
 *********************************************************************************/
		

/** <p>Two factor lattice for a basket of constant volatility assets.
 *  You must have at least two assets.
 *  For more details see the file reference for the file BasketLattice.h.
 */
class BasketLattice2F : public BasketLattice {

	
public:
	
// CONSTRUCTOR
	
/** Use index base zero throughout.
 *  @param fl the asset factor loading.
 *  @param T number of time steps.
 *  @param dt size of time step.
 *  @param S0 initial asset prices.
 */
BasketLattice2F
(ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0);


// SAMPLE LATTICE AND TEST
			
/** Sample lattice with n assets and T time steps of size dt=0.1.
 *  All asset prices start at S_j(0)=100.0 and all volatilities are 0.3.
 *  The correlation of returns are given by 
 *  \f$\rho_{ij}=exp(0.2(i-j))\ i\leq j\f$.
 */
static BasketLattice2F* sample(int n,int T);


/** Builds a lattice in in dimension n (number of Libor accrual periods)
 *  build with T time steps (one per accrual period) and runs the selfTest().
 */
static void test(int n, int T);
		

private:
	
// build lattice with m time steps.
void buildLattice(int m);
	

}; // end BasketLattice2F




/**********************************************************************************
 *
 *            THREE FACTOR BASKET LATTICE
 *
 *********************************************************************************/
		

/** <p>Three factor lattice for a basket of constant volatility assets.
 *  You must have at least three assets.
 *  For more details see the file reference for the file BasketLattice.h.
 */
class BasketLattice3F : public BasketLattice {
	

public:
	
// CONSTRUCTOR
	
/** Use index base zero throughout.
 *  @param fl the asset factor loading.
 *  @param T number of time steps.
 *  @param dt size of time step.
 *  @param S0 initial asset prices.
 */
BasketLattice3F
(ConstantFactorLoading* fl, int T, Real dt, RealArray1D S0);


// SAMPLE LATTICE AND TEST

/** Sample lattice with n assets and T time steps of size dt=0.1.
 *  All asset prices start at S_j(0)=100.0 and all volatilities are 0.3.
 *  The correlation of returns are given by 
 *  \f$\rho_{ij}=exp(0.2(i-j))\ i\leq j\f$.
 */
static BasketLattice3F* sample(int n, int T);


/** Builds a lattice in in dimension n (number of Libor accrual periods)
 *  build with T time steps (one per accrual period) and runs the selfTest().
 */
static void test(int n, int T);


private:
	
// build lattice with m time steps.
void buildLattice(int m);
	

}; // end BasketLattice3F

	



MTGL_END_NAMESPACE(Martingale)

#endif
 
