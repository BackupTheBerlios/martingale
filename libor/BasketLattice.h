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

MTGL_BEGIN_NAMESPACE(Martingale)




// forward declarations
class Realvector;
class UTRRealMatrix;
class BasketNode2F;
class BasketNode3F;




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
 *  with \f$V_j(t)=u_j\cdot dW(t)\f$. Taking exponentials we get
 *  \f[S_j(t)=S_j(0)exp[-\sigma_j^2t/2+V_j(t)].\f]
 *  The unit vectors \f$u_j\f$ determine the instantaneous
 *  correlations of the return increments \f$dY_j\f$ as \f$\rho_{ij}=u_i\cdot u_j\f$,
 *  that is 
 *  \f[\rho=UU',\f]
 *  where U is the matrix with rows \f$r_i(U)=u_i\f$. Typically the unit vectors 
 *  \f$u_i\f$ are not given but the correlation matrix \f$\rho\f$ is a model parameter.
 *  From this the matrix \f$U\f$ can be computed as a pseudo square root of \f$\rho\f$.
 *  See book, Appendix A.1. If we are willing to work with the full set of n factors
 *  we can compute U using the Cholesky factorization. This gives us some speed gain
 *  because of the triangular nature of the resulting matrix U. However we can also 
 *  reduce the number of factors, that is, the dimension of the Brownian motion W.
 *  If we want to run an r-factor model we compute U as an approximate r-factor
 *  pseudo square root of \f$\rho\$. Then U is an \f$n\times r\f$ matrix and the 
 *  \f$V_j(t)\f$ are approximated as
 *  \f[V_j(t)\simeq\sigma_j\left[U_{j1}Z_1(t)+\dots+U_{jr}Z_r(t),\right]
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
 * <p><a name="lmm-lattice-3f">Three factor basket lattice:</a>
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
 * only on the number of time steps and the number r of factors. However the matrix of 
 * correlations of asset returns will be approximated with a rank r matrix. If we have a large 
 * portfolio of assets it needs to be determined wether such an approximation is reasonable,
 * equivalently how dominant are the r largest eigenvalues of the correlation matrix.
 * See book, Appendix A.1.
 * 
 * <p>Each node has four edges in the case of a two factor lattice and eight edges in the 
 * case of a three factor lattice. We want to have a large number of nodes but not too large. 
 * With 30 time steps a 2 factor lattice allocates about 10000 nodes while a three factor 
 * lattice allocates 250000 nodes.
 */
 
 

/**********************************************************************************
 *
 *         LATTICE FOR A BASKET OF ASSETS WITH CONSTANT VOLATILITIES
 *
 *********************************************************************************/

/** Lattice for a basket of constant volatility assets.
 *  No decision on the number of factors.
 *  For more details see the file reference BasketLattice.h.
 */
template<typename BasketNode>
class BasketLattice : public Lattice<BasketNode> {
	
protected:
	
	int n,       // number of assets
	    T;       // number of time steps to the horizon
	    
	Real dt,     // size of the time step
	     a;      // a=sqrt(dt) tick size of standard Brownian motion
	
	// all indices start at zero
	RealVector S0;         // initial asset prices
	RealVector sg;         // volatilities sg[j]=sigma_j
	RealVector driftunit;  // drift sigma_j^2*dt/2 of Y_j over a single time step
	UTRRealMatrix rho;     // correlation matrix.	
	
public:

// ACCESSORS

    /** Number of assets
	 */
    int getDimension() const { return n; }
	
    /** The vector of initial asset prices.
	 */
    RealVector& getS0() { return S0; }
	
	
// CONSTRUCTOR
	
	/** Use index base zero throughout.
	 *  @param n number of assets
	 *  @param T number of time steps.
	 *  @param dt size of time step.
	 *  @param S0 initial asset prices.
	 *  @param sg asset volatilities.
	 *  @param rho correlation of returns.
	 */
	BasketLattice
	(int _n, int _T, Real _dt, RealVector _S0, RealVector _sg, UTRRealMatrix& _rho) : 
	Lattice<BasketNode>(_T),
	n(_n), T(_T), dt(_dt), a(sqrt(dt)), S0(_S0), sg(_sg), driftunit(_n), rho(_rho)
    {  
		for(int j=0;j<n;j++) driftunit[j]=-sg[j]*sg[j]*dt/2.0;
	}
	
  
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
class BasketLattice2F : public BasketLattice<BasketNode2F> {
	
	
	// approximate rank two root of the correlation matrix rho
	RealMatrix& R;
	
public:


	
// CONSTRUCTOR
	
	/** Use index base zero throughout.
	 *  @param n number of assets
	 *  @param T number of time steps.
	 *  @param dt size of time step.
	 *  @param S0 initial asset prices.
	 *  @param sg asset volatilities.
	 *  @param rho correlation of returns.
	 */
	BasketLattice2F
	(int _n, int _T, Real _dt, RealVector& _S0, RealVector& _sg, UTRRealMatrix& _rho);
	
			
	/** Sample lattice with n assets and T time steps of size dt=0.1.
	 *  All asset prices start at S_j(0)=100.0 and all volatilities are 0.3.
	 *  The correlation of returns are given by 
	 *  \f$\rho_{ij}=exp(0.2(i-j))\ i\leq j\f$.
	 */
	static BasketLattice2F* sample(int n,int T);
	


// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
    /** Computes the relative errors in the trae norm of the approximate rank two 
     *  factorization rho=RR' of the correlation matrix rho. See book, Appendix A.1.
     */
    void testFactorization();

	
	
	
// TEST

   /** Builds sample lattice with n assets and T time steps and runs the self test.
	*/
	static void test(int n, int T);	
		

private:
	
	
	/** Computes the vector H of accrual factors on the node node.
	 */
	void setAssetPrices(BasketNode2F* node);
	
	// build lattice with m time steps
	void buildLattice();
	

}; // end BasketLattice2F



/**********************************************************************************
 *
 *            TWO FACTOR BASKET LATTICE
 *
 *********************************************************************************/
		

/** <p>Three factor lattice for a basket of constant volatility assets.
 *  You must have at least three assets.
 *  For more details see the file reference for the file BasketLattice.h.
 */
class BasketLattice3F : public BasketLattice<BasketNode3F> {
	
	// approximate rank three root of the correlation matrix rho
	RealMatrix& R;
	
public:


	
// CONSTRUCTOR
	
	/** Use index base zero throughout.
	 *  @param n number of assets
	 *  @param T number of time steps.
	 *  @param dt size of time step.
	 *  @param S0 initial asset prices.
	 *  @param sg asset volatilities.
	 *  @param rho correlation of returns.
	 */
	BasketLattice3F
	(int _n, int _T, Real _dt, RealVector& _S0, RealVector& _sg, UTRRealMatrix& _rho);
	
	
	/** Sample lattice with n assets and T time steps of size dt=0.1.
	 *  All asset prices start at S_j(0)=100.0 and all volatilities are 0.3.
	 *  The correlation of returns are given by 
	 *  \f$\rho_{ij}=exp(0.2(i-j))\ i\leq j\f$.
	 */
	static BasketLattice3F* sample(int n, int T);


// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
    /** Computes the relative errors in the trae norm of the approximate rank two 
     *  factorization rho=RR' of the correlation matrix rho. See book, Appendix A.1.
     */
    void testFactorization();

	
	
	
// TEST

   /** Builds sample lattice with n assets and T time steps and runs the self test.
	*/
	static void test(int n, int T);	
		

private:
	
	
	/** Computes the vector H of accrual factors on the node node.
	 */
	void setAssetPrices(BasketNode3F* node);
	
	// build lattice with m time steps
	void buildLattice();
	

}; // end BasketLattice3F


	



MTGL_END_NAMESPACE(Martingale)

#endif
 
