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




#ifndef martingale_factorloading_h    
#define martingale_factorloading_h

#include "TypedefsMacros.h"
//#include "Matrices.h"
#include "Matrix.h"

MTGL_BEGIN_NAMESPACE(Martingale)



// we are using
class std::ostream;





/*******************************************************************************
 *
 *                     Class: FactorLoading
 * 
 ******************************************************************************/


/** <p>This class provides access to the factor loadings \f$\nu_i(s)\f$
 *  of an Ito integral \f$Y(t)=\int_0^t\nu(s)dW(s)\f$ and related covariation
 *  integrals. See book, 3.11.</p>
 *  
 *  <p> This Ito integral is the unbounded variation part ("volatility part")
 *  of several interesting stochastic processes or their component wise logarithms, 
 *  see {@link GaussianMartingale, LiborMarketModel, BlackScholesBasket}.</p>
 *  
 *  <p>We assume that  \f$\nu\f$ is a matrix with rows \f$\nu_i(s)\f$ of the form
 *  \f$\nu_i(s)=\sigma_i(s)u_i\f$, where \f$\sigma_i(s)\f$ is a scalar function
 *  and \f$u_i\f$ a constant unit vector. Setting \f$\rho_{ij}=u_i\cdot u_j\f$
 *  we have 
 *  \f[\nu_i(s)\cdot\nu_j(s)=\sigma_i(s)\sigma_j(s)\rho_{ij}.\f] 
 *  Here \f$\sigma_i(s)\f$ is the instantaneous volatility of \f$Y_i(s)\f$ and
 *  \f$\rho_{ij}\f$ the instantaneous correlation of 
 *  \f$dY_i(s),dY_j(s),\quad i,j=0,\ldots,n-1\f$.</p>
 *
 *  <p>With \f$t\leq T\f$ continuous times, set
 *  \f$cv_{ij}(s)=\sigma_i(s)\sigma_j(s)\rho_{ij}=\nu_i(s)\cdot\nu_j(s)\f$
 *  and
 *  \f[CV_{ij}(t,T)=\int_t^T cv_{ij}(s)ds=\left\langle Y_i,Y_j\right\rangle_t^T\f]
 *  and let \f$CV(t,T)\f$ be the matrix
 *  \f[CV(p,n,t,T)=( CV_{ij}(t,T) )_{p\leq i,j<n}.\f]
 *
 *  <p>Here angular brackets denote the covariation process as usual and
 *  \f$CV(t,T)\f$ is the <i>covariation matrix</i> of Y on
 *  the interval \f$[t,T]\f$. This matrix and its upper triangular root  
 *  are used in the simulation of the time step t->T  for the process
 *  Y and related processes.
 *
 * <p>The class <code>FactorLoading</code> has abstract methods to compute 
 * various covariation integrals (which obviously depend on the concrete
 * volatility and correlation structure) and uses these to compute the various
 * matrices needed for path simulation.</p>
 *
 * @author  Michael J. Meyer
 */
class FactorLoading {

protected:
	
    int n;         // dimension of Y
    
	
public:
    
// ACCESSORS
    
    /** Number <code>n</code> of forward Libors including 
	 *  <code>L_0</code>.
     */
    int getDimension() const { return n; }

    
    
// CONSTRUCTOR
	
 
	
    /** Constructor.
     * @param dim dimension n of Y.
     */
    FactorLoading(int dim) : n(dim) {  }
	
	
// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

   /** Instantaneous correlations \f$\rho_{ij}\f$ of \f$dY_i\f$ increments.
    */
   virtual Real rho(int i, int j) const = 0;
   
   /** Volatility \f$\sigma_i(t)\f$ of \f$Y_i(t)\f$.
    *
    *@param i components of \f$Y\f$.
    *@param t <i>continuous</i> time.
    */
   virtual Real sigma(int i, Real t) const = 0;
   
   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)\rho_{ij}ds=\langle Y_i,Y_j\rangle_t^T.\f]
    *
    *@param i,j components of \f$Y\f$.
    *@param t,T continuous times.
    */
   virtual Real 
   integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const = 0;
   
   
	/** Message and fields.*/
	virtual std::ostream& printSelf(std::ostream& os) const = 0;
   

// LOG-COVARIATION MATRICES AND DRIFT LINEARIZATION MATRICES

   
  /** The matrix CV of covariations \f$\langle Y_i,Y_j\rangle_t^T\f$
   *  \f[CV_{ij}=\int_t^T\sigma_i(s)\sigma_j(s)\rho_{ij}(s)ds,\f]
   * where \f$p\leq i\leq j<q\f$. Index base p, natural indices.
   */ 
   const UTRRealMatrix&
   covariationMatrix(int p,int q, Real t, Real T) const;
   
   
  /** The full matrix CV of covariations \f$\langle Y_i,Y_j\rangle_t^T\f$
   *  \f[CV_{ij}=\int_t^T\sigma_i(s)\sigma_j(s)\rho_{ij}(s)ds,\f]
   * where \f$0\leq i\leq j<n\f$. Index base 0.
   */ 
   const UTRRealMatrix&
   covariationMatrix(Real t, Real T) const;
   

    
   /** The factor loading \f$\nu(t)\f$ computed as the upper triangular root 
    *  of the matrix \f$C=(\sigma_i(t)\sigma_j(t)\rho_{ij})\f$, where \f$i,j\geq p\f$.
    *  This matrix satisfies \f$RR'=C\f$. The index base is p, that is, natural indices 
    *  when using the subscripting operator.
    */
   const UTRRealMatrix& nu(Real t, int p) const;

   
}; // end FactorLoading


// GLOBAL INSERTION
std::ostream& operator << (std::ostream& os, const FactorLoading& fl);



/*******************************************************************************
 *
 *            Constant volatility-correlation factor loading
 *
 ******************************************************************************/


/** {@link FactorLoading} with constant volatilities \f$\sigma_i\f$ and
 *  correlations \f$\rho_{ij}\f$.
 */
class ConstantFactorLoading : public FactorLoading {
	
   UTRRealMatrix corr;          // the correlations rho_ij
   RealVector sg;               // the volatilities sigma_i
	
public:
	
	/** @param dim dimension
	 *  @param vols constant volatilities \f$\sigma_i\f$.
	 *  @param rho  constant instantaneous correlations \f$\rho_{ij}\f$.of \f$dY_i\f$..
	 */
	ConstantFactorLoading(int dim, const RealVector& vols, const UTRRealMatrix& rho) :
	FactorLoading(dim), corr(rho), sg(vols) {  }

// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

   /** Instantaneous correlation \f$\rho_{ij}\f$ of \f$dY_i\f$ increments.
    */
   Real rho(int i, int j) const;
	 
   
   /** Volatility \f$\sigma_i(t)\f$ of \f$Y_i(t)\f$.
    *
    *@param i components of \f$Y\f$.
    *@param t <i>continuous</i> time.
    */
   Real sigma(int i, Real t) const;
   

// COVARIATION INTEGRALS CONTROLLING THE TIME STEPS
   
   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)\rho_{ij}ds=\langle Y_i,Y_j\rangle_t^T.\f]
    *
    *@param i,j components of \f$Y\f$.
    *@param t,T continuous times.
    */
   Real 
   integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const;
   

   
// PSEUDO SQUARE ROOTS OF THE CORRELATION MATRIX
   
   /** Upper triangular root of the matrix of instantaneous correlations,
    *  Cholesky factorization.
    */
   const UTRRealMatrix& correlationMatrixRoot();
   
   
   /** Rank r approximate pseudo square root of the matrix of instantaneous correlations.
    *  See book, Appendix A.1. This is used to run a low factor approximation to a dynamics
    *  based on this factor loading.
    */
   const RealMatrix& correlationMatrixRankReducedRoot(int r);
   
   
   /** String containing a message indicating what type of factor loading
    *  it is, all the parameter values.
    */
   std::ostream& printSelf(std::ostream& os) const;
	
}; // end ConstantFactorLoading





MTGL_END_NAMESPACE(Martingale)

#endif