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




#ifndef martingale_liborfactorloading_h    
#define martingale_liborfactorloading_h

#include "TypedefsMacros.h"
#include "Matrices.h"
#include "VolatilityAndCorrelation.h"
#include <string>
#include <iostream>

MTGL_BEGIN_NAMESPACE(Martingale)


// we are using
extern Real exp(Real);




/*******************************************************************************
 *
 *                     LiborFactorLoading
 * 
 ******************************************************************************/


/** <p>{@link FactorLoading} for the log-Libors \f$Y_i(t)=log(L_i(t))\f$
 *  of a Libor Market Model. Volatilities \f$\sigma_j(t)=k_j\sigma(t,T_j)\f$
 *  provided by the {@link VolSurface} member \f$\sigma(t,T)\f$. The factor loading
 *  contains the vector k of volatility scaling factors. Log-Libor correlations
 *  \f$\rho_{ij}\f$ provided by the {@link Correlations} member.
 * 
 *  <p><a name="CVpqtT"></a>With \f$t\leq T\f$ continuous times, set
 *  \f$cv_{ij}(u)=\sigma_i(u)\sigma_j(u)\rho_{ij}=\nu_i(u)\cdot\nu_j(u)\f$
 *  and 
 *  \f[CV_{ij}(s,t)=\int_s^t cv_{ij}(u)ds=\left\langle Y_i,Y_j\right\rangle_s^t\f]
 *  and let \f$CV(p,q,s,t)\f$ be the matrix
 *  \f[CV(p,q,s,t)=( CV_{ij}(s,t) )_{p\leq i,j<q}.\f]
 *
 *  <p>Here angular brackets denote the covariation process as usual and
 *  \f$CV(p,q,s,t)\f$ is the <i>covariation matrix</i> of the \f$Y_i,\ p\leq i<q\f$ 
 *  on the interval \f$[s,t]\f$. This matrix and its upper triangular root  
 *  are used in the simulation of the time step s->t  for the process
 *  Y and related processes.
 *
 * <p><a name="CVt"></a>
 * More precisely to simulate the time step \f$Y(T_t)\to Y(T_{t+1})\f$ where t is now 
 * discrete time we need the matrix \f$CV(t)=CV(t+1,n,T_t,T_{t+1})\f$,
 * \f[CV(t)_{ij}=\int_{T_t}^{T_{t+1}}sigma_i(s)\sigma_j(s)\rho_{ij}ds,\quad t+1\leq i,j<n,\f]
 * for the drift step as well as its upper triangular root R satisfying
 * \f$CV(t)=R(t)R(t)'\f$ for the volatility step.
 *
 * <p><a name="low-rank-R(t)"></a>
 * If a low factor model with r factors is desired we need an approximation
 * \f$CV(t)\simeq R(t)R(t)'\f$, where R(t) has rank at most r and more precisely dimension
 * \f$(n-t-1)\times r\f$. The best approximation is obtained by diagonalizing the matrix
 * \f$CV(t)\f$ as \f$CV(t)=H\,diag(\lambda_j)H'\f$ where H is a unitary matrix and 
 * \f$\lambda_1\geq\lambda_2\geq\dots\f$ the eigenvalues of H. One then sets
 * \f[R(t)=(\sqrt{\lambda_1}h_1,\dots,\sqrt{\lambda_r}h_r),\f]
 * where the \f$h_j\f$ are the columns of H.
 *  
 * @author  Michael J. Meyer
 */
class LiborFactorLoading {	

protected:
	
    int n;                     // dimension of Y	
    
    RealArray1D delta;         // delta[j]=T[j+1]-T[j]
    RealArray1D T;             // tenor structure T[t]=T_t
	RealArray1D l;             // l[j]=L_j(0)
	RealArray1D x;             // x[j]=X_j(0)=delta_jL_j(0)
	RealArray1D k;             // volatility scaling factors
	
    VolSurface* vol;           // volatility surface
	Correlations* corr;        // log-Libor correlations
	Correlations&  rho;         // correlations *corrs
		
public:
	
    
// ACCESSORS
	
	std::string factorLoadingType();
    
    /** Number <code>n</code> of forward Libors including \f$L_0(t)\f$.
     */
    int getDimension() const { return n; }
	
    /** Array of accrual intervals \f$\delta_j\f$.
     */
    const RealArray1D& getDeltas() const { return delta; }
    
    /** Array of continuous Libor reset times \f$T_j\f$.
     */
    const RealArray1D& getTenorStructure() const { return T; }
	
	/** Libor reset date \f$T_i\f$. */
	Real getT(int i) const { return T[i]; }
	
	/** The array of initial Libors \f$L_j(0)\f$. 
     */
    const RealArray1D& getInitialLibors() const { return l; }
	
    /** The array of initial XLibors \f$X_j(0)=\delta_jL_j(0)\f$. 
     */
    const RealArray1D& getInitialXLibors() const { return x; }
	
	/** The instantaneous log-Libor correlations.
	 */
	const UTRRealMatrix& getRho() const;
	
	/** Type flag for the volatility surface: VolSurface::JR,M,CONST.
	 */
	int getVolSurfaceType() const;
	
	/** The {@link VolSurface} of the factor loading.*/
	VolSurface* getVolSurface() { return vol; }
	
	/** Type flag for the volatility surface: VolSurface::JR,M,CONST.
	 */
	int getCorrelationType() const;
	
	/** The {@link Correlations} of the factor loading.*/
	Correlations* getCorrelations() { return corr; }
	
	
    
// CONSTRUCTOR

	
    /** 
	 * @param L0 initial libors \f$L0[j]=L_j(0)\f$.
     * @param deltas array of accrual interval lengths \f$\delta_j\f$.
     * @param _k scaling factors for Libor vols.
	 * @param vols volatility surface.
	 * @param corrs correlations.	
     */
    LiborFactorLoading
	(const RealArray1D& L0, const RealArray1D& deltas, const RealArray1D& _k, 
	 VolSurface* vols, Correlations* corrs);
	
	~LiborFactorLoading(){ delete vol; delete corr; }
	
	
// SET THE PARAMETERS (CALIBRATION)
	
   /** Set the parameters of the factorloading from the vector X.
    *  The first coordinates of u populate the scaling factors k,
	*  the rest goes to VolSurface and Correlations.
	*/
   void setParameters(Real* u);
	
   
	
	
// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

   
   /** Volatility \f$\sigma_i(t)\f$ of forward Libor \f$L_i(t)\f$.
	*  See book, 6.4.
    *
    *@param i Libor index.
    *@param t continuous time.
    */
   Real sigma(int i, Real t) const;
	   
   
   /** The integral
    *  \f[\int_s^t\sigma_i(u)\sigma_j(u)\rho_{ij}du=\langle log(L_i),log(L_j)\rangle_s^t\f]
    *  neeeded for the distribution of time step increments s->t. See book, 6.5.1.
    *
    *@param i,j forward Libor indices.
    *@param s,t continuous integration bounds (times).
    */
   Real integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const;
	   
   
   /** Annualized volatility of Libor \f$L_i\f$. */
   Real annualVol(int i) const;
   
   

// LOG-LIBOR COVARIATION MATRICES 

   
  /** The upper triangular half of the covariation matrix 
   *  <a href="CVpqtT">CV(p,q,s,t)</a>.
   *
   * @param p,q Libors \f$L_j, j=p,...,q-1\f$.
   * @param s,t continuous times \f$0\leq s<t\f$. 
   */ 
   const UTRRealMatrix&
   logLiborCovariationMatrix(int p,int q, Real s, Real t) const;
  
  
  /** The upper triangular half of the covariation matrix 
   *  <a href="CVt">CV(t)</a>.
   *  This is the matrix needed for the drift part of the Libor
   *  process time step \f$T_t\rightarrow T_{t+1}\f$.
   *
   * @param t discrete time.
   */ 
   const UTRRealMatrix& 
   logLiborCovariationMatrix(int t) const;

   
   
  /** <p>The upper triangular root R of the covariation matrix 
   *  <a href="CVt">CV(t)</a>.
   *  This matrix is needed for the volatility part of the 
   *  Libor process time step \f$T_t\rightarrow T_{t+1}\f$.
   *
   * @param t discrete time.
   */ 
   const UTRRealMatrix& 
   logLiborCovariationMatrixRoot(int t) const;


  /** <p>Rank r approximate root R of the covariation matrix
   *  <a href="CVt">CV(t)</a>.
   *  Best approximation \f$CV(t)\simeq RR'\f$, with R of rank 
   *  at most r.
   *  This matrix is needed for the volatility part of the 
   *  Libor process time step \f$T_t\rightarrow T_{t+1}\f$.
   *
   * @param t discrete time.
   * @param r rank (number of factors).
   */ 
   const RealMatrix& 
   reducedRankLogLiborCovariationMatrixRoot(int t, int r) const;
   

//  EIGEN ANALYSIS OF THE COVARIATION MATRICES

	 /** Examines how dominant are the first r eigenvalues of the covariation 
	  *  matrix CV(p,q,s,t)=(C_ij) with entries
	  *  \f[C_{ij}=\int_s^t\sigma_i(u)\sigma_j(u)\rho_{ij}du,\quad p\leq i,j<q.\]
	  *  This matrix is needed for the Libor time step s->t.
	  */
     void factorAnalysis(int p, int q, Real s, Real t, int r) const;
	 
	 
	 /** Examines how dominant are the first r eigenvalues of the covariation 
	  *  matrix CV(t)=(C_ij) with entries
	  *  \f[C_{ij}=\int_{T_t}^{T_{t+1}}\sigma_i(u)\sigma_j(u)\rho_{ij}du,\quad t+1\leq i,j<n.\]
	  *  This matrix is needed for the Libor time step \f$T_t\rightarrow T_{t+1}\f$.
	  */
     void factorAnalysis(int t, int r) const;

 
   
   

// TEST PROGRAM 

		 
/** Test the roots of all {@link logLiborCovariationMatrix(int t)}.
 */
void selfTest() const;


/** Computes approximate rank r factorizations C(t)=R(t)R(t)' for all the 
 *  matrices {@link logLiborCovariationMatrix(int t)} and prints the relative error.
 */
void factorizationTest(int r) const;


/** Returns sample factor loading in dimension n.
 *  @param volType type of VolSurface (VolSurface::CONST, JR, M)
 *  @param corrType type of Correlation (Correlation::JR, CS)
 */
static LiborFactorLoading* 
sample(int n, int volType=VolSurface::JR, int corrType=Correlations::CS);

   
/** Message and fields.*/
std::ostream& printSelf(std::ostream& os) const;
   

}; // end LiborFactorLoading



// GLOBAL INSERTION
std::ostream& operator << (std::ostream& os, const LiborFactorLoading& fl);



MTGL_END_NAMESPACE(Martingale)

#endif