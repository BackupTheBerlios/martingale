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


#include "Lattice.h"
#include "LiborFactorLoading.h"
#include "Utils.h"
#include "Array.h"
#include <math.h>


MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file LmmLattice.h
 *  Lattices for the driftless Libor market model {@link DriftlessLMM}. 
 *  Two and three factor lattices are implemented. It is assumed that the
 *  volatilities \f$\sigma_j(t)\f$ of the forward transported Libors \f$U_j\f$
 *  (see book, 6.8) are <i>constant</i> and that all Libor accrual intervals 
 *  have equal length \f$\delta\f$.
 *
 * <p>Each {@link Node} stores the vector of accrual factors
 *  \f[H=(H_t,...,H_n)\f]
 *  as the most suitable way to store the state of the underlying Libor process since
 *  all other quantities can be recovered from these with minimal computational effort.
 *
 * <p><a name="lmm-lattice-3f">Three factor LMM lattice:</a>
 * <p>Nodes in a 3 factor lattice for the Libor market model {@link LmmLattice3F}
 *  compute the \f$H_j\f$ from the volatility parts \f$V_j\f$ of the forward transported
 *  Libors \f$U_j\f$ approximated as
 *  \f[V_j(s)\simeq\sigma_j\left[R_{j1}Z_1(s)+R_{j2}Z_2(s)+R_{j3}Z_3(s)\right],\f]
 *  where the matrix R is the approximate root of rank 3 of the correlation matrix 
 *  \f$\rho\f$ of the underlying driftless LMM: \f$\rho\simeq RR'\f$.
 *  See book, 8.1.2, for details and notation. 
 *
 * <p>The \f$Z_j(t)\f$ are independent standard Brownian motions which evolve from the state
 *  \f$Z_j(0)=0\f$ and then tick up or down in ticks of size \f$a=\sqrt{dt}\f$ where dt is the 
 *  size of the time step. The state at any node is then given by the triple of integers
 *  (i,j,k) with
 *  \f[Z_1=ia,\quad Z_2=ja,\quad Z_3=ka.\f]
 *  Two factor lattices are completely similar but the implementation is kept separate for
 *  greater clarity.
 */




/**********************************************************************************
 *
 *            GENERAL LMM LATTICE
 *
 *********************************************************************************/
 

/** <p>Lattice of the driftless Libor market model. This is a lattice with nodes 
 *  of type LmmNode, see {@link Node}. Time steps equal Libor accrual periods. 
 *  See book 6.8, 8.1.1 for the details and notation. 
 *  Interface and partial implementation. 
 *  No decision is made on number of factors and state variables.
 *  For more details see the file reference for the file LmmLattice.h (click on
 *  "File List").
 */
template<typename LmmNode>
class LmmLattice : public Lattice<LmmNode> {
	
protected:
	
	int n,            // dimension (number of accrual intervals)
	    nSteps;       // number of time steps in each accrual interval
	
	vector<Real> U0;     // U0[j]=U_j(0)
	vector<Real> H0;     // H0[j]=H_j(0)

	// factor loading of the underlying driftless LMM
	LiborFactorLoading* factorLoading; 
	
	// deterministic drifts mu(s,j)=-0.5*integral_0^{\tau(s)}sigma_j^2(s)ds, 
	// natural indexation. tau(s) is the continuous time reached at time step s=0,1,...
	LiborArray2D<Real> mu;
	
	
public:

// ACCESSORS

    /** The factor loading of the underlying LMM.
	 */
    LiborFactorLoading* getFactorLoading(){ return factorLoading; } 

    /** Dimension of underlying LMM (number of accrual periods).
	 */
    int getDimension() const { return n; }
	
	/** The matrix of deterministic drifts
	 *  \f[\mu(j,t)=\mu_j(0,T_t)=-{1\over2}\int_0^{T_t}\sigma_j^2(s)ds,\q t\leq j,\f]
	 *  of the logarithms \f$Y_j=log(U_j)\f$.
	 */
    LTRMatrix<Real>& getDrifts() { return mu; }
	
    /** The vector U0[j]=U_j(0), initial values.
	 */
    vector<Real>& getU0() { return U0; }
	
	/** The vector H0[j]=H_j(0)), initial values.
	 */
    vector<Real>& getH0() { return H0; }
	
	
// CONSTRUCTOR
	
	/** The lattice can only be built up to time \f$t=n-2\f$ at most.
	 *  The last step must be handled differently since only one state variable
	 *  remains to make this step. Will fix this later.
	 *
	 *  @param fl factor loading of the underlying LMM.
	 *  @param s lattice is built for s time steps from time zero.
	 *  @param steps number of equal sized time steps in each Libor accrual interval.
	 */
	LmmLattice(LiborFactorLoading* fl, int s, int steps=1) : Lattice<LmmNode>(s),
	n(fl->getDimension()), nSteps(steps),
	U0(n),
	H0(n+1),
	factorLoading(fl),
	mu(n,nSteps)
    {  
		// set log(U_j(0)), j=0,1,...,n-1
        Real* x=factorLoading->getInitialXLibors();     // x[j]=X_j(0)
        for(int j=0;j<n;j++){ 
			
			// U_j(0)=X_j(0)(1+X_{j+1}(0))...(1+X_{n-1}(0))
			Real Uj0=x[j]; for(int k=j+1;k<n;k++) Uj0*=1+x[k]; 
			U0[j]=Uj0;
		}
		
		// set H_j(0)
        H0[n]=1.0;
	    for(int j=n-1;j>=0;j--) H0[j]=H0[j+1]+U0[j];

		// write the deterministic drifts mu_j(t)=mu_j(0,T_t)
		for(int t=0;t<n-1;t++)
		for(int u=0;u<nSteps;u++)
		for(int j=t+1;j<n;j++){
			
			int s=t*nSteps+u;     // number of time step
			mu(s,j)=-0.5*factorLoading->integral_sgi_sgj_rhoij(j,j,0.0,tau(s));
		}
		
	} // end constructor
	
	
private:
	
	// continuous time reached after time step s=0,1,...
	Real tau(int s)
    {
	    Real* Tc=factorLoading->getTenorStructure();
		int t=s/nSteps;
		Real delta_t=Tc[t+1]-Tc[t];
		
		return Tc[t]+(delta_t*(s%nSteps))/nSteps;
	}
		
		
}; // end LmmLattice





/**********************************************************************************
 *
 *            CONSTANT VOLATILITY TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
		

/** <p><a href="lmm-lattice-3f">Two factor lattice</a> for driftless Libor market model  
 *  {@link DriftlessLMM} with constant volatility functions \f$\sigma_j(t)=\sigma_j\f$. 
 *
 * <p>It is assumed that all Libor accrual intervals have equal length \f$\delta\f$..
 *  The lattice makes nSteps time steps of equal length \f$dt=\delta/nSteps\f$ in 
 *  each accrual interval. Consequently discrete time t corresponds to continuous time
 *  \f[t*dt=T_{t/nSteps}+(t\;mod\;nSteps)*dt.\f]
 *  For more details see the file reference for the file LmmLattice.h (click on
 *  "File List").
 */
class ConstVolLmmLattice2F : public LmmLattice<LmmNode2F> {

	
    int nSteps;         // number of equal sized time steps in each accrual interval
	
	Real  delta,        // constant Libor accrual periods
	      dt,           // size of time steps
	      a;            // tick size of standard Brownian motion
	
	vector<Real> sg;    // constant Y_j - volatilities
	
	
	// approximate rank two root of the correlation matrix rho
	// row index base 1, column indices 0,1, natural indexation.
	Matrix<Real>& R;
	
public:

	
// ACCESSORS
	

	
// CONSTRUCTOR
	
	/** The lattice can only be built up to time \f$t=n-2\f$ at most.
	 *  The last step must be handled differently since only one state variable
	 *  remains to make this step. Will fix this later.
	 *
	 *  @param fl factor loading of the underlying LMM.
	 *  @param t lattice is built up to time \f$T_t\f$.
	 *  @param steps number of time seps in each accrual interval.
	 */
	ConstVolLmmLattice2F(ConstVolLiborFactorLoading* fl, int t, int steps=1);	
	


// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
    /** Computes the relative errors in the trae norm of the approximate rank two 
     *  factorization rho=RR' of the correlation matrix rho. See FastLibor.pdf 8.3.
     */
    void testFactorization();

	
	
	
// TEST

/** Sets up a ConstVolLmmLattice2F in dimension n and runs the {@link #selfTest}.
 */
static void test(int n);
		
		

private:
	
	
	/** Computes the vector H of accrual factors on the node node.
	 */
	void setStateVariables(LmmNode2F* node);
	
	// build lattice with m time steps
	void buildLattice(int m);
	

}; // end LmmLattice2F



/**********************************************************************************
 *
 *            CONSTANT VOLATILITY THREE FACTOR LMM LATTICE
 *
 *********************************************************************************/
		

/** <p><a href="lmm-lattice-3f">Three factor lattice</a> for driftless Libor market model  
 *  {@link DriftlessLMM} with constant volatility functions \f$\sigma_j(t)=\sigma_j\f$. 
 *
 * <p>It is assumed that all Libor accrual intervals have equal length \f$\delta\f$..
 *  The lattice makes nSteps time steps of equal length \f$dt=\delta/nSteps\f$ in 
 *  each accrual interval. Consequently discrete time t corresponds to continuous time
 *  \f[t*dt=T_{t/nSteps}+(t\;mod\;nSteps)*dt.\f]
 */
class ConstVolLmmLattice3F : public LmmLattice<LmmNode3F> {

	
    int nSteps;         // number of equal sized time steps in each accrual interval
	
	Real  delta,        // constant Libor accrual periods
	      dt,           // size of time steps
	      a;            // tick size of standard Brownian motion
	
	vector<Real> sg;    // constant Y_j - volatilities
	
	
	// approximate rank two root of the correlation matrix rho
	// row index base 1, column indices 0,1,2; natural indexation.
	Matrix<Real>& R;
	
public:

	
// ACCESSORS
	

	
// CONSTRUCTOR
	
	/** The lattice can only be built up to time \f$t=n-2\f$ at most.
	 *  The last step must be handled differently since only one state variable
	 *  remains to make this step. Will fix this later.
	 *
	 *  @param fl factor loading of the underlying LMM.
	 *  @param t lattice is built up to time \f$T_t\f$.
	 *  @param steps number of time seps in each accrual interval.
	 */
	ConstVolLmmLattice3F(ConstVolLiborFactorLoading* fl, int t, int steps=1);	
	


// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
    /** Computes the relative errors in the trae norm of the approximate rank two 
     *  factorization rho=RR' of the correlation matrix rho. See FastLibor.pdf 8.3.
     */
    void testFactorization();

	
	
	
// TEST

/** Sets up a ConstVolLmmLattice3F in dimension n and runs the {@link #selfTest}.
 */
static void test(int n);
		
		

private:
	
	
	/** Computes the vector H of accrual factors on the node node.
	 */
	void setStateVariables(LmmNode3F* node);
	
	// build lattice with m time steps.
	void buildLattice(int m);
	



}; // end LmmLattice3F


	



MTGL_END_NAMESPACE(Martingale)

#endif
 
