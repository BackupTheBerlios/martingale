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
#include <math.h>


MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file LmmLattice.h
 *  Lattices for the Libor market model.
 */



/**********************************************************************************
 *
 *            GENERAL LMM LATTICE
 *
 *********************************************************************************/
 

/** <p>Lattice of the driftless Libor market model. This is a lattice with nodes 
 *  of type {@link LmmNode}. Time steps equal Libor accrual periods. 
 *  See book 6.8, 8.1.1 for the details and notation. 
 *  Interface and partial implementation. 
 *  No decision is made on number of factors and state variables.
 */
template<typename LmmNode>
class LmmLattice : public Lattice<LmmNode> {
	
protected:
	
	int n;            // dimension (number of accrual intervals)
	
	vector<Real> U0;     // U0[j]=U_j(0)
	vector<Real> H0;     // H0[j]=H_j(0)

	// factor loading of the underlying driftless LMM
	LiborFactorLoading* factorLoading; 
	
	// deterministic drifts mu(j,t)=-0.5*integral_0^{T_t}sigma_j^2(s)ds, 
	// 1<=j<=t<=n-1, index base 1, natural indexation
	LTRMatrix<Real> mu;
	
	
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
	 *  @param s lattice is built up to discrete time \f$t=s\leq n-2\f$ inclusively.
	 */
	LmmLattice(LiborFactorLoading* fl, int s) : Lattice<LmmNode>(s),
	n(fl->getDimension()), 
	U0(n),
	H0(n+1),
	factorLoading(fl),
	mu(n-1,1)
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

		Real* Tc=factorLoading->getTenorStructure();
		// write the deterministic drifts mu_j(t)=mu_j(0,T_t)
		for(int j=1;j<n;j++)
		for(int t=1;t<=j;t++) 
			mu(j,t)=-0.5*factorLoading->integral_sgi_sgj_rhoij(j,j,0.0,Tc[t]);
		
	} // end constructor
	
		
		
}; // end LmmLattice




/**********************************************************************************
 *
 *            TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
 	

/** <p>Two factor lattice of the driftless Libor market model.
 *  See book 6.8, 8.1.1 for the details and notation. The variables
 *  evolved in the tree are \f$V1=V_{n-1},\ V2=V_{n-2}-V_{n-1}\f$ and 
 *  \f$t=0,1,\dots,n-2\f$ indexes the time steps \f$T_t\to T_{t+1}\f$.
 *  Current implementation constructs the lattice only up to discrete
 *  time \f$t=n-2\f$.
 *
 * <p>The lattice evolves the state variables \f$U_j\f$ of the driftless 
 * Libor market model but the nodes store the accrual factors \f$H_j\f$
 * since these allow us to compute Libors, swaprates, etc with the minimum effort.
 * See book, 6.8. 
 * 
 * <p>The variables \f$V_j\f$ are the unbounded variation parts of the
 * logarithms \f$log(U_j)\f$. Only these need to be modelled after the standard
 * simplifications See book 3.11, 8.1.1.
 */
class LmmLattice2F : public LmmLattice<LmmNode2F> {
	
// <------------- TO DO ------------> 
// work with 2 independent tick sizes	
	Real a;            // minimum tick size a_1=a_2=a
	
	vector<int> k1, k2;      // k_j(t)    
	vector<int> K1, K2;      // K_j(t)
	
	vector<Real> D11, D22, D12;  // D_{ij(t)
	vector<Real> Q11, Q22, Q12;  // Q_{ij(t)
	
	
	// the sequence of rank 2 roots R(t) of the covariation matrices C(t) 
	// of the vectors V(t)=(V_t(T_t),...,V_{n-1}(T_t)).
	vector< Matrix<Real>* > R;
	
    // RQ[t] is the 2 by 2 inverse of the last two rows of R[t].
	vector< Matrix<Real>* > RQ;
	
	
public:

// ACCESSORS
	
	/** The rank 2 approximate root R(t) of the covariation matrix C(t) 
	 *  of the vector V(t)=(V_t(T_t),...,V_{n-1}(T_t)). See book 8.1.1.
	 */
    Matrix<Real>& getR(int t) { return *R[t]; }
	
    /** The 2 by 2 inverse of the last two rows of
	 *  {@link #getR(int t)}.
	 */
    Matrix<Real>& getRQ(int t) { return *RQ[t]; }
	

	
// CONSTRUCTOR
	
	/** The lattice can only be built up to time \f$t=n-2\f$ at most.
	 *  The last step must be handled differently since only one state variable
	 *  remains to make this step. Will fix this later.
	 *
	 *  @param fl factor loading of the underlying LMM.
	 *  @param s lattice is built up to discrete time \f$t=s\leq n-2\f$ inclusively.
	 */
	LmmLattice2F(LiborFactorLoading* fl, int s);	
	


// ERROR IN THE RANK 2 FATCORIZATION OF THE COVARIANCE MATRICES
	
    /** Computes the relative errors in the trae norm of the	approximate rank two 
     *  factorization C(t)=R(t)R(t)' of the rlevant covariance matrices C(t). See
     *  Book, 8.1.2.
     */
    void testFactorization();

	
	
	
// TEST

/** Sets up an LmmLattice2F in dimension n and runs the {@link #selfTest}.
 */
static void test(int n);
		
		

private:
	
	/** The covariation matrix of the vector
	 *  \f[\tilde V(t)=(V_t,V_{t+1},\dots,V_{n-3},V_{n-2}-V_{n-1},V_{n-1})\f]
	 * evaluatd at time \f$T_t\f$. The last two coordinates have been decorrelated.
	 */
	UTRMatrix<Real>& decorrelatedCovariationMatrix(int t);	
			
				
	/** <p>Computes the vector of matrices R=(R(0),R(1),...,R(m)), where R(t) is 
	 *  the rank 2 approximate root of the covariation matrix C(t) of the vector
	 *  \f[\tilde V(t)=(V_t,...,V_{n-3},V_{n-2}-V_{n-1},V_{n-1})\f]
	 *  at time \f$T_t\f$. This matrix is needed to compute
	 *  the vector V(t) from the variables V1, V2 at the nodes at time t.
	 *
	 *  <p>Also computes the sequence of matrices RQ[t] (the 2 by 2 inverse of the last 
	 *  two rows of R[t]). Assumes \f$m\leq n-2\f$.
	 */
	void setRank2CovariationMatrixRoots(int m);	


	/** Computes the tick size a>0 and sets the corresponding quantities Q_ij(t).
	 * @param m lattice is built up to time t=m inclusively.
	 */
	void setTickSize(int m);
	
	/** Computes the vector H of accrual factors on the node node.
	 */
	void setStateVariables(LmmNode2F* node);
	
	// build lattice up to discrete time t=m (inclusive)
	void buildLattice(int m);
	
	
	// the transition probability p_{ij}
	Real probability(int t, int i, int j);
	


}; // end LmmLattice2F



/**********************************************************************************
 *
 *            CONSTANT VOLATILITY TWO FACTOR LMM LATTICE
 *
 *********************************************************************************/
		

/** <p>Two factor lattice for driftless Libor market model with constant 
 *  volatility functions \f$\sigma_j(t)=\sigma_j\f$. In this case the variables
 *  \f$V_j\f$ are given as
 *  \f[V_j(t)=\sigma_jv_j\cdot Z(t),\f]
 *  where the vectors $v_j$ are the rows $v_j=row_j(R)$ of the approximate root R of
 *  rank 2
 *  \f[\rho\simeq RR'\f]
 *  of the correlation matrix \f$\rho\f$ of the LiborFactorLoading of the underlying
 *  DriftlessLMM and \f$Z(t)=(Z_1(t),Z_2(t))\f$ is a standard two dimensional Brownian motion.
 *  The lattice evolves the variables \f$Z_1(t),Z_2(t)\f$ (called V1,V2 in the nodes).
 *  See FastLibor.pdf  8.1, 8.2, 8.3 for the details and notation.
 *
 * <p>This class is the basic proof of concept for our attempts to evolve only
 * two state variables and to recover the remaining ones according to the joint
 * distribution. If we can't make it happen in this case it won't happen at all.
 */
class ConstVolLmmLattice2F : public LmmLattice<LmmNode2F> {
	

	Real  delta,        // constant Libor accrual periods
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
	 *  @param s lattice is built up to discrete time \f$t=s\leq n-2\f$ inclusively.
	 */
	ConstVolLmmLattice2F(ConstVolLiborFactorLoading* fl, int s);	
	


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
	
	// build lattice up to discrete time t=m (inclusive)
	void buildLattice(int m);
	



}; // end LmmLattice2F


	



MTGL_END_NAMESPACE(Martingale)

#endif
 
