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


#include <string>
#include <sstream>
#include "Matrices.h"
#include "TnT.h"
#include <cmath>

MTGL_BEGIN_NAMESPACE(Martingale)


/*
 * LiborFactorLoading.h
 * 
 * Classes: FactorLoading, CS_FactorLoading
 *
 * Created on March 9, 2003, 6:00 AM
 */


/*******************************************************************************
 *
 *                     Class: FactorLoading
 * 
 ******************************************************************************/


/** <p>{@link FactorLoading} for the log-Libors \f$Y_i(t)=log(L_i(t))\f$
 *  of a Libor Market Model. Differs from a FactorLoading since the 
 *  volatilities \f$\sigma_i(t)\f$ are defined only for \f$t\leq T_i\f$.
 *  The correlation matrix drops the first row and column
 *  (L_0 is constant).Provides additional functionality to compute the 
 *  lognormal approximation.
 * 
 *  <p><a name="CVpqtT"></a>With \f$t\leq T\f$ continuous times, set
 *  \f$cv_{ij}(s)=\sigma_i(s)\sigma_j(s)\rho_{ij}=\nu_i(s)\cdot\nu_j(s)\f$
 *  and 
 *  \f[CV_{ij}(t,T)=\int_t^T cv_{ij}(s)ds=\left\langle Y_i,Y_j\right\rangle_t^T\f]
 *  and let \f$CV(p,q,t,T)\f$ be the matrix
 *  \f[CV(p,q,t,T)=( CV_{ij}(t,T) )_{p\leq i,j<q}.\f]
 *
 *  <p>Here angular brackets denote the covariation process as usual and
 *  \f$CV(p,q,t,T)\f$ is the <i>covariation matrix</i> of the \f$Y_i,\ p\leq i<q\f$ 
 *  on the interval \f$[t,T]\f$. This matrix and its upper triangular root  
 *  are used in the simulation of the time step t->T  for the process
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


	
    int n;         // dimension of Y
    
    Real* delta;   // delta[j]=Tc[j+1]-Tc[j]
    Real* Tc;      // tenor structure Tc[t]=T_t
	Real* l;       // l[j]=L_j(0)
	Real* x;       // x[j]=X_j(0)=delta_jL_j(0)
	
	// drift linearization constants alpha_k, beta_k
	Real* alpha;
	Real* beta;
    
	
public:
	
		
	/** type id */
	static const int 
	/** Coffee-Shoenmakers */
	CS=0, 
	/** Jaeckel-Rebonato */
	JR=1,
	/** Constant Volatility (Jaeckel-Rebonato correlations) */
	CV=2;   
    
// ACCESSORS
    
    /** Number <code>n</code> of forward Libors including \f$L_0(t)\f$.
     */
    int getDimension() const { return n; }
	
    /** Array of accrual intervals \f$\delta_j\f$.
     */
    Real* getDeltas() const { return delta; }
    
    /** Array of continuous Libor reset times \f$T_j\f$.
     */
    Real* getTenorStructure() const { return Tc; }
	
	/** The array of initial Libors \f$L_j(0)\f$. 
     */
    Real* getInitialTermStructure() const { return l; }
	
    /** The array of initial XLibors \f$X_j(0)=\delta_jL_j(0)\f$. 
     */
    Real* getInitialXLibors() const { return x; }
	
	/** Type ID of the default type (Coffee-Shoenmakers). 
	 *  Override as appropriate.
	 */
	virtual int getTypeID(){ return CS; }
	
	/** Drift linearization constants, see book, 6.5.2.
     */
    Real* getAlphas() const { return alpha; }
    
	/** Drift linearization constants, see book, 6.5.2.
     */
    Real* getBetas() const { return beta; }
	

    
    
// CONSTRUCTOR
	
    /** Function \f$f(y)=exp(y)/(1+exp(y))\f$ which is linearized
     *  for the <code>X1</code>-dynamics and needed to define the linearization
     *  constants \f$\alpha_k, \beta_k\f$.
     */
    Real f(Real y)
    {
        Real x=exp(y);
        return x/(1+x);
    }
	
    /** Constructor.
     *
     * @param dim dimension n of Libor process
	 * @param L0 initial libors \f$L0[j]=L_j(0)\f$.
     * @param deltas array of accrual interval lengths \f$\delta_j\f$.
     */
    LiborFactorLoading(int dim, Real* L0, Real* deltas);
	
	
// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

   /** Instantaneous correlation <code>rho_ij</code> of log-Libor increments for
    *  <code>0&lt;=i,j&lt;n</code>. This includes <code>L_0</code> and so takes
    *  the view that all Libors live to the horizon.
    *  See document <i>LiborProcess.ps.</i></p>
    *
    * @param Libor indices \f$i,j\geq 1\f$.
    */
   virtual Real rho(int i, int j) const = 0;
   
   /** Volatility \f$\sigma_i(t)\f$ of forward Libor \f$L_i(t)\f$.
	*  See book, 6.4.
    *
    *@param i Libor index.
    *@param t continuous time.
    */
   virtual Real sigma(int i, Real t) const = 0;
   
   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)\rho_{ij}ds=\langle log(L_i),log(L_j)\rangle_t^T\f]
    *  neeeded for the distribution of time step increments. See book, 6.5.1.
    *
    *@param i,j forward Libor indices
    *@param t,T continuous integration bounds
    */
   virtual Real 
   integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const = 0;
   
   
   /** String containing a message indicating what type of factor loading
    *  it is, all the parameter values.
    */
   virtual string toString() const = 0;
   

// LOG-LIBOR COVARIATION MATRICES AND DRIFT LINEARIZATION MATRICES


  /** The upper triangular half of the correlation matrix 
   *  \f$\rho_{ij}=u_i\cdot u_j,\quad 1\leq i,j<n.\f$. See book, 6.4.
   *  Index base 1, natural indexation.
   */ 
   UTRMatrix<Real>& getRho() const;

   
  /** The upper triangular half of the covariation matrix 
   *  <a href="CVpqtT">CV(p,q,t,T)</a>.
   *
   * @param p,q Libors \f$L_j, j=p,...,q-1\f$.
   * @param t,T time interval \f$[t,T]\f$, continuous time.
   */ 
   UTRMatrix<Real>&
   logLiborCovariationMatrix(int p,int q, Real t, Real T) const;
  
  
  /** The upper triangular half of the covariation matrix 
   *  <a href="CVt">CV(t)</a>.
   *  This is the matrix needed for the drift part of the Libor
   *  process time step \f$T_t\rightarrow T_{t+1}\f$.
   *
   * @param t discrete time.
   */ 
   UTRMatrix<Real>& 
   logLiborCovariationMatrix(int t) const;

   
   
  /** <p>The upper triangular root R of the covariation matrix 
   *  <a href="CVt">CV(t)</a>.
   *  This matrix is needed for the volatility part of the 
   *  Libor process time step \f$T_t\rightarrow T_{t+1}\f$.
   *
   * @param t discrete time.
   */ 
   UTRMatrix<Real>& 
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
   Matrix<Real>& 
   reducedRankLogLiborCovariationMatrixRoot(int t, int r) const;
   

//  EIGEN ANALYSIS OF THE COVARIATION MATRICES

    /** Prints how much variability is captured by the 5 largest
	 *  eigenvalues of the covariation matrix <a href="CVpqtT">CV(p,q,t,T)</a>.
	 */
	void factorAnalysis(int p, int q, Real t, Real T)
    {
		TnT::factorAnalysis(logLiborCovariationMatrix(p,q,t,T));
	}
    
    
	/** Prints how much variability is captured by the 5 largest
	 *  eigenvalues of the covariation matrix <a href="CVt">CV(t)</a>.
	 */
	void factorAnalysis(int t)
    {
		TnT::factorAnalysis(logLiborCovariationMatrix(t));
	}

   
//  FUNCTIONS NEEDED FOR THE LOGNORMAL lIBOR APPROXIMATION
   
  /** <p>The upper triangular matrix <code>A(t)</code>, see
   *  LMM.ps, section 4, eqn.(11). Indices i,j=p,...,n-1
   *  with subscripting operators using index base p (ie.
   *  natural indexing).
   *
   * @param t time.
   * @param p indices i,j=p,...,n-1.
   */ 
   UTRMatrix<Real>& A(Real t, int p) const;

   
   
  /** <p>The upper triangular matrix \f$B=int_0^t A(s,p)ds\f$,
   *  where <code>A</code> is {@link A}. See book, 6.5.2. 
   *  Indices i,j=p,...,n-1 with subscripting operators using index base p 
   *  (ie. natural indexing).</p>
   *
   * @param t time.
   * @param p indices i,j=p,...,n-1.
   */ 
   UTRMatrix<Real>& B(Real t, int p) const;

   

   /** The matrix exponential exp(B(t,p)), see book, 6.5.2.
    */
   UTRMatrix<Real>& eB(Real t, int p) const;
 
 
   /** The matrix inverse exp(B(t,p))^{-1}=exp(-B(t,p)),
    *  where B is {@link B}. See book, 6.5.2.
    */
   UTRMatrix<Real>& eBInverse(Real t, int p) const;

	   
   /** The vector \f$(u_i(t),u_{i+1}(t),...,u_{n-1}(t))\f$. See book, 6.5.2.  
    *  Index base b=i, natural indices when using the subscripting operator.
    */
   vector<Real>& u(Real t, int i);

   

   /** The matrix \f$C=(\sigma_i(t)\sigma_j(t)\rho_{ij})\f$, where \f$i,j>=p\f$. 
    *  The index base is p, that is, natural indices when using the subscripting operator.
    */
   UTRMatrix<Real>& C(Real t, int p);

   
   /** The matrix \f$\nu(t)\f$ of factor loadings. Computed as the upper triangular root 
    *  R of the matrix \f$C=(\sigma_i(t)\sigma_j(t)\rho_{ij})\f$, where \f$i,j>=p\f$.
    *  Satisfies RR'=C. Index base is p, natural indices when using the subscripting operator.
    */
   UTRMatrix<Real>& nu(Real t, int p);

   
   

// TEST PROGRAM 

	
/** Diagnostic function, default implementation does nothing,
 *  handle this from concrete subclasses as desired.
 */
virtual void printFields(){ }

		 
/** Test the roots of all {@link logLiborCovariationMatrix(int t)}
 *  the factor loadings {@link nu(t)}, and the matrix exponentials
 *  {@link eB(int t)} and {@link eBInverse(int t)}, t=0,...,n-1.
 */
void selfTest();


/** Computes approximate rank r factorizations C(t)=R(t)R(t)' for all the 
 *  matrices {@link logLiborCovariationMatrix(int t)} and prints the relative error.
 */
void factorizationTest(int r);


/** Returns sample factor loading in dimension n.
 *  Default implementation returns null.
 */
virtual LiborFactorLoading* sample(int n)
{
	LiborFactorLoading* null=0;
	return null;
}

 
   
}; // end LiborFactorLoading




/*******************************************************************************
 *
 *                       Class: CS_FactorLoading
 *
 ******************************************************************************/



/** Implements the correlation and volatility structure based on ideas of 
 *  B. Coffey and J. Schoenmakers (CS), see book, 6.7.1.
 *
 * @author Michael J. Meyer
 */
class CS_FactorLoading : public LiborFactorLoading {
    
	
	/** Parameters describing the volatility and correlation structure
     *  as described in <i>LiborProcess.ps</i>.
     *
     */
    Real A,D,alpha,beta,rho00;
    
    /** Array of factors <code>c_j</code> defining the volatilities 
     *  <code>sigma_j</code> as <code>sigma_j=c_jg(1-t/T_j)</code>.
     *  See <i>LiborProcess.ps</i>.
     */
    Real* c;
	
	/** Upper triangular matrix of instantaneous Libor correlations 
	 *  corr(i,j)=rho_{ij}.
	 */
	UTRMatrix<Real> corr;
    

public:

// CONSTRUCTOR
    
    /** For the meaning of the parameters A,D,alpha,beta,rho00 see 
     *  book, 6.7.1. 
     *
     * @param dim number n of forward Libors including <code>L_0</code>.
     * @param a,d,Alpha,Beta,Rho00 the parameters A,D,alpha,beta,rho00 
     * defining the volatility and correlation structure.
     * @param C \f$c_j\f$ factors calibrating volatilities to caplet prices
	 * (including the superfluous \f$c_0\f$ to preserve natural indexation).
     * @param l the initial Libors \f$l[j]=L_j(0)\f$.
     * @param deltas array of accrual periods \f$\delta_j\f$.
     */
    CS_FactorLoading
    (int dim, Real a, Real d, Real Alpha, Real Beta, Real Rho00, Real* C, Real* l, Real* deltas);	
	
                                                                                                        
       
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS 

   int getTypeID(){ return CS; }
  
   /** <p>Instantaneous log-Libor correlations <code>rho_ij</code>
    *  for <code>i,j&gt=1</code>.</p>
    *
    *@param i,j Libor indices \f$i,j\geq1\f$.
    */
   Real rho(int i, int j) const;

  
   /** Volatility \f$\sigma_i(t)\f$ of \f$L_i(t)\f$ defined on \f$[0,T_i]\f$.
    *
    *@param i Libor index
    *@param t continuous time.
    */
   Real sigma(int i, Real t) const;


   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)rho_ijds=\langle log(L_i),log(L_j)\rangle_t^T\f]
    *  neeeded for the distribution of time step increments. See book, 6.5.1.
    *
    *@param i,j Libor indices.
    */
   Real integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const;



   
// SAMPLE FACTOR LOADING
    
    /** Provides a sample <code>CS_FactorLoading</code> object of dimension n.
     *  Parameter values are as follows:
     *  \f[\delta_j=0.25, c_j=0.25, A=1.3, D=2, alpha=0.9, beta=0.01, rho00=0.9*exp(-0.02*n).\f]
     * This choice implies annualized Libor volatilities of about 33%
     * and a humped volatility peaking halfway out to the reset date.
     * 
     * @param n dimension of the Libor process.
     * @param delta constant accrual periods (default: 0.25).
     */
    static CS_FactorLoading* sample(int n, Real delta=0.25);   

	
// STRING REPRESENTATION
  
    /** A message what type of factor loading it is, all the parameter values.
     */
    string toString()  const;
	
		
	/** Diagnostic function, prints the values of most fields to
	 *  check proper setup.
	 */
    void printFields();
	

private:
	

// AUXILLIARY FUNCTIONS
 
// Exponential Integrals
    
    // int exp(s/D)ds 
    Real F(Real D, Real s) const { return D*std::exp(s/D); }  
    
    // int s*exp(s/D)ds
    Real G(Real D, Real s) const { return D*std::exp(s/D)*(s-D); }   

    // int s^2*exp(s/D)ds 
    Real H(Real D, Real s) const { return D*std::exp(s/D)*((s-D)*(s-D)+D*D); }  
    
    
    
    /** The convex function f(x) from which the log-Libor correlations
     *  are derived. See book, 6.11.1.
     */
    Real f(Real x) const 
    {
        Real a=alpha/2, b=(beta-alpha)/6,
               c=log(rho00)-(a+b);
              
        return x*(c+x*(a+b*x));  //cx+ax^2+bx^3
    }
    
     
    /** The function g(x) defining the volatilities <code>sigma_j</code> 
     *  as <code>sigma_j=c_jg(1-t/T_j)</code>. See <i>LiborProcess.ps</i>.
     */
    Real g(Real x) const { return 1+A*x*std::exp(-x/D); }
    
    
    /** The integral int_0^T g(s)^2ds 
     */
    Real integral_g_squared(Real T) const;
    
    
    /** int_t^T g(1-s/a)g(1-s/b)ds, TeX document LMM.
     */
    Real integral_gg(Real t, Real T, Real a, Real b) const;        
 

}; // end CS_FactorLoading




              
/*******************************************************************************
 *
 *                     JR_FactorLoading
 *
 *******************************************************************************/        




/** Implements the correlation and volatility structure from Jaeckel's book
 *  <i>Monte Carlo Methods in Finance</i> with constant log-Libor correlations
 *  \f$\rho_{ij}=exp(\beta*(T_i-T_j)),\ i\leq j\f$ 
 *  and deterministic log-Libor volatilities
 *  \f[\sigma_i(t)=k_i[d+(a+bs e^{-cs})],\ s=T_i-t,\ t\in[0,T_i].\f]
 *
 * @author Michael J. Meyer
 */
class JR_FactorLoading : public LiborFactorLoading {
    
    
    /** Parameters describing the volatility and correlation structure.
     */
    Real a,b,c,d, Beta;
 
    
    /** Array of factors <code>k_j</code> defining the volatilities 
     *  <code>sigma_j</code> as <code>sigma_j=k_jg(1-t/T_j)</code>.
     *  See <i>LiborProcess.ps</i>.
     */
    Real* k;
    
	/** Upper triangular matrix of instantaneous Libor correlations 
	 *  corr(i,j)=rho_{ij}.
	 */
	UTRMatrix<Real> corr;

    
public:

   
   int getTypeID(){ return JR; } 
   

// CONSTRUCTOR

    
    /** For the meaning of the parameters a,b,c,d,Beta,k see Jaeckel's book
     *  <i>Monte Carlo Methods in Finance</i>. 
     *
     * @param dim number n of forward Libors including \f$L_0\f$.
     * @param L0 array of initial Libors \f$L0[j]=L_j(0)\f$.
     * @param deltas array of accrual periods \f$\delta_j\f$.
     * @param factors \f$k_j\f$ calibrating the Libor volatilities to caplet prices.
     */
    JR_FactorLoading
    (int dim, Real* L0, Real* deltas, Real a0, Real b0, Real c0, Real d0, Real Beta0, Real* k0);
	

	
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS

	
   /** Instantaneous log-Libor correlations \f$\rho_{ij}\f$
    *
    *@param i,j Libor indices, \f$i,j\geq 1\f$.
    */
   Real rho(int i, int j) const { if(i<=j) return corr(i,j); return corr(j,i); }
 

  
   /** Volatility \f$\sigma_i(t)\f$ of \f$L_i(t)\f$,
    *  \f[\sigma_i(t)=k_i[d+(a+bs e^{-cs})],\ s=T_i-t,\ t\in[0,T_i].\f]
    *
    *@param i Libor index.
    *@param t continuous time.
    */
   Real sigma(int i, Real t) const
   { 
       Real* Tc=getTenorStructure();
	   Real s=(Tc[i]-t);
       return k[i]*(d+(a+b*s)*std::exp(-c*s));
   }
  
   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)rho_ijds=\langle log(L_i),log(L_j)\rangle_t^T\f]
    *  neeeded for the distribution of time step increments. See book 6.5.1.
    *
    * @param i,j Libor indices.
    * @param t,T continous times $t\leq T$.
    */
   Real integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const
   {
      return Fij(i,j,T)-Fij(i,j,t);
   }
   
   
   

// STRING MESSAGE

    /** A message what type of factor loading it is, all the parameter values.
     */
    string toString() const;
   
   
   	/** Diagnostic function, prints the values of most fields to
	 *  check proper setup.
	 */
    void printFields();
	

   
// SAMPLE FACTOR LOADING

    
    /** Provides a sample <code>CS_FactorLoading</code> object of dimension n.
     *  Parameter values are as follows:
     *  \f[\delta_j=0.25, k_j=0.25, a=1.3, b=2, c=0.5, d=1.0, \beta=0.1\f]
	 * Initial Libors are set to \f$L_j(0)=0.04\f$.
     * 
     * @param n dimension of the Libor process.
	 * @param delta constant accrual interval length.
     */
    static JR_FactorLoading* sample(int n, Real delta=0.25);
	


private:
	 
// AUXILLIARY FUNCTIONS

    
    /** Indefinite integral \f$\int\sigma_i(t)\sigma_j(t)\rho_{ij}(t)dt\f$.
     */
   Real Fij(int i, int j, Real t) const;
       

}; // end JR_FactorLoading



              
/*******************************************************************************
 *
 *          Constant Volatility FactorLoading
 *
 *******************************************************************************/        




/** Libor factor loading with constant volatility functions. Implements both
 *  the Cofee-Shoenmakers and Jaeckel-Rebonato log-Libor correlation structures.
 *
 * @author Michael J. Meyer
 */
class ConstVolLiborFactorLoading : public LiborFactorLoading {
    
    
	/** Flags for the CS and JR correlation structures. */
	static const int CS=0, JR=1;
	
    /** Parameters describing the CS and JR correlations.
     */
    Real alpha,beta,r_oo;           
	              
    /** Array of constant log-Libor volatilities \f$sg[j]=\sigma_j\f$.
     */
    Real* sg;
    
	/** Upper triangular matrix of instantaneous Libor correlations 
	 *  corr(i,j)=rho_{ij}.
	 */
	UTRMatrix<Real> corr;

    
public:

 
   int getTypeID(){ return CV; }
   

// CONSTRUCTOR

    
    /** Constant volatility factor loading with JR correlations
	 *  \f$\rho_{ij}=exp(\beta*(T_i-T_j)),\ i\leq j\f$.
     *
     * @param dim number n of forward Libors including \f$L_0\f$.
     * @param L0 array of initial Libors \f$L0[j]=L_j(0)\f$.
     * @param deltas array of accrual periods \f$\delta_j\f$.
     * @param sg log-Libor volatilities.
     */
    ConstVolLiborFactorLoading
    (int dim, Real* L0, Real* deltas, Real _beta, Real* sg);
	
	
	/** Constant volatility factor loading with CS correlations
	 *  described by the parameters alpha, beta, r_oo. See book, 6.11.1.
     *
     * @param dim number n of forward Libors including \f$L_0\f$.
     * @param L0 array of initial Libors \f$L0[j]=L_j(0)\f$.
     * @param deltas array of accrual periods \f$\delta_j\f$.
     * @param sg log-Libor volatilities.
     */
    ConstVolLiborFactorLoading
    (int dim, Real* L0, Real* deltas, Real _alpha, Real _beta, Real _r_oo, Real* sg);
	

	
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS

	
   /** Instantaneous log-Libor correlations \f$\rho_{ij}\f$
    *
    *@param i,j Libor indices, \f$i,j\geq 1\f$.
    */
   Real rho(int i, int j) const { if(i<=j) return corr(i,j); return corr(j,i); }
 

  
   /** Volatility \f$\sigma_i(t)\f$ of \f$L_i(t)\f$,
    *  \f[\sigma_i(t)=k_i[d+(a+bs e^{-cs})],\ s=T_i-t,\ t\in[0,T_i].\f]
    *
    *@param i Libor index.
    *@param t continuous time.
    */
   Real sigma(int i, Real t) const { return sg[i]; }

  
   /** The integral
    *  \f[\int_t^T\sigma_i(s)\sigma_j(s)rho_ijds=\langle log(L_i),log(L_j)\rangle_t^T\f]
    *  neeeded for the distribution of time step increments. See book 6.5.1.
    *
    * @param i,j Libor indices.
    * @param t,T continous times $t\leq T$.
    */
   Real integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const
   {
      return (T-t)*sg[i]*sg[j]*rho(i,j);
   }
   
   
   

// STRING MESSAGE

    /** A message what type of factor loading it is, all the parameter values.
     */
    string toString() const;
   
   
   	/** Diagnostic function, prints the values of most fields to
	 *  check proper setup.
	 */
    void printFields();
	

   
// SAMPLE FACTOR LOADING

    
    /** Provides a sample <code>CS_FactorLoading</code> object of dimension n
     *  with volatilities \f$\sigma_j=0.4\f$. Initial Libors are set to \f$L_j(0)=0.04\f$.
     * 
     * @param n dimension of the Libor process.
	 * @param delta constant accrual interval length.
	 * @param corrs correlation type, must be CS of JR.
     */
    static ConstVolLiborFactorLoading* 
	sample(int n, Real delta=0.25, int corrs=CS);
	


private:
	 
// AUXILLIARY FUNCTIONS

    /** The convex function f(x) from which the log-Libor correlations
     *  are derived. See book, 6.11.1.
     */
    Real f(Real x) const 
    {
        Real a=alpha/2, b=(beta-alpha)/6,
               c=log(r_oo)-(a+b);
              
        return x*(c+x*(a+b*x));  //cx+ax^2+bx^3
    }
    

       

}; // end ConstVolLiborFactorLoading






MTGL_END_NAMESPACE(Martingale)

#endif