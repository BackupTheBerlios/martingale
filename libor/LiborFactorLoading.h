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
	Correlation&  rho;         // correlations *corrs
		
public:
	
    
// ACCESSORS
    
    /** Number <code>n</code> of forward Libors including \f$L_0(t)\f$.
     */
    int getDimension() const { return n; }
	
    /** Array of accrual intervals \f$\delta_j\f$.
     */
    RealArray1D& getDeltas() const { return delta; }
    
    /** Array of continuous Libor reset times \f$T_j\f$.
     */
    RealArray1D& getTenorStructure() const { return T; }
	
	/** The array of initial Libors \f$L_j(0)\f$. 
     */
    RealArray1D& getInitialLibors() const { return l; }
	
    /** The array of initial XLibors \f$X_j(0)=\delta_jL_j(0)\f$. 
     */
    RealArray1D& getInitialXLibors() const { return x; }
	
	/** The instantaneous log-Libor correlations.
	 */
	UTRMatrix<Real>& getRho() const { return corr->getCorrelationMatrix(); }
	    
    
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
	
	
// SET THE PARAMETERS (CALIBRATION)
	
   /** Set the parameters of the factorloading from the vector x.
    *  The first coordinates of x populate the scaling factors k,
	*  the rest goes to VolSurface and Correlations.
	*/
   void setParameters(const RealArray1D& x) = 0;
	
   
	
	
// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

   
   /** Volatility \f$\sigma_i(t)\f$ of forward Libor \f$L_i(t)\f$.
	*  See book, 6.4.
    *
    *@param i Libor index.
    *@param t continuous time.
    */
   Real sigma(int i, Real t) const { return k[i]*(vol->sigma(t,T[i])); }
	   
   
   /** The integral
    *  \f[\int_s^t\sigma_i(u)\sigma_j(u)\rho_{ij}du=\langle log(L_i),log(L_j)\rangle_s^t\f]
    *  neeeded for the distribution of time step increments s->t. See book, 6.5.1.
    *
    *@param i,j forward Libor indices.
    *@param s,t continuous integration bounds (times).
    */
   Real integral_sgi_sgj_rhoij(int i, int j, Real s, Real t) const;
	   
   
   /** Annualized volatility of Libor \f$L_i\f$. */
   Real annualVol(int i);
   
   
   /** Message and fields.*/
   std::ostream& printSelf(std::ostream& os);
   
   

// LOG-LIBOR COVARIATION MATRICES 

   
  /** The upper triangular half of the covariation matrix 
   *  <a href="CVpqtT">CV(p,q,s,t)</a>.
   *
   * @param p,q Libors \f$L_j, j=p,...,q-1\f$.
   * @param s,t continuous times \f$0\leq s<t\f$. 
   */ 
   UTRMatrix<Real>&
   logLiborCovariationMatrix(int p,int q, Real s, Real t) const;
  
  
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
	 *  eigenvalues of the covariation matrix <a href="CVpqtT">CV(p,q,s,t)</a>.
	 */
	void factorAnalysis(int p, int q, Real s, Real t)
    {
		TnT::factorAnalysis(logLiborCovariationMatrix(p,q,s,t));
	}
    
    
	/** Prints how much variability is captured by the 5 largest
	 *  eigenvalues of the covariation matrix <a href="CVt">CV(t)</a>.
	 */
	void factorAnalysis(int t)
    {
		TnT::factorAnalysis(logLiborCovariationMatrix(t));
	}

 
   
   

// TEST PROGRAM 

		 
/** Test the roots of all {@link logLiborCovariationMatrix(int t)}.
 */
void selfTest();


/** Computes approximate rank r factorizations C(t)=R(t)R(t)' for all the 
 *  matrices {@link logLiborCovariationMatrix(int t)} and prints the relative error.
 */
void factorizationTest(int r);


/** Returns sample factor loading in dimension n.
 *  @param volType type of VolSurface (VolSurface::CONST, JR, M)
 *  @param corrType type of Correlation (Correlation::JR, CS)
 */
static LiborFactorLoading* 
sample(int n, int volType=VolSurface::JR, int corrType=Correlation::CS)

   
}; // end LiborFactorLoading


// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, LiborFactorLoading* fl){ return fl->printSelf(os); }





/*******************************************************************************
 *
 *                    VOLATILITY SURFACE
 *
 ******************************************************************************/


/** <p>Deterministic volatility surface \f$\sigma(t,T)\f$ for a LiborFactorLoading.
 *  The actual volatility function of Libor \f$X_j\f$ is then given as
 *  \f[\sigma_j(t)=k_j\sigma(t,T_j),\quad j=1,\dots,n-1.\f]
 *  where \f$K-j\f$ is a scaling factor to calibrate the annualized volatility
 *  of \f$X_j\f$. In the case of a DriftlessLMM this function is the volatility of 
 *  forward transported Libor \f$U_j\f$.
 *
 *  <p>To simplify calibration we assume that all volatility surface types
 *  depend on four Real parameters a,b,c,d. This suffices for all types which 
 *  are implemented and the maximal parameter set becomes part of the interface.
 *  In general this is objectionable but here it simplifies the calibration of
 *  factorloadings. The calibration routines use the maximal (and hence the same)
 *  set of parameters for all types of surfaces even though any particular type 
 *  may not depend on all parameters. 
 *  
 * <p>The entire Libor volatility structure then depends on the vector k
 *  of scaling factors and the parameters a,b,c,d, an n+3 dimensional
 *  vector. The vector k is part of the {@link LiborFactorLoading}.
 */
class VolSurface {
	
protected:
	
	/** Flag for the type of vol surface: M,JR,CT.
	 */
	int volType;
	
	Real a,b,c,d;           // the parameters

	
public:
	
		
	/** Flags identifying the volatility surface type:<br>
	 *  M: book, 6.11.1.<br>
	 *  JR: Jaeckel-Rebonato, Jaeckel book, p164, 12.11.<br>
	 *  CT: constant.
	 */
	static const int M=0, JR=1, CONST=2; 
	
	
	/** @param a,b,c,d parameters of the vol surface
	 *  @param type vol surface type M, JR or CONST.
	 */
	VolSurface
	(Real _a, Real _b, Real _c, Real _d, int type) :
	a(_a), b(_b), c(_c), d(_d), volType(type)
    {    }
	
	/** Type ID, implemented: M, JR, CONST.
	 */
	int getType() const { return volType; }
	
	
	/** A sample surface of the type = M, JR, CONST. 
	 */
	static VolSurface* sample(int type);
	

	
	/** The volatility surface \f$\sigma(t,T)\f$. 
	 */
	virtual Real sigma(Real t, Real T) const = 0;

	
	/** The indefinite integral 
	 *  \f[\int\sigma(t,T_1)\sigma(t,T_2)ds.\f]
	 *  Needed for covariation integrals.
	 */
	virtual Real integral_sgsg(Real t, Real T_1, Real T_2) const = 0;
	
	
	/** The definite integral 
	 *  \f[\int_s^t\sigma(u,T_1)\sigma(u,T_2)du.\f]
	 *  Needed for covariation integrals.
	 */
	Real integral_sgsg(Real s, Real t, Real T_1, Real T_2)
    {  return integral_sgsg(t,T1,T2)-integral_sgsg(s,T1,T2); }
	
	
					
	/** Sets the parameters (calibration).
	 */
	void setParameters(Real _a, Real _b, Real _c, Real _d)
	{ a=_a; b=_b; c=_c; d=_d; }
	
	/** Message and fields.*/
	virtual std::ostream& printSelf(std::ostream& os) = 0;
	

}; // end VolSurface



// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, VolSurface* vols){ return vols->printSelf(os); }
	





// JAECKEL-REBONATO VOL SURFACE

/** Jaeckel-Rebonato volatility surface. See Jaeckel book, p164, 12.11.
 */
class JR_VolSurface : public VolSurface {
	
public:
	
	JR_VolSurface
	(Real _a, Real _b, Real _c, Real _d) :
	VolatilitySurface(_a,_b,_c,_d,JR)
    {    }

  
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
	
	
   /** Volatility surface 
    *  \$[\sigma(t,T)=d+(a+b(T-t)\,e^{-c(T-t)})\f$.
    */
   Real sigma(Real t, Real T) const 
   { 
	   Real s=(T-t);
       return d+(a+b*s)*exp(-c*s));
   }

  
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2) const
   {
      return Fij(i,j,s)-Fij(i,j,t);
   }
   
   
   	/** Message and fields.*/
	std::ostream& printSelf(std::ostream& os)
    {
   		return
		os << "\n\nVolSurface, type Jaeckel-Rebonato, " 
		   << "a=" << a << ", b=" << b << ", c=" << c << ", d=" << d 
		   << endl;
    }
	
	/** Sample surface.*/
	static VolSurface* sample()
	{ 
		Real _a=-0.05, _b=0.5, _c=1.5, _d=0.15; 
		return new JR_VolSurface(_a,_b,_c,_d); 
	}

   
 
private:
   
    /** Indefinite integral \f$\int\sigma_i(t)\sigma_j(t)dt\f$.
     */
   Real Fij(int i, int j, Real t) const;
   
	
}; // end JR_VolatilitySurface




// OUR OWN VOL SURFACE

/** Our own volatility surface. See book, 6.11.1.
 *  \f[\sigma(t,T)=g(1-t/T), \quad g(s)=1+ase^{-s/d}.\f]
 */
class M_VolSurface : public VolSurface {
	
public:
	
	M_VolSurface
	(Real _a, Real _b, Real _c, Real _d) :
	VolatilitySurface(_a,_b,_c,_d,M)
    {    }

  
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
	
	
   /** Volatility surface (book 6.11.1)
    *  \f[\sigma(t,T)=g(1-t/T), \quad g(s)=1+ase^{-s/d}.\f]
    */
   Real sigma(Real t, Real T) const { return g(1-t/T); }

  
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2);
   
    
   /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os)
    {
   		return os << "\n\nVolSurface, type book, 6.11.1, " 
		          << "a=" << a << ", d=" << d << endl;
    }
	
	/** Sample surface.*/
	static VolSurface* sample()
	{ 
		Real _a=1.5, _d=2.0; 
		return new M_VolSurface(_a,0.0,0.0,_d); 
	}

   
 
private:
   
	

// AUXILLIARY FUNCTIONS
 
// Exponential Integrals
    
    // int exp(s/D)ds 
    Real F(Real D, Real s) const { return D*exp(s/D); }  
    
    // int s*exp(s/D)ds
    Real G(Real D, Real s) const { return D*exp(s/D)*(s-D); }   

    // int s^2*exp(s/D)ds 
    Real H(Real D, Real s) const { return D*exp(s/D)*((s-D)*(s-D)+D*D); }  
         
    /** The function g(x) defining the volatilities <code>sigma_j</code> 
     *  as \f$\sigma_j(t)=c_jg(1-t/T_j)\f$. See book, 6.11.1
     */
    Real g(Real x) const { return 1+a*x*exp(-x/d); }
    
		
		
}; // end M_VolatilitySurface



// CONSTANT VOLATILITY SURFACE


/** Constant volatility surface. 
 */
class CONST_VolSurface : public VolSurface {
	
public:
	
	CONST_VolSurface(Real _a, Real _b, Real _c, Real _d) :
	VolatilitySurface(_a,_b,_c,_d,CONST)
    {    }

  
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
	
	
   /** Volatility surface \f$\sigma(t,T)=1.\f$.
    */
   Real sigma(Real t, Real T) const { return 1; }

   
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2) const { return t; }
   
   
    /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os)
    {
   		return os << "\n\nVolSurface, type: constant." << endl;
    }
	
	/** Sample (there is only one such surface).
	 */
	static VolSurface* sample(){ return new CONST_VolSurface(0.0,0.0,0.0,0.0); }


	
}; // end CONST_VolSurface





/*******************************************************************************
 *
 *                    CORRELATION
 *
 ******************************************************************************/

/** <p>Constant instantaneous correlations of asset returns, log-Libors, ...
 *  To simplify calibration (the same parameters for all types) we assume that
 *  the correlations depend on three Real paramters 
 *  \f[\alpha, \beta, r_\infty.\f]
 *
 * If more are needed the calibration routines will have to be rewritten in trivial 
 * ways. Only Libor correlations are implemented. To cater to this case the
 * correlations are indexed
 * \f[\rho_{ij},\quad 1\leq i,j<n\f]
 * with index base 1. All matrices, vectors respect this convention.
 * The reason is that Libor \f$L_0\f$ is constant and is not evolved.
 */
class Correlations {
	
protected:
	
	/** Index range [1,n). */
	int n;
	
	/** The type of correlations JR, CS
	 */
	int corrType;
	
	Real alpha, beta, r_oo;
	
	UTRMatrix<Real> correlationMatrix;
	
public:
	
	/** Correlation type JR, CS. */
	static const int JR=0, CS=1;
	
	/** CS or JR. */
	int getCorrelationType() const { return corrType; }
	
	/** n: index range [1,n). */
	int getDimension() const { return n; }
	
	
	/** Correlation matrix is allocated but not initialized
	 *  Do this from the concrete subclasses.
	 *  @param correlationType JR or CS.
	 */
	Correlations(int _n, Real _alpha, Real _beta, Real _r_oo, int correlationType) : 
	n(_n), corrType(correlationType),
	alpha(_alpha), beta(_beta), r_oo(_r_oo), correlationMatrix(n-1,1) 
	{   }
	
	
	/** A sample correlation in dimension n of type = JR, CS.
	 */
	static Correlations* sample(int n, int type);
	
		
	/** Update correlation matrix using the current parameters.
	 */
	virtual void setCorrelations() = 0;
		
	/** Correlations \f$\rho_{ij}\f$ for \f$1\leq i,j<n\f$. 
	 */
	Real& operator()(int i, int j){ return correlationMatrix(i,j); }
	
	
	/** The upper half of the correlation matrix
	 */
	UTRMatrix<Real>& getCorrelationMatrix(){ return correlationMatrix; }
	
				
	/** Sets the parameters (calibration).
	 */
	void setParameters(Real _alpha, Real _beta, Real _r_oo)
	{ alpha=_alpha; beta=_beta; r_oo=_r_oo; }
	
	
	/** Message and fields.
	 */
	virtual std::ostream& printSelf(std::ostream& os) = 0;
	
	
}; // end Correlations



// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, Correlations* corrs){ return corrs->printSelf(os); }




// JAECKEL-REBONATO CORRELATIONS

/** Jaeckel-Rebonato correlations of the form 
 *  \f[\rho_{ij}=exp(-\beta(T_j-T_i)),\quad 1\leq i\leq j<n.\f]
 */
class JR_Correlations : public Correlations {
	
	/** Tenor Structure of Liborprocess. */
	RealArray1D T;
	
public:
	
	/** @param T tenor structure of Libor process.
	 *  @param beta see class description.
	 */
	JR_Correlations(const RealArray1D& T, Real beta) : 
	Correlations(T.getDimension(),0.0,beta,0.0,JR) 	
	{ 
		setCorrelations(); 
	}

	   
	/** Updates correlation matrix using the current parameters.
	 */
    void setCorrelations()
	{
       for(int i=1;i<n;i++)
       for(int j=i;j<n;j++) rho(i,j)=exp(beta*(T[i]-T[j])); 
	}
	
    
    /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os)
    {
		return
		os << "\n\nCorrelations: Jaeckel-Rebonato, " 
		   << "beta=" << beta << endl;
	}
	
	/** Sample correlations.
	 *  @param n dimension.
	 *  @param delta Libor accrual interval length.
	 */
	static Correlations* sample(int n, Real delta=0.25)
	{ 
		RealArray1D T(n+1); T[0]=0.0;
		for(int i=0;i<n;i++) T[i+1]=T[i]+delta;
		
		Real beta=0.1;
		return new JR_Corrleations(T,beta); 
	}


}; // end JR_Correlations





// COFFEE-SHOENMAKERS CORRELATIONS

/** Coffee-Shoenmakers correlations. See book, 6.11.1.
 */
class CS_Correlations : public Correlations {

	
public:
	
	CS_Correlations(int _n, Real _alpha, Real _beta, Real _r_oo) : 
	Correlations(_n,_alpha,_beta,_r_oo,CS) 	
	{ 
		setCorrelations(); 
	}

	   
	/** Updates correlation matrix using the current parameters.
	 */
    void setCorrelations();
	
    
    /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os)
    {
		return
		os << "\n\nCorrelations: Coffee-Shoenmakers, " 
		   << "alpha=" << alpha << "beta=" << beta << "r_oo=" << r_oo 
		   << endl;
	}
	
	
	/** Sample correlations.
	 *  @param n dimension.
	 */
	static Correlations* sample(int n)
	{ 
		Real _alpha=1.8, _beta=0.1, _r_oo=0.4;
		return new CS_Correlations(int _n,_alpha,_beta,_r_oo);
	}
	
	
private:
	
	
	/** The convex function f(x) from which the log-Libor correlations
     *  are derived. See book, 6.11.1.
     */
    Real f(Real x) const;


}; // end CS_Correlations







MTGL_END_NAMESPACE(Martingale)

#endif