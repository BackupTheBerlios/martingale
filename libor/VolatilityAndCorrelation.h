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




#ifndef martingale_volatilityandcorrelation_h    
#define martingale_volatilityandcorrelation_h

#include "TypedefsMacros.h"
#include "Matrix.h"
#include <string>
#include <iostream>

MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file VolatilityAndCorrelations.h
 * Interface and implementation of deterministic volatility surfaces and
 * constant correlations. These are used in the factor loadings of Ito processes
 * (such as asset or Libor logarithms (returns)), see {@link FactorLoading} and 
 * {@link LiborFactorLoading}.
 */




/*******************************************************************************
 *
 *                    VOLATILITY SURFACE
 *
 ******************************************************************************/


/** <p>Deterministic volatility surface \f$\sigma(t,T)\f$ for a FactorLoading.
 *  The actual volatility function of the process \f$Y_j\f$ (in applications to finance 
 *  usually asset returns, that is asset logarithms or Libor logarithms) is then given as
 *  \f[\sigma_j(t)=k_j\sigma(t,T_j),\quad j=1,\dots,n-1.\f]
 *  where \f$k_j\f$ is a scaling factor to calibrate the annualized volatility
 *  of \f$Y_j\f$. In the case of a {@link DriftlessLMM} this function is the volatility 
 *  of the logarithms forward transported Libor \f$Y_j=log(U_j)\f$.
 *
 *  <p>To simplify calibration we assume that all volatility surface types
 *  depend on four Real parameters a,b,c,d. This suffices for all types which 
 *  are implemented and the maximal parameter set becomes part of the interface.
 *  In general this is objectionable but here it simplifies the calibration of
 *  factorloadings. The calibration routines use the maximal (and hence the same)
 *  set of parameters for all types of surfaces even though any particular type 
 *  may not depend on all parameters. 
 *  
 * <p>For simplicity k[0] is included in the calibration (although it has 
 *  no effect in the case of Libors). The entire volatility structure then depends on 
 *  the vector k of scaling factors and the parameters a,b,c,d, an n+4 dimensional
 *  vector. In the case of Libors the vector k is part of the {@link LiborFactorLoading}.
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
	volType(type), a(_a), b(_b), c(_c), d(_d)
    {    }
	
	/** Type ID, implemented: M, JR, CONST.
	 */
	int getType() const { return volType; }
	
	
	/** A sample surface of the type = M, JR, CONST. 
	 */
	static VolSurface* sample(int type);
	
	/** "CONST", "M" or "JR", converts integer type ID to string. */
	static string volSurfaceType(int type);
	

	
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
	Real integral_sgsg(Real s, Real t, Real T_1, Real T_2) const
    {  return integral_sgsg(t,T_1,T_2)-integral_sgsg(s,T_1,T_2); }
	
	
	
    /** Diagnostic. Monte Carlo test of the volatility integrals 
     *  {@link integral_sgsg(Real,Real,Real,Real)}.
     *  Analytic formulas tested against Monte Carlo over a range of
     *  integration intervals.
     *
	 * @param N number of Monte Carlo sample points.
     * @param precision maximal acceptable relative error in percent.
     */
    void testVolSurfaceIntegrals(int N, Real precision);
	
	
					
	/** Sets the parameters (calibration).
	 */
	void setParameters(Real _a, Real _b, Real _c, Real _d)
	{ a=_a; b=_b; c=_c; d=_d; }
	
	/** Message and fields.*/
	virtual std::ostream& printSelf(std::ostream& os) const = 0;
	
	
	/** Type of volatility surface "CONST, "JR", "M".
	 */
	string volSurfaceType();
	

}; // end VolSurface




// JAECKEL-REBONATO VOL SURFACE

/** Jaeckel-Rebonato volatility surface. See Jaeckel book, p164, 12.11.
 */
class JR_VolSurface : public VolSurface {
	
public:
	
	JR_VolSurface
	(Real _a, Real _b, Real _c, Real _d) :
	VolSurface(_a,_b,_c,_d,JR)
    {    }

  
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
	
	
   /** Volatility surface 
    *  \$[\sigma(t,T)=d+(a+b(T-t)\,e^{-c(T-t)})\f$.
    */
   Real sigma(Real t, Real T) const;

  
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2) const;
   
   
   	/** Message and fields.*/
	std::ostream& printSelf(std::ostream& os) const;
	
	/** Sample surface.*/
	static VolSurface* sample();

   
	
}; // end JR_VolatilitySurface




// OUR OWN VOL SURFACE

/** Our own volatility surface. See book, 6.11.1.
 *  \f[\sigma(t,T)=g(1-t/T), \quad g(s)=1+ase^{-s/d}.\f]
 */
class M_VolSurface : public VolSurface {
	
public:
	
	M_VolSurface
	(Real _a, Real _b, Real _c, Real _d) :
	VolSurface(_a,_b,_c,_d,M)
    {    }

  
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
	
	
   /** Volatility surface (book 6.11.1)
    *  \f[\sigma(t,T)=g(1-t/T), \quad g(s)=1+ase^{-s/d}.\f]
    */
   Real sigma(Real t, Real T) const { return g(1-t/T); }

  
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2) const;
   
    
   /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os) const;
	
	/** Sample surface.*/
	static VolSurface* sample();

  
private:

// Exponential Integrals
    
    // integral of exp(s/D)ds 
    Real F(Real D, Real s) const;
    
    // integral of  s*exp(s/D)ds
    Real G(Real D, Real s) const;

    // integral of  s^2*exp(s/D)ds 
    Real H(Real D, Real s) const;
         
    /** The function g(x) defining the volatilities <code>sigma_j</code> 
     *  as \f$\sigma_j(t)=c_jg(1-t/T_j)\f$. See book, 6.11.1
     */
    Real g(Real x) const;
    
		
		
}; // end M_VolatilitySurface



// CONSTANT VOLATILITY SURFACE


/** Constant volatility surface. 
 */
class CONST_VolSurface : public VolSurface {
	
public:
	
	CONST_VolSurface(Real _a, Real _b, Real _c, Real _d) :
	VolSurface(_a,_b,_c,_d,CONST)
    {    }
		
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS
		
   /** Volatility surface \f$\sigma(t,T)=1.\f$.
    */
   Real sigma(Real t, Real T) const { return 1; }
   
   /** See {@link VolatilitySurface#integral_sgsg}.
    */
   Real integral_sgsg(Real t, Real T_1, Real T_2) const { return t; }
      
   /** Message and fields.*/
   std::ostream& printSelf(std::ostream& os) const;
	
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
	
	UTRRealMatrix correlationMatrix;
	
public:
	
	/** Correlation type JR, CS. */
	static const int JR=0, CS=1;
	
	/** CS or JR. */
	int getType() const { return corrType; }
	
	/** "CS" or "JR", converts integer type ID to string.*/
	static string correlationType(int type);
	
	/** n: index range [1,n). */
	int getDimension() const { return n; }
	
	
	/** Correlation matrix is allocated but not initialized
	 *  Do this from the concrete subclasses.
	 *  @param correlationType JR or CS.
	 */
	Correlations(int _n, Real _alpha, Real _beta, Real _r_oo, int correlationType);
	
	
	/** A sample correlation in dimension n of type = JR, CS.
	 */
	static Correlations* sample(int n, int type);
	
		
	/** Update correlation matrix using the current parameters.
	 */
	virtual void setCorrelations() = 0;
		
	/** Correlations \f$\rho_{ij}\f$ for \f$1\leq i,j<n\f$. 
	 */
	Real& operator()(int i, int j);
	
	
	/** The upper half of the correlation matrix
	 */
	const UTRRealMatrix& getCorrelationMatrix() const
	{ return correlationMatrix; }
	
				
	/** Sets the parameters (calibration).
	 */
	void setParameters(Real _alpha, Real _beta, Real _r_oo);
	
	/** Message and fields.
	 */
	virtual std::ostream& printSelf(std::ostream& os) const = 0;
	
	
	/** Converts integer type ID to string "JR", "CS".
	 */
	std::string correlationType();
	
	
}; // end Correlations




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
	JR_Correlations(const RealArray1D& _T, Real beta);

	   
	/** Updates correlation matrix using the current parameters.
	 */
    void setCorrelations();
	
    
    /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os) const;
	
	/** Sample correlations.
	 *  @param n dimension.
	 *  @param delta Libor accrual interval length.
	 */
	static Correlations* sample(int n, Real delta=0.25);


}; // end JR_Correlations





// COFFEE-SHOENMAKERS CORRELATIONS

/** Coffee-Shoenmakers correlations. See book, 6.11.1.
 */
class CS_Correlations : public Correlations {

	
public:
	
	CS_Correlations(int _n, Real _alpha, Real _beta, Real _r_oo);

	   
	/** Updates correlation matrix using the current parameters.
	 */
    void setCorrelations();
	
    
    /** Message and fields.*/
	std::ostream& printSelf(std::ostream& os) const;
	
	
	/** Sample correlations.
	 *  @param n dimension.
	 */
	static Correlations* sample(int _n);
	
	
private:
	
	
	/** The convex function f(x) from which the log-Libor correlations
     *  are derived. See book, 6.11.1.
     */
    Real f(Real x) const;


}; // end CS_Correlations




MTGL_END_NAMESPACE(Martingale)

#endif