 /***************************************************************************
 *            gpr.h
 *
 *  Sat May 29 16:37:19 2004
 *  Copyright  2004  cpp
 *  cpp@linux
 ****************************************************************************/


#ifndef gpr_h
#define gpr_h


#include "Matrix.h"
#include "BasisFunctions.h"
#include "FunctionExamples.h"
#include <string>
using namespace Martingale;
using std::string;




/** Gaussian process regression on [-1,+1] using expansion in
 *  basis functions. See gprs.ps.
 **/
class GPR {
	
public:

   /** Flag, center of Gaussian prior (origin or empirical coefficients).
    **/
   enum RegressionType {EMPIRICAL, GAUSSIAN};

   /** Set the type of regression (EMPIRICAL, GAUSSIAN).
    **/
   void setRegressionType(RegressionType& rt){ regrType=rt; }

   /** Set the mean of the Gaussian prior P equal to mean.
    **/
   void setPriorMean(RealArray1D& mean)
   { mu=mean; recomputeGaussianCoefficients(); }

   /** Set the mean of the Gaussian prior P equal to mean.
    **/
   void setPriorMeanToOrigin()
   {
      RealArray1D origin(n+1); mu=origin;
      recomputeGaussianCoefficients();
   }

   /** Set the mean of the Gaussian prior P equal to mean.
    **/
   void setPriorMeanToEmpiricalCoefficients()
   { mu=empCoeff; recomputeGaussianCoefficients(); }

   void setBasisFunctions(BasisFunctions* bFcns)
   { delete basis; basis=bFcns; fullInitialization(); }

   /** The sequence \f$(\psi_0(t),\psi_1(t),\dots,\psi_m(t))\f$ of the basis
    *  functions evaluated at t.
    **/
   RealArray1D basisFunctionValues(Real t, int m)
   { return basis->values(t,m); }

   /** Designation of basis functions. **/
   string basisName(){ return basis->name(); }

   /** "empirical" or "Gaussian". **/
   string regressionType();

   /** The lower triangular Cholesky root of the kernel matrix 
    *  \f$(K(s_i,s_j))\f$.
    **/
   LTRRealMatrix& getR(){ return R; }
	
   /** The regression type is set to GAUSSIAN and the
    *  prior mean to the origin. Use public methods to change this
    *  setup. 
    *
    *  @param dmax maximum degree N of expansion.
    *  @param t data abscissas \f$t_j,\ 0<=j<=n\f$.
	 *  @param w function data \f$w_j=f(t_j),\ 0<=j<=n\f$.
    *  @param bFcns object containing the basis functions.
	 **/
	GPR(int dmax, RealArray1D& t, RealArray1D& w, BasisFunctions* bFcns);
	
	/** The expansion \f$f_q(t)\f$ using Gaussian or empirical coefficients \f$a_k\f$
	 *  depending on the current state of <code>this</code>.
	 **/
	Real expansion(Real t, int q);
	
	/** Writes data files for the expansions expansions \f$f_0,f_1,\dots,f_q\f$ of
	 *  f in Legendre polynomials on [-1,+1] with coefficients computed by
	 *  Gaussian regression for potting with gnuplot.
	 *  
	 *  The expansions are evaluated at 801 evenly spaced points \f$t_j\in[-1+1]\f$
	 *  and written to a file "ExpansionData.txt" in gnuplot data format.
    *  The first column contains the points \f$t_j\f$, the next column the true
	 *  function values \f$f(t_j)\f$ and the subsequent columns the expansions
	 *  \f$f_0(t_j),f_1(t_j),\dots,f_q(t_j)\f$.
	 *
	 *  See "doc/gnuplot_Readme.html" for instructions how to plot such data
	 *  with gnuplot.
	 *
    *  The data points are written to the file "FunctionData.txt" to be overlaid
	 *  with the graph of the function f in the plot. The first columns
	 *  contains the points \f$s_j\f$ and the second column the data points \f$y_j\f$.
	 *
	 *  Note: the parameter f must be the function that set up the array
	 * <code>this::y</code>, that is, we must have \f$y_j=f(s_j)\f$.
    *
    * @param f pointer to pointer to function to be expanded.
    **/	
   void expansionData(RealFunction f, int q);

   /** Writes the data files (see {@link expansionData(RealFunction,int)})
	 *  for plotting the regressors \f$f_0,f_1,\dots,f_N\f$ of the function \f$f(t)\f$
    *  with gnuplot. User chooses from several sample functions.
    *  The mean of the Gaussian prior P can be set to the origin, the vector of
    *  empirical coefficients or the vector (r,r,...,r), with r a real number.
    *  This last option is included to investigate the dependence of the
    *  algorithm on the mean of P.
	 *
	 *  The coefficients are computed from the data points \f$y_j=f(s_j)\f$,
	 *  \f$0\leq j\leq n\f$. The points \f$s_j\f$ can be evenly spaced in [-1,+1]
	 *  with \f$s_0=-1\f$, \f$s_n=+1\f$ or they can be uniformly random.
	 *  The function values \f$y_j=f(s_j)\f$ can be exact or corrupted by
	 *  independent Gaussian noise. Noisy data are of the form
	 *  \f[y_j=f(s_j)+\sigma Z_j,\f]
	 *  where the \f$Z_j\f$ are independent and standard normal. That is the standard 
	 *  deviation of the error is independent of the size of \f$f(s_j)\f$.
	 *  The user supplies all parameters.
    **/
   static void regressionExample();

   /** Tests the basis functions \f$\psi_0,\dots,\psi_q\f$ for orthonormality
    *  using Monte Carlo integration on m random points in [-1,+1].
	 *  Prints the matrix of inner products \f$(\psi_i,\psi_j)\f$.
    **/
   void orthoTest(int q, int m);

   /** Prints the values \f$\psi_i(t_j)\f$ of the basis functions
    *  \f$\psi_i\f$ on [-1,+1] for \f$i=0,1,\dots,q\f$ at evenly spaced points
    *  \f$t_j\in[-1,+1]\f$, \f$j=0,1,...,m\f$ to the file BasisFunctionsData.txt in
	 *  gnuplot data format.
	 *
	 *  The first column contains the points \f$t_j\f$ and the next columns
	 *  contain \f$\psi_0(t_j),\psi_1(t_j),\dots,\psi_q(t_j)\f$. See
	 *  "doc/gnuplot_Readme.html" for instructions to plot these data
	 *  with gnuplot.
    **/
   void basisFunctionData(int q, int m);

   /** Runs {@link basisFunctionData} after user supplies parameters. **/
   static void printBasisFunctions();	


private:

   RegressionType regrType; // flag, type of regression (EMPIRICAL, GAUSSIAN).
   int N;                   // expansions cut off at psi_N
	int n;                   // points s_j, j=0,...,n.
   RealArray1D s;           // the points s_j, j<=n
   RealArray1D y;           // the values y_j, j<=n
   BasisFunctions* basis;   // the basis functions psi_k
   string basis_name;       // name of basis
	RealMatrix psi;          // basis functions psi(k,j)=psi_k(s_j), 0<=k<=N, 0<=j<=n.
	UTRRealMatrix K;         // kernel matrix K(s_i,s_j), 0<=i,j<=n.
	LTRRealMatrix R;         // lower triangular Cholesky root of K.
   RealArray1D mu;          // center of the Gaussian prior
   RealArray1D a;           // Gaussian coefficients
   RealArray1D empCoeff;    // empirical coefficients

   Real EA(int k);   // coefficient a_k=E(A_k) of psi_k
   Real EC(int k);   // empirical coefficient of psi_k
	void writeCoefficients();
   void recomputeGaussianCoefficients();
   void fullInitialization();


   		
}; // end GPR

 

#endif
