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

#ifndef martingale_randomobject_h    
#define martingale_randomobject_h

#include "TypedefsMacros.h"
#include "Matrices.h"
#include "Utils.h"
#include "Random.h"  
#include <string>


// The whole implementation is in the header!

/** <p><code>RangeType</code>-valued random variable <code>X</code>. 
 *  <code>RangeType</code> is assumed to be a vector type with components of 
 *  type <code>ScalarType</code>. These are accessed via the subscripting operator 
 *  <code>[](int i)</code> and defined for \f$0\leq i<dim\f$ where 
 *  <code>dim</code> is the dimension of <code>RangeType</code>.
 *  Both types default to the type Real and the dimension <code>dim</code>
 *  defaults to 1.</p>
 *
 *  <p>The subscripting operators only come into play if covariances are computed.
 *  Need not be defined otherwise. The dimension of <code>RangeType<code> only comes into 
 *  play if the covariance matrix is computed. If no covariances are computed the ScalarType 
 *  does not come into play. In each case the defaults can be used.</p>
 *
 * <p><a name="expectation"><b>Expectation:</b></a> 
 *  For sample means to be computed the type <code>RangeType</code> must support 
 *  addition and division by integers in the form of overloaded operators 
 *  <code>+=(const RangeType& x</code> and <code>/=(int N)</code> (or some type that 
 *  <code>int N</code> can be converted to). This is the minimum. More 
 *  structure in <code>RangeType</code> is required for variances, covariances and 
 *  covariance matrices.
 * </p>
 *
 *  <p><a name="variance"><b>Variance:</b></a> 
 *  Computed as \f$Var(X)=E(X^2)-E(X)^2\f$ this is defined only if 
 *  the type <code>RangeType</code> supports a multiplication operator * in the 
 *  form of an overloaded <code>operator *=(const RangeType& x)</code> and an
 *  operator <code>-=(const RangeType& x)</code>. Typically this
 *  will be true in the scalar case only.
 *  </p>
 *  
 *  <p>NOTE: the denominator "N" is used instead of the more common "N-1".
 *  The reason is a slight simplification in the code. The corresponding 
 *  estimator for the variance is biased. However we have applications in 
 *  mind for which N is 1000 or greater in which case the bias is 0.1% or
 *  smaller. This applies to covariances as well.<br>
 *  Standard deviations are not implemented  because of the absence of square 
 *  roots in general.</p>
 *
 *  <p><a name="covariance"><b>Covariance:</b></a> 
 *  the covariance of the components computed as 
 *  \f[Cov(X_i,X_j)=E(X_iX_j)-E(X_i)E(X_j)\f]
 *  The type <code>ScalarType</code> must support operators
 *  <code>+=,-=,*=(const ScalarType& u), /=(int N)</code>.
 *  </p>
 *
 *  <p><a name="covariance-matrix"><b>Covariance matrix:</b></a> 
 *  the (symmetric) square matrix of all covariances 
 *  \f[Cov(X_i,X_j), 0\leq i,j<dim.\f] 
 *  Returned as an upper triangular <code>UTRMatrix</code>.
 *  </p>
 *
 * @author Michael J. Meyer
 */



MTGL_BEGIN_NAMESPACE(Martingale)	



template<typename RangeType=Real, typename ScalarType=Real>
class RandomObject {


    int dim;        // the range dimension

	
 public:
	
    int getDimension() const { return dim; }


/*******************************************************************************    
    
                           CONSTRUCTOR
 
*******************************************************************************/      


  /** Constructor
   *
   * @param d dimension of the range type, default = 1.
   */
  RandomObject(int d=1): dim(d) { }
	  
 
    
/*******************************************************************************    
    
                           SAMPLING 
 
*******************************************************************************/    
    
    /**<p>The next observation from the distribution of <code>this</code>.
     * Returning an object of type <code>RangeType</code> involves copying and thus
	 * carries some overhead. But all the alternativs are very akward in case
	 * <code>RangeType</code> is a built in numeric type. Quite often the computation
	 * of the next observation is vastly more expensive than copying the return
	 * value.
     *
     * <p>This is the crucial method defining the random object. Typically 
     * this method will be called tens of thousands of times.
     *
     */
    virtual RangeType nextValue() = 0;
    



    

/*******************************************************************************    
    
    MEAN 
    Fixed sample size.
 
*******************************************************************************/
    
    
    /** Monte Carlo expectation computed from sample of size N.
	 *
     * @param N Sample size.
     */
    RangeType expectation(int N)
    {
		// initialization n=0 (RangeType need not have a zero element)
        RangeType EX=nextValue();
        for(int n=1;n<N;n++) EX+=nextValue(); EX/=N;
			
        return EX;

    } // end conditionalExpectation

    
 
    
    /** <p>Expectation computed from sample of size N.
     * Progress of the computation is reported to the console.</p>
     * 
     * @param N Sample size (must be a multiple of 100).
     * @param message message indicating what computation is in progress.
     */
    RangeType expectation(int N, string message)
    {
        int m=N/100;
        LoopStatus LS(message);
        RangeType EX=nextValue();
        for(int n=1;n<N;n++){

			if(n%m==0)LS.consoleReport(n,N);
            EX+=nextValue();
        }
        
        EX/=N;
        return EX;

    } // end conditionalExpectation

 

/*******************************************************************************    
    
    MEAN AND STANDARD DEVIATION
    Fixed sample size, all samples independent
 
*******************************************************************************/


    
    /** <p>Mean (return[0]) and <a href="#variance">variance</a> 
     *  (return[1]) computed from sample of size N. Default constructor for 
	 *  <code>RangeType</code> is used for returning.</p>
     *
     * @param N sample size
     */
    RangeType* meanAndVariance(int N)
    {
        RangeType x=nextValue(),
                  EX=x,            // X_1+X_2+...+X_n
                  EXX=x*x,         // X_1^2+X_2^2+...+X_n^2
	              Var_X;

        for(int n=1;n<N;n++){
            
            x=nextValue();
            EX+=x;
            EXX+=x*x;
        }
        
        EX/=N;
        EXX/=N;
		Var_X=EXX;
        Var_X-=EX*EX;
        
        RangeType* v = new RangeType[2]; v[0]=EX; v[1]=Var_X;
        return v;

    } // end meanAndVariance
	
	
	/** <p>Mean (return[0]) and <a href="#variance">variance</a> 
     *  (return[1]) computed from sample of size N. Default constructor for 
	 *  <code>RangeType</code> is used for returning.
	 *  Progress of the computation is reported to the console.</p>
     * 
     * @param N Sample size (must be a multiple of 100).
     * @param message message indicating what computation is in progress.
     */
    RangeType* meanAndVariance(int N, string message)
    {
        RangeType x=nextValue(),
                  EX=x,            // X_1+X_2+...+X_n
                  EXX=x*x,         // X_1^2+X_2^2+...+X_n^2
	              Var_X;
		
		int m=N/100;
		LoopStatus LS(message);
        for(int n=1;n<N;n++){
            
            if(n%m==0)LS.consoleReport(n,N);
			x=nextValue();
            EX+=x;
            EXX+=x*x;
        }
        
        EX/=N;
        EXX/=N;
		Var_X=EXX;
        Var_X-=EX*EX;
        
        RangeType* v = new RangeType[2]; v[0]=EX; v[1]=Var_X;
        return v;

    } // end meanAndVariance
    
    


/*******************************************************************************    
    
              VARIANCE, COVARIANCE and CORRELATION
 
*******************************************************************************/


    
    
   /** <p><a href="#variance">Variance</a>
	*  computed from a sample of size N.</p>
    *
    * @param N Size of sample used to estimate the variance.
    */
   RangeType variance(int N)
   {
        RangeType x=nextValue(),
                  EX=x,            // X_1+X_2+...+X_n
                  EXX=x*x;         // X_1^2+X_2^2+...+X_n^2

        for(int n=1;n<N;n++){
            
            x=nextValue();
            EX+=x;
            EXX+=x*x;
        }
        
        EX/=N;
        EXX/=N;
		EXX-=(EX*EX);
        return EXX;

   } // end Variance
    

    

   /** <p><a href="#covariance">Covariance</a> 
    *  computed from a sample of size N.</p>
    *
    * @param N Size of sample used to estimate the variance.
    */
   ScalarType covariance(int i, int j, int N)
   {
        RangeType  x=nextValue();
        ScalarType xi=x[i], xj=x[j],
	               EXi=xi, EXj=xj, EXiXj=xi*xj;         

        for(int n=1;n<N;n++){
            
	       x=nextValue(); xi=x[i]; xj=x[j];
	       EXi+=xi; EXj+=xj, EXiXj+=xi*xj;
        }
        
        EXi/=N;
        EXj/=N;
        EXiXj/=N;
		EXiXj -= EXi*EXj;
        return EXiXj;

   } // end covariance
   
   
   
   /**<p>Correlation Corr(X_i,X_j) computed from a sample of size N.</p> 
    *
    * <p>Recall that Corr(X,Y)=Cov(X,Y)/sqrt(Var(X))sqrt(Var(Y))
    * However efficieny calls for a different implementation.
    * There must be a global function <code>sqrt(const ScalarType&)</code>.
    * </p>
    *
    * @param i index of X-component.
    * @param j index of X-component.
    * @param N Size of sample used to estimate the correlation.
    */     
    ScalarType correlation(int i, int j, int N)
    {
        RangeType  x=nextValue();
		ScalarType sum_Xi=x[i], sum_Xj=x[j],
                   sum_XiXi=x[i]*x[i], sum_XjXj=x[j]*x[j], sum_XiXj=x[i]*x[j];
        
        for(int n=0;n<N;n++)
        {
            x=nextValue();
            sum_Xi+=x[i]; sum_Xj+=x[j];
            sum_XiXi+=x[i]*x[i]; sum_XjXj+=x[j]*x[j]; sum_XiXj+=x[i]*x[j];
        }
                
        //when divided by N^4: 
        ScalarType 
		Var_XiVar_Xj=((N*sum_XiXi-sum_Xi*sum_Xi)*(N*sum_XjXj-sum_Xj*sum_Xj));
        
        return (N*sum_XiXj-sum_Xi*sum_Xj)/sqrt(Var_XiVar_Xj);
        
    } // end correlation
  
    


/*******************************************************************************    
    
                      COVARIANCE MATRIX
 
*******************************************************************************/



   /** <p><a href="#covariance-matrix">Covariance matrix</a> 
    *  computed from a sample of size N.</p>
    *
    * @param N Size of sample used to estimate the variance.
    */
   UTRMatrix<ScalarType> covarianceMatrix(int N)
   {
       ScalarType E[dim];                       // E[i]=EX_i
       UTRMatrix<ScalarType> C(dim);            // C(i,j)=Cov(X_i,X_j)

       // initialization n=0 (the type S need not have a zero element)
       RangeType x=nextValue();
       for(int i=0;i<dim;i++){

	       E[i]=x[i];
           for(int j=i;j<dim;j++) C(i,j)=x[i]*x[j];
       }

       // accumulate samples for the component means, averaging later
       for(int n=1;n<N;n++){

	       x=nextValue();
           for(int i=0;i<dim;i++){

	           E[i]+=x[i];
               for(int j=i;j<dim;j++) C(i,j)+=x[i]*x[j];
           }
       }
   
       // average
       for(int i=0;i<dim;i++){

	       E[i]/=N;
           for(int j=i;j<dim;j++) C(i,j)/=N;
       }

       // now C(i,j)=E(X_iX_j) change this to C(i,j)=E(X_iX_j)-EX_i*EX_j
       for(int i=0;i<dim;i++)
       for(int j=i;j<dim;j++) C(i,j)-=E[i]*E[j];

       return C;

   } // end conditionalCovarianceMatrix

	
 
   /** <p><a href="#covariance-matrix">Covariance matrix</a> 
    *  computed from a sample of size N. Progress of the computation is 
    *  reported to the console.</p>
    *
    * @param N Size of sample used to estimate the variance.
    * @apram message string describing computation.
    */
   UTRMatrix<ScalarType> covarianceMatrix(int N, string message)
   {
       ScalarType E[dim];                  // E[i]=EX_i
       UTRMatrix<ScalarType> C(dim);       // C(i,j)=Cov(X_i,X_j)
       RangeType x=nextValue();

       LoopStatus LS(message); 

       // accumulate samples for the component means, averaging later
       // initialization n=0 (the type S need not have a zero element)
       for(int i=0;i<dim;i++){

	       E[i]=x[i];
           for(int j=i;j<dim;j++) C(i,j)=x[i]*x[j];
       }

       int m=N/100;
       for(int n=1;n<N;n++){

           x=nextValue();
	       if(n%m==0)LS.consoleReport(n,N);
           for(int i=0;i<dim;i++){

	           E[i]+=x[i];
               for(int j=i;j<dim;j++) C(i,j)+=x[i]*x[j];
           }
       }
   
       // average
       for(int i=0;i<dim;i++){

	       E[i]/=N;
           for(int j=i;j<dim;j++) C(i,j)/=N;
       }

       // now C(i,j)=E(X_iX_j) change this to C(i,j)=E(X_iX_j)-EX_i*EX_j
       for(int i=0;i<dim;i++)
       for(int j=i;j<dim;j++) C(i,j)-=E[i]*E[j];

       return C;

    } // end conditionalCovarianceMatrix




}; // end RandomObject



/*******************************************************************************    
    
                  RANDOMVECTORS and RANDOMVARIABLES
 
*******************************************************************************/



/** Default random variables (based on ScalarType Real).
 */
typedef RandomObject<> RandomVariable;
  

/** Default random vectors (based on ScalarType Real).
 */
typedef RandomObject< RealVector > RandomVector;
  







MTGL_END_NAMESPACE(Martingale)


#endif	
