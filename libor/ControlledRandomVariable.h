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


// the entire implementation is in the header


#ifndef martingale_controlled_randomvariable_h    
#define martingale_controlled_randomvariable_h 


#include "RandomObject.h"

MTGL_BEGIN_NAMESPACE(Martingale)

/**<p>Class improving the convergence of expectations of the
 * underlying random variable X by the use of <i>control variates</i>.</p>
 *
 * 
 * 
 * @author  Michael J. Meyer
 * 
 */
class ControlledRandomVariable : private RandomVector {
    
    // number of samples used to compute the {@link #betaCoefficient}.
    static const int nBeta=2000;
	Real beta;
    
public:   
	

// ACCESSORS

	Real getBeta() const { return beta; }
	
	/** This function MUST be called from the concrete subclass defining
	 *  <code>nextValue()</code>. The beta coefficient cannot be initialized
	 *  by the current constructor since in a constructor the virtual
	 *  function call <code>nextValue()</code> resolves to the local 
	 *  (here: pure virtual) version.
	 */
	void setBeta() { beta=betaCoefficient(); }
	
	
// CONSTRUCTOR: WE HAVE A DEFAULT CONSTRUCTOR
	
	/** <code>beta</code> must be initialized from each concrete subclass 
	 *  by calling <code>setBeta()</code>.
	 */
	ControlledRandomVariable() : RandomVector(2), beta(0.0) { }

	
// RANDOM DEVIATE - CONTROL VARIATE PAIR, CONTROL VARIATE MEAN
	
    /** random deviate - control variate pair from base class RandomVector
	 *  so a derived class can define this.
	 */
    virtual vector<Real> nextValue() = 0;

    
    /** <p>The mean of the control variate.
     * It is either known or derived from a simulation which is significantly 
     * faster than the simulation of X (<code>this</code>).</p>
     */
    virtual Real getControlVariateMean() = 0;
	
	
     /** Auxilliary class: the random variable \f$X-\beta(Y-E(Y))\f$, 
	  *  see book, 2.8.
      */
     class ControlledVariable : public RandomVariable {
		
	     ControlledRandomVariable* XC;
	
	 public:
		 
	      Real nextValue() 
	      {
               Real mean_y=XC->getControlVariateMean(),     // control variate mean E(Y)
		            beta=XC->getBeta();                     // beta coefficient
               vector<Real> v=XC->nextValue();
               Real  x=v[0], y=v[1], xc=x-beta*(y-mean_y);
   	           return xc;
           } //end nextValue
	
	       ControlledVariable(ControlledRandomVariable* xc) : XC(xc) { }
      };
	
			
   /** <p>The random variable \f$X-\beta(Y-E(Y))\f$, where Y is the control 
     * variate of X and \f$\beta\f$ is the beta coefficient.</p>
	 */
    RandomVariable* controlled()
	{ 		 
		 return new ControlledVariable(this);
    }
		

 
              
/*******************************************************************************    
    
    MEAN 
    Fixed sample size, all samples independent
 
*******************************************************************************/
    
    
    /** <p>Expectation of X computed from a sample of size N.
	 * It's simply the ordinary Monte Carlo expectation of the controlled 
     * version <code>controlled()</code> of X.</p>
     *
     * @param N sample size.
     */
    Real expectation(int N)
    {
		return controlled()->expectation(N);
    }
    



    

    
/*******************************************************************************    
    
    MEAN, computation time reported as count down
 
*******************************************************************************/
   
 
    
    
    /** <p>Same as {@link #conditionalExpectation(int)} but with 
     * computational progress reported as count down.</p>
     *
     * @param N Sample size.
	 * @param message string descriptive of computation.
     */
    Real expectation(int N, string message)
    {        
		return controlled()->expectation(N,message);
    }
    
    

    

		

     
/*******************************************************************************

     TEST OF CONTROL VARIATE MEAN FORMULA
     
*******************************************************************************/
 
     /** <p>Tests if the method for computing the mean of the
      *  control variate at time zero is correct by comparing the returned
      *  value against a Monte Carlo mean of the control variate.</p>
	  *
	  * @param N sample size for the Monte Carlo control variate mean 
      */
     void controlVariateMeanTest(int N)
     {
         cout << "\nTesting control variate mean:\n";
         
         // analytic control variate mean
         cout << "analytic: " << getControlVariateMean() << endl;
         
         // Monte Carlo control variate mean
         cout << "Monte Carlo: " << RandomVector::expectation(N)[1] << endl;
         
     } // end controlVariateMeanTest
     
     
     
    

/*******************************************************************************

     BETA COEFFICIENT, CORRELATION WITH CONTROL VARIATE
     
*******************************************************************************/
     
     
    /** <p>Computes the coefficient beta=Cov(X,Y)/Var(X), where Y is the 
     *  control variate of X (<code>this</code>); nBeta=2000 samples are 
     *  generated to estimate this quantity.</p>
     */  
    Real betaCoefficient()
    {
        // Recall that Cov(X,Y)=E(XY)-E(X)E(Y) and Var(X)=E(X^2)-E(X)^2.
		int N=nBeta;
		Real sum_X=0, sum_Y=0,
             sum_XX=0, sum_XY=0;
        
        for(int n=0;n<N;n++){
			
            vector<Real> v=nextValue();
            Real x=v[0],y=v[1];
            
            sum_X+=x; sum_Y+=y;
            sum_XX+=x*x; sum_XY+=x*y;
        }
        
        return (N*sum_XY-sum_X*sum_Y)/(N*sum_XX-sum_X*sum_X);
        
    } //end betaCoefficient
    
    
     
     
     /** <p>The correlation of the control variate with random variable X 
      * (<code>this</code>) computed from a sample of size N.</p> 
      *
      * <p>This routine is only called to test the quality of a prospective 
      * control variate. Thus we forgo efficiency and simply reduce this to
      * variance, covariance computations.</p>
      *
      * @param N Sample size.
      */
      Real correlationWithControlVariate(int N)
      {
          return correlation(0,1,N);
         
      } //end correlationWithControlVariate
        
 
}; //end ControlledRandomVariable



MTGL_END_NAMESPACE(Martingale)

#endif
 