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

#include "TypedefsMacros.h"
#include "RandomObject.h"                 // base class
#include "Matrix.h"                       // typedef problem in forward declarations


MTGL_BEGIN_NAMESPACE(Martingale)


// dependencies: 
// class RealVector;



/**<p>Class improving the convergence of expectations of the
 * underlying random variable X by the use of <i>control variates</i>.</p>
 * 
 * @author  Michael J. Meyer 
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
	ControlledRandomVariable();

	
// RANDOM DEVIATE - CONTROL VARIATE PAIR, CONTROL VARIATE MEAN
	
    /** random deviate - control variate pair from base class RandomVector
	 *  so a derived class can define this.
	 */
    virtual RealVector nextValue() = 0;

    
    /** <p>The mean of the control variate.
     * It is either known or derived from a simulation which is significantly 
     * faster than the simulation of X (<code>this</code>).</p>
     */
    virtual Real getControlVariateMean() = 0;
	
	
	
// AUXILLIARY CLASS
	
     /** Auxilliary class: the random variable \f$X-\beta(Y-E(Y))\f$, 
	  *  see book, 2.8.
      */
     class ControlledVariable : public RandomVariable {
		
	       ControlledRandomVariable* XC;
	
	 public:
		 
	       Real nextValue();
	
	       ControlledVariable(ControlledRandomVariable* xc) : XC(xc) { }
      };
	
			
   /** <p>The random variable \f$X-\beta(Y-E(Y))\f$, where Y is the control 
     * variate of X and \f$\beta\f$ is the beta coefficient.</p>
	 */
    RandomVariable* controlled();
		

 
// MEAN
    
    
    /** <p>Expectation of X computed from a sample of size N.
	 * It's simply the ordinary Monte Carlo expectation of the controlled 
     * version <code>controlled()</code> of X.</p>
     *
     * @param N sample size.
     */
    Real expectation(int N);
    

    
    /** <p>Same as {@link #conditionalExpectation(int)} but with 
     * computational progress reported as count down.</p>
     *
     * @param N Sample size.
	 * @param message string descriptive of computation.
     */
    Real expectation(int N, string message);
    
    

// TEST OF CONTROL VARIATE MEAN FORMULA
 
 
     /** <p>Tests if the method for computing the mean of the
      *  control variate at time zero is correct by comparing the returned
      *  value against a Monte Carlo mean of the control variate.</p>
	  *
	  * @param N sample size for the Monte Carlo control variate mean 
      */
     void controlVariateMeanTest(int N);
     
     
// BETA COEFFICIENT, CORRELATION WITH CONTROL VARIATE

     
    /** <p>Computes the coefficient beta=Cov(X,Y)/Var(X), where Y is the 
     *  control variate of X (<code>this</code>); nBeta=2000 samples are 
     *  generated to estimate this quantity.</p>
     */  
    Real betaCoefficient();    
    
     
     
     /** <p>The correlation of the control variate with random variable X 
      * (<code>this</code>) computed from a sample of size N.</p> 
      *
      * <p>This routine is only called to test the quality of a prospective 
      * control variate. Thus we forgo efficiency and simply reduce this to
      * variance, covariance computations.</p>
      *
      * @param N Sample size.
      */
      Real correlationWithControlVariate(int N);
        
 
}; //end ControlledRandomVariable



MTGL_END_NAMESPACE(Martingale)

#endif
 