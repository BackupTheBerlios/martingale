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




#ifndef martingale_randomvariables_h    
#define martingale_randomvariables_h 

#include "TypedefsMacros.h"
#include "RandomObject.h"                     // base class
#include "ControlledRandomVariable.h"       // base class


MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file RandomVariables.h
 *  A collection of some random variables.
 */


/*******************************************************************************    
    
                        SOME CONCRETE RANDOMVECTORS 
 
*******************************************************************************/

/** Standard normal randomvector.
 */
class StandardNormalVector : public RandomVector {
	
public:
	
	StandardNormalVector(int d) : RandomVector(d) {}
	
	virtual ~StandardNormalVector(){}
	
	/** next observation. */
	RealVector nextValue();

}; // end StandardNormalVector
	


 

/*******************************************************************************    
    
                      Standard Normal Variable             
	
*******************************************************************************/

/** Test class with maningless control variate 0, used only to test 
 *  the class <code>ControlledRandomVariable</code>.
 */ 
class StandardNormalVariable : public ControlledRandomVariable {
	
public:
	
	// sample - control variate pair
	RealVector nextValue();
	
	Real getControlVariateMean() { return 0; }
	
	StandardNormalVariable() { setBeta(); }

}; // end StandardNormalVariable





/*******************************************************************************    
    
                      Empirical Random Variable             
	
*******************************************************************************/


/** <p>A random variable X distributed according to the empirical distribution 
 *  associated with a data sample. This distribution assumes that the data points 
 *  exhaust all possible values and all values are equally likely. Drawing samples 
 *  of X amounts to sampling from the data set with replacement.</p>
 *
 * <p>This is the type of random variable to use if you have a raw data set.
 * Samples drawn from the empirical distribution do not need to have the same size as 
 * the data set and indeed can be arbitrarily large (sampling with replacement).</p>
 *
 * <p>Class maintains a reference to the data set array but does not construct
 * a copy of this array. The data are assumed to be of scalar type <code>Real</code>. 
 * </p>
 *
 * @author  Michael J. Meyer
 */
class EmpiricalRandomVariable : public RandomVariable {
    
    int sampleSize;                          // number of data points
	Real* dataSet;                           // array of data points
	
public:
    
    /** Number of data points.
     */
    int getSampleSize() const { return sampleSize; }
    
    /** The array containing the sample data.
     */
    Real* getDataSet() const { return dataSet; }
    
    
    /** Constructor.
     *
	 * @param sample_size number of data points
     * @param data_set array containing the data sample.
     */
    EmpiricalRandomVariable(Real* data_set, int sample_size) :
	sampleSize(sample_size), 
	dataSet(data_set)
    {   }
     
    
    /** <p>Sampling from the distribution of X, 
     *  Samples from the underlying data set with replacement.</p>
     */
    Real nextValue();
	
	
// TEST PROGRAM
	
	/** Allocates the data set 0,1,...,n-1, then prints a sample of size 200  
	 * from the empirical distribution.
	 */
	static void test(int n);
   
    
}; // end EmpiricalRandomvariable




MTGL_END_NAMESPACE(Martingale)


#endif