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

#ifndef martingale_stochasticgenerator_h
#define martingale_stochasticgenerator_h

#include "TypedefsMacros.h"
#include "Random.h"
#include "Matrices.h"
#include <iostream>

MTGL_BEGIN_NAMESPACE(Martingale) 




// we are using
class LowDiscrepancySequence;




/*! \file StochasticGenerator.h
 *  <p>A stochastic generator is a driver generating independent standard normal
 *  deviates and writing these deviates to vectors or matrices. The deviates
 *  are used to compute the next observation of a random quantity \f$X\f$.
 *  If \f$X\f$ is a deterministic function of the path of a stochastic process,
 *  observations of \f$X\f$ amount to path computation in which case the 
 *  stochastic generator "drives" the path of the process.</p>
 *
 * <p>The effective dimension \f$d\f$ of the simulation is the number of deviates needed
 * to compute a single observation of the random quantity. If as usual the deviates 
 * are obtained by application of the inverse cumulative normal distribution to 
 * independent uniform deviates in \f$(0,1)\f$, these uniform deviates must be
 * equidistributed in \f$(0,1)^d\f$ for Monte Carlo methods to be theoretically sound.</p>
 *
 * <p>If the dimension is high random number generators might not (do not) have this
 * property. In this case one is well advised to rely on low discrepancy sequences
 * which have this property in all dimensions.</p>
 *
 * <p>This file declares the interface for a stochastic driver and some drivers
 * which we need for our stochastic processes. We have drivers generating
 * pseudorandom numbers based on the Mersenne Twister and drivers generating 
 * quasirandom numbers based on the Sobol sequence.</p>
 */


/*********************************************************************************
 *
 *         Stochastic Generator (Driver)       
 *
 *********************************************************************************/

/** <p>Interface. Provides routines to fill a vector, respectively 
 *  matrix with standard normal deviates needed to drive a simulation.
 *  Empty default implementations provided to admit selective redefinition.</p>
 *
 * <p>The purpose of this is to switch between pseudo random dynamics based on a
 * random number generator (MC) and quasi random dynamics based on a low discrepancy 
 * sequence (QMC) by simply assigning a corresponding generator object.</p>
 *
 * @author Michael J. Meyer
 */
class StochasticGenerator {
	
protected:
	
	int n;  // dimension of the process, 
	        // size of the vector needed to drive one time step
	
public:
	
	/** @param m size of vector needed to drive one time step,
	 *  default = 1.
	 */
	StochasticGenerator(int m=1) : n(m) {  }
	
// ADAPTERS
	
	/** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the vector Z.
	 *  Z must have zero based index.
	 */
	virtual void newWienerIncrements(int t, int T, RealVector& Z) const {  }
	
	/** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the vector Z.
	 *  Z must have zero based index.
	 */
	virtual void newWienerIncrements(int t, int T, Real* Z) const {  }
    
	/** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the upper triangular 
	 *  matrix Z. Z must have zero based indices.
	 */
    virtual void newWienerIncrements(int t, int T, UTRRealMatrix& Z) const  {  }
	
	/** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the square matrix Z.
	 *  Z must have zero based indices.
	 */
    virtual void newWienerIncrements(int t, int T, RealMatrix& Z) const {  }
	
	/** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the square matrix Z.
	 *  Z must have zero based indices.
	 */
    virtual void newWienerIncrements(int t, int T, RealArray2D& Z) const {  }
	
    /** Writes standard normal deviates needed to drive a simulation 
	 *  from discrete time t to discrete time T into the matrix Z.
	 *  Z must have zero based index.
	 */
	virtual void newWienerIncrements(int t, int T, Real** Z) const  {  }
	
	/** Restarts the generator. Necessary for low discrepancy sequences.
	 *  Default implementation: empty, suitable for random number generators.
	 */
	virtual void restart() const  {  }
	
	/** String identifying the generator.
	 */
   /** Message and fields.*/
   virtual std::ostream& printSelf(std::ostream& os) const;
	
}; // end StochasticGenerator



// GLOBAL INSERTION
std::ostream& operator << (std::ostream& os, const StochasticGenerator& sg);


/*********************************************************************************
 *
 *         Libor Drivers   
 *
 *********************************************************************************/


/** Stochastic generator for a LiborProcess based on the Mersenne Twister.
 */
class MonteCarloLiborDriver : public StochasticGenerator {
	
public:
	
	MonteCarloLiborDriver(int n) : StochasticGenerator(n) {  }
	
	/** Writes standard normal deviates needed to drive a Libor path 
	 *  simulation from discrete time t to discrete time T into the upper
	 *  triangular matrix Z. Z must have zero based indices.
	 */
    void newWienerIncrements(int t, int T, UTRRealMatrix& Z);
	
	/** String identifying the generator.
	 */
   std::ostream& printSelf(std::ostream& os) const;


}; // end MonteCarloLiborDriver



/** Stochastic generator for a LiborProcess based on the Sobol sequence.
 */
class SobolLiborDriver : public StochasticGenerator {
	
	LowDiscrepancySequence* lds;
	
public:
	
	SobolLiborDriver(int n) : StochasticGenerator(n), lds(0) {  }
	
	/** Writes standard normal deviates needed to drive a Libor path 
	 *  simulation from discrete time t to discrete time T into the upper
	 *  triangular matrix Z. Z must have zero based indices.
	 */
    void newWienerIncrements(int t, int T, UTRRealMatrix& Z);
	
	void restart();
	
    /** String identifying the generator.
	 */
	std::ostream& printSelf(std::ostream& os) const;


}; // end SobolLiborDriver



/*********************************************************************************
 *
 *         VectorProcess Drivers       
 *
 *********************************************************************************/


/** Stochastic generator for a {@link VectorProcess} based on the Mersenne Twister.
 */
class MonteCarloVectorDriver : public StochasticGenerator {
	
public:
	
	MonteCarloVectorDriver(int n) : StochasticGenerator(n) {  }
	
	/** Writes standard normal deviates needed to drive one 
	 *  path from discrete time t to discrete time s into the 
	 *  matrix Z. 
	 */
    void newWienerIncrements(int t, int s, Real** Z);
	
	/** Writes standard normal deviates needed to drive one path 
	 *  path from discrete time t to discrete time s into the 
	 *  matrix Z. 
	 */
	void newWienerIncrements(int t, int s, RealMatrix& Z);
	
	/** String identifying the generator.
	 */
	std::ostream& printSelf(std::ostream& os) const;
	

}; // end MonteCarloVectorDriver



/** Stochastic generator for a {@link VectorProcess} based on the Sobol sequence.
 */
class SobolVectorDriver : public StochasticGenerator {
	
	int T;     // number of time steps to horizon
	LowDiscrepancySequence* lds;
	
public:
	
	/** @param n size of Z-vector needed to drive one time step.
	 *  @param T_oo number of time steps to horizon.
	 */
	SobolVectorDriver(int n, int T_oo) : StochasticGenerator(n), 
	T(T_oo), lds(0) {  }
	
	/** Writes standard normal deviates needed to drive one path 
	 *  path from discrete time t to discrete time s into the 
	 *  matrix Z. 
	 */
    void newWienerIncrements(int t, int s, Real** Z);
	
    
	/** Writes standard normal deviates needed to drive one path 
	 *  path from discrete time t to discrete time s into the 
	 *  matrix Z. 
	 */
	void newWienerIncrements(int t, int s, RealMatrix& Z);
	
	void restart();
	
    /** String identifying the generator.
	 */
	std::ostream& printSelf(std::ostream& os) const;

}; // end SobolVectorDriver



/*********************************************************************************
 *
 *         ScalarProcess Drivers    
 *
 *********************************************************************************/


/** Stochastic generator for a {@link VectorProcess} based on the Mersenne Twister.
 */
class MonteCarloScalarDriver : public StochasticGenerator {
	
public:
	
	MonteCarloScalarDriver() : StochasticGenerator() {  }
	
	/** Writes standard normal deviates needed to drive one path 
	 *  path from discrete time t to discrete time s into the vector Z. 
	 */
    void newWienerIncrements(int t, int s, Real* Z);
	
	/** String identifying the generator.
	 */
	std::ostream& printSelf(std::ostream& os) const;

}; // end MonteCarloScalarDriver



/** Stochastic generator for a {@link ScalarProcess} based on the Sobol sequence.
 */
class SobolScalarDriver : public StochasticGenerator {
	
	int T;   // number of time steps to horizon
	LowDiscrepancySequence* lds;
	
public:
	
	/** @param T_oo number of time steps to horizon. */
	SobolScalarDriver(int T_oo) : StochasticGenerator(), T(T_oo), lds(0) {  }
	
	/** Writes standard normal deviates needed to drive one path 
	 *  path from discrete time t to discrete time s into the vector Z. 
	 */
    void newWienerIncrements(int t, int s, Real* Z);
	
	void restart();
	
    /** String identifying the generator.
	 */
	std::ostream& printSelf(std::ostream& os) const;

}; // end SobolScalarDriver






MTGL_END_NAMESPACE(Martingale) 


#endif

