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


#ifndef martingale_examples_h    
#define martingale_examples_h

#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Examples)

/** A collection of short example and test programs.
 *  Each member function <code>foo</code> represents an example
 *  and is intended to be called in <code>main</code>.
 *  The function definitions are included in the header as well
 *  since they are mostly fairly short.
 */
 
/*******************************************************************************    
    
             TIMING THE GENERATION OF STANDARD NORMAL DEVIATES
			 
	
*******************************************************************************/

   
/** Times (countdown) the generation of 10 million standard normal increments
 *  using the class <code>LoopStatus</code> 
 */ 
void timeSTN();


/*******************************************************************************    
    
              TIMING UPPER TRIANGULAR MATRIX TO MATRIX MULTIPLY   
	
*******************************************************************************/

/** Times N multiplications of an square matrix by an upper triangular matrix in 
 *  dimension dim.
 */
void timeMatrixMultiply(int dim, int N);


/*******************************************************************************    
    
              MATRIX EXPONENTIALS exp(At), exp(-At)  
	
*******************************************************************************/

/**<p>We compute the matrix exponentials <code>H(t)=exp(At), K(t)=exp(-At)</code>
 * for <code>t=1</code>, where A is a (constant) upper triangular dim by dim matrix 
 * and time t is a real number.</p>
 *
 * The computation simulates the paths <code>t -> H(t),K(t)</code> until time
 * <code>t=1</code> using the dynamics <code>dH(t)=H(t)Adt, dK(t)=-K(t)Adt</code>. 
 * The purpose is to determine what step sizes are necessary to obtain reasonable 
 * accuracy and how long such a computation takes in various dimensions.</p>
 *
 * <p>This path dynamics is discretized as
 * <center>H(t+dt)-H(t)=H(t)Adt, ie. H(t+dt)=H(t)(I+Adt)</center>
 * <center>K(t+dt)-K(t)=-K(t)Adt, ie. K(t+dt)=K(t)(I-Adt)</center>
 * and computes the exponential <code>exp(At)</code> as
 * <code>(1+Adt)^{t/dt}=(1+At/m)^m</code> where <code>m=t/dt</code> in 
 * accordance with the limit formula <code>(1+x/m)^m -> e^x, m->oo</code>.
 * So if we want the exponentials <code>exp(At)</code> for only one value
 * of <code>t</code> then there are faster methods to do this (let m be a 
 * power of 2 and use repeated squaring). However in the simulation of 
 * forward Libors we actually have to carry out a path computation as above.
 * </p>
 *
 * <p>The accuracy is checked by checking the relationship 
 * <code>K(t)=H(t)^{-1}</code>, that is, <code>H(t)K(t)=K(t)H(t)=I</code>.</p> 
 *
 * @param dim matrix dimension
 * @param dt step size in the matrix path dynamics
 */
void matrixExponentials();


             
/*******************************************************************************
 *
 *                    LMM PATH TIMING
 *
 ******************************************************************************/


// PATH TIMING
    
    /** Allocates sample Libor process in dimension n
     *  then times the computation of N full Libor paths.
     *  Note the dramatic speedup if SUBSCRIPT_CHECK is undefined
     *  in Array.h and Matrices.h.
     */   
    void liborPathTiming(int n, int N);	
	
            
/*******************************************************************************
 *
 *                    STOCHASTIC PROCESSES
 *
 ******************************************************************************/
	
	
	/** Allocates a Brownian motion X in dimension d and computes the mean time
	 *  until X hits the boundary of a ball centered at the origin with radius 1.
	 *
	 * @param dim dimension.
	 * @param T time steps to horizon 
	 * (must be large enough to ensure hit on the boundary).
	 * @param dt size of time step.
	 * @param N number of paths.
	 */
    void brownianMotionInBall(int dim, int T, Real dt, int N);	



MTGL_END_NAMESPACE(Examples)
MTGL_END_NAMESPACE(Martingale)

#endif
 