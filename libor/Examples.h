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

/** Timing of N multiplications of square matrices in dimension dim.
 *  Compares row dot column (suboptimal) to row dot row (optimal) memory
 *  access pattern.
 */
void timeMatrixMultiply(int dim, int N);


/*******************************************************************************    
    
              MATRIX EXPONENTIALS exp(At), exp(-At)  
	
*******************************************************************************/

/**<p>We compute the matrix exponentials <code>H(t)=exp(At), K(t)=exp(-At)</code>
 * for <code>t=1</code>, where A is a (constant) upper triangular dim by dim matrix 
 * and time t is a real number.
 *
 * The computation simulates the paths <code>t -> H(t),K(t)</code> until time
 * <code>t=1</code> using the dynamics 
 * \f[dH(t)=H(t)Adt,\quad dK(t)=-K(t)Adt.\f]
 * The purpose is to determine what step sizes are necessary to obtain reasonable 
 * accuracy and how long such a computation takes in various dimensions.
 *
 * <p>This path dynamics is discretized as
 * \f[H(t+dt)-H(t)=H(t)Adt,\quad\hbox{ie.}\quad H(t+dt)=H(t)(I+Adt).\f]
 * \f[K(t+dt)-K(t)=-K(t)Adt,\quad\hbox{ie.}\quad K(t+dt)=K(t)(I-Adt).\f]
 * The accuracy is checked by checking the relationship 
 * <code>K(t)=H(t)^{-1}</code>, that is, <code>H(t)K(t)=K(t)H(t)=I</code>.
 * Matrix dimension <code>dim</code> and step size <code>dt</code> are user supplied.
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
     *  in Array.h.
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
 