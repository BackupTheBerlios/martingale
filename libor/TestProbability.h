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

#ifndef martingale_testprobability_h    
#define martingale_testprobability_h


#include "TypedefsMacros.h"

MTGL_BEGIN_NAMESPACE(Martingale)

/** A collection of free standing short test programs.
 */
 

/*******************************************************************************    
    
                      RANDOM NUMBER GENERATOR       
	
*******************************************************************************/



// STANDARD NORMAL COVARIANCE MATRIX    

/** Times the computation of a standard normal covariance matrix.
 *  Then prints the largest and smallest diagonal and off diagonal elements.
 *  This is a test of the random number generation. The true matrix is the
 *  identity matrix.
 */
void testSTNCovarianceMatrix(int d, int N);




/*******************************************************************************    
    
                             RANDOM OBJECTS  
	
*******************************************************************************/


/** Defines the path functional \f$H=X_1(T)+X_2(T)\f$ of a two dimensional
 *  Brownian motion X conditioned on the state at time t. The conditional distribution 
 *  is normal with mean X_1(t)+X_2(t) and variance 2*(T-t)*dt.
 *  The conditional mean and variance are computed from a sample of size 
 *  nPath and compared with the analytic values.
 */
void testPathFunctional(int t, int T, int nPath);
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 