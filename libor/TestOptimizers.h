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

#ifndef martingale_testoptimizers_h    
#define martingale_testoptimizers_h

#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)

/** A collection of free standing short test programs.
 */
 


/** Test program in dimension n. Very nasty test function 
  *  \f[f(x)=1\bigg/\sum\nolimits_jexp(-x_j^2)\f]
  *  consists of three very narrow valleys intersecting at right angles. 
  *  The optimizer must crawl along each valley minimizing each variable separately. 
  *  Takes a huge number of steps. Minimum is assumed at the origin.
  *
  * @param n dimension.
  * @param steps maximum number steps.
  */
void testDownhillSimplex(int n, int steps);


/** <p>Test program in dimension n. Very nasty test function 
 *  \f[f(x)=1\bigg/\sum\nolimits_jexp(-x_j^2)\f]
 *  consists of n very narrow valleys intersecting at right angles. 
 *  The optimizer must crawl along each valley minimizing each variable separately. 
 *  Takes a huge number of steps. Minimum is assumed at the origin. The search starts 
 *  the point \f$x_j=1.6+j*0.6\f$.
 *
 * <p>For \f$n\geq 3\f$ this fails miserably. It seems that we have to rescale the 
 * objective function so that the gradient has coordinates of roughly the same
 * magnitude. 
 *
 * @param n dimension.
 * @param nVals maximum number of function evaluations.
 */
void testBFGS(int n, int nVals);
	

/** Test program in dimension n. Very nasty test function 
 *  \f[f(x)=1\bigg/\sum\nolimits_jexp(-x_j^2)\f]
 *  consists of three very narrow valleys intersecting at right angles. 
 *  The optimizer must crawl along each valley minimizing each variable separately. 
 *  Takes a huge number of steps. Minimum is assumed at the origin.
 *
 * @param n dimension.
 * @param nVals number of function evaluations.
 */
void testSobolSearch(int n, int nVals);

	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 