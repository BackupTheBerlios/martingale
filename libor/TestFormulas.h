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

#ifndef martingale_testformulas_h    
#define martingale_testformulas_h


#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)


/*******************************************************************************
 *
 *               FORMULAS
 *
 ******************************************************************************/

/** Monte Carlo test of the analytic expectations 
 *  book::Appendix::C.12, C.13, C.16.
 */
void testExponentialIntegralFormulas();
	


/** Monte Carlo test of the volatility integrals for the samples of the 
 *  2 volatility surface types JR and M.
 *
 * @param N number of Monte Carlo sample points.
 * @param precision maximal acceptable relative error in percent.
 */
void testVolSurfaceIntegrals(int N, Real precision);
	
	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 