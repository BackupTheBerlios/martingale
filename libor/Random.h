/* WARANTY NOTICE AND COPYRIGHTThis program is free software; you can redistribute it and/ormodify it under the terms of the GNU General Public Licenseas published by the Free Software Foundation; either version 2of the License, or (at your option) any later version.This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty ofMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See theGNU General Public License for more details.You should have received a copy of the GNU General Public Licensealong with this program; if not, write to the Free SoftwareFoundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.Copyright (C) Michael J. Meyermatmjm@mindspring.comspyqqqdia@yahoo.com*/

#ifndef martingale_random_h    
#define martingale_random_h

#include "TypedefsMacros.h"


/*! \files Random.h
 *  Random number generators: uniform and standard normal.
 */


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Random)

    
//typedef boost::mt19937 MTW;     
//typedef boost::uniform_01<MTW,Real> U_01_RNG;
                           
//extern U_01_RNG uniform01;             // uniform rng in [0,1) based on Mersenne twister mt19937
   
/**  generates a random draw from {1,-1} */
int sign();

/**  cumulative normal distribution function */
Real N(Real x);
	
/**  inverse normal cumulative distribution function 
 */
Real nInverse(Real x);
 
/** simple uniform random number generator */
Real U01();

/** simple uniform random integers in [0,n-1) */
int Uint(int n);
   
/** standard normal deviate based on Mersenne twister mt19937 and
 *  nInverse().
 */
Real sTN();
   

MTGL_END_NAMESPACE(Random)   
MTGL_END_NAMESPACE(Martingale)

#endif
 
