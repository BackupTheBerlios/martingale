/* WARANTY NOTICE AND COPYRIGHTThis program is free software; you can redistribute it and/ormodify it under the terms of the GNU General Public Licenseas published by the Free Software Foundation; either version 2of the License, or (at your option) any later version.This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty ofMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See theGNU General Public License for more details.You should have received a copy of the GNU General Public Licensealong with this program; if not, write to the Free SoftwareFoundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.Copyright (C) Michael J. Meyermatmjm@mindspring.comspyqqqdia@yahoo.com*/

#ifndef martingale_random_h    
#define martingale_random_h

#include "TypedefsMacros.h"


// Constants used by the Mersenne Twister
/* Period parameters */  
#define MT_N 624
#define MT_M 397
#define MT_D 221               // MT_N-MT_M
#define MATRIX_A 0x9908b0df    // constant vector a 
#define UPPER_MASK 0x80000000  // most significant w-r bits 
#define LOWER_MASK 0x7fffffff  // least significant r bits 

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


/*! \file Random.h
 *  Random number generators: uniform and standard normal.
 *  Uniform reals are based on the Mersenne Twister taken from
 *  <a href="http://www.math.keio.ac.jp/~nisimura/random/real2/mt19932-2.c">
 *  A C-program for MT19937</a>. See also 
 *  <a href="http://www.math.keio.ac.jp/matumoto/emt.html">Matsumoto</a>.
 */


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Random)


/** Mersenne Twister MT19937. Minimal rewrite to C++ of a
 *  <a href="http://www.math.keio.ac.jp/~nisimura/random/real2/mt19932-2.c">
 *  A C-program for MT19937</a> by Makoto
 *  <a href="http://www.math.keio.ac.jp/matumoto/emt.html">Matsumoto</a>
 *  and  Takuji Nishimura. 
 *  License and copyright statement at the bottom of the file 
 *  MersenneTwister.h.
 */
class MersenneTwister {
	
    typedef unsigned long ulong;
    ulong* mt;                     // the array for the state vector 
    int mti;                       // mti==N+1 means mt[N] is not initialized 
	
	// mag01[x] = x * MATRIX_A  for x=0,1 
	static unsigned long mag01[2];
		
public : 

/** @param seed seed used to initialize state vector
 */
MersenneTwister(ulong seed=4357);

/** Next uniform real in (0,1). */
Real u01();

}; // end MersenneTwister



/**  Random draw from {1,-1} */
int sign();

/**  Cumulative normal distribution function */
Real N(Real x);
	
/**  Inverse normal cumulative distribution function */
Real nInverse(Real x);
 
/** Uniform (0,1) real based on MT19937. */
Real U01();

/** Uniform random integers in [0,n-1) */
int Uint(int n);
   
/** Standard normal deviate based on MT19937 and nInverse().*/
Real sTN();
   

MTGL_END_NAMESPACE(Random)   
MTGL_END_NAMESPACE(Martingale)


/* A C-program for MT19937: Real number version((0,1)-interval) */
/* (2001/9/28)                                                  */
/*   genrand() generates one pseudorandom real number (double)  */
/* which is uniformly distributed on (0,1)-interval, for each   */
/* call. sgenrand(seed) sets initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be       */
/* called once. (seed is any 32-bit integer.)                   */
/* Integer generator is obtained by modifying two lines.        */
/*   Coded by Takuji Nishimura, considering the suggestions by  */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.            */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997, 1999, 2001                                  */
/*    Makoto Matsumoto and  Takuji Nishimura.                      */
/*                                                                 */
/* Any feedback is very welcome. For any question, comments,       */
/* see http://www.math.keio.ac.jp/matumoto/emt.html or email       */
/* matumoto@math.keio.ac.jp                                        */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */


#endif
 
