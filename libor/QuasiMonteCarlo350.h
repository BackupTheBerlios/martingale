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


/*
 * QuasiMonteCarlo.h
 *
 * Created on March 21, 2003, 12:45 PM
 */
 
 
#ifndef martingale_lowdiscrepancysequence_h    
#define martingale_lowdiscrepancysequence_h    

#include "TypedefsMacros.h"
#include "Array.h"

MTGL_BEGIN_NAMESPACE(Martingale)



/*! \file QuasiMonteCarlo.h
 *  This file declares the basic interface to low discrepancy sequences and
 *  two implementations: Sobol and Niederreiter-Xing sequences. Dimension is 
 *  limited to 300 in the case of the Sobol generator and 20 in the case of the
 *  Niederreiter-Xing generator. The Sobol generator can very easily be extended
 *  to unlimited dimensions.
 */

/*******************************************************************************
 *
 *                     LOW DISCREPANCY SEQUENCE
 *
 ******************************************************************************/

/** <p>Interface and methods to compute points, quasinormal vectors and the 
 *  L2-discrepancy for low discrepancy sequences.</p>
 *
 * @author  Michael J. Meyer
 */
class LowDiscrepancySequence {
	
protected:
    
    int dim;           // dimension
    int index;         // index of current point in sequence
    
    RealArray1D x;     // current uniform point
    RealArray1D z;     // current quasi normal transform of x
	
public:
	
// ACCESSORS
	
	/** the dimension of the sequence 
	 */
	int getDimension(){ return dim; }
	
	/** empty default implementation */
	virtual void restart(){ };
	
	
// CONSTRUCTOR
    
    /** Constructor
     * @param d dimension
     */
    LowDiscrepancySequence(int d);
	
	/** destructor */
	virtual ~LowDiscrepancySequence(){ }
    
    

// POINT GENERATION

    /** Returns the next point in the sequence.
     */
    virtual const RealArray1D& nextPoint() = 0;
    
    /** <p>Writes the next point of the sequence into 
     *  the array r. Useful when discrepancy is computed.</p>
     */
    void writeNextPoint(RealArray1D& r);
     
    
    
// L^2-DISCREPANCY 
    
    
    /** <p>The \f$L^2\f$-discrepancy of the first N points
     *  r[j][.], j=0,...,N-1.</p>
     *
     * @param N number of points.
     * @param r array of the first N points r[j][.], j=0,...,N-1.
     */    
     Real l2Discrepancy(int N, const RealArray2D& r);
	 
        
    /** <p>The \f$L^2\f$-discrepancy of the first N points. This
     *  computes the \f$L^2\f$-discrepancy \f$T_{N+1}\f$ of the first N+1 points 
     *  r[j][.], j=0,...,N from the discrepancy \f$x=T_N\f$ of the first N points
     *  based on the recursion followed by the quantities \f$h_N=(NT_N)^2/2\f$.</p> 
     *
     * <p>Use this method if the \f$L^2\f$-discrepancy is to be computed for all
     * n=1,2,3,... (number of points).</p>
     *
     * @param N number of points.
     * @param r array of the first N points r[j][], j=0,...,N of the sequence.
     * @param T_N discrepancy \f$T_N\f$.
     */
   
     Real l2Discrepancy(int N, const RealArray2D& r, Real T_N);     
     

// TRANSFORM UNIFORM -> MULTINORMAL

	 
     /** The transform of the next uniform point in the sequence to a
      *  quasinormal vector. Method of transform: coordinatewise inverse
      *  normal CDF.</p>
      *
      */
     const RealArray1D& nextQuasiNormalVector();
	 
	 
private:
	 
	 /** Helper function for \f$L^2\f$-discrepancy */
     // r[j][] is the jth point in the low discrepancy sequence
     Real product(int i, int j, const RealArray2D& r);
         
     
	 /** Helper function for \f$L^2\f$-discrepancy */
     // r[j][] is the jth point in the low discrepancy sequence
     Real product(int i, const RealArray2D& r);
     
     
	 /** Helper function for \f$L^2\f$-discrepancy */
     // r[j][] is the jth point in the low discrepancy sequence
     Real productSQ(int i, const RealArray2D& r);


}; // end LowDiscrepancySequence



/*******************************************************************************
 *
 *                     NIEDERREITER-XING SEQUENCE
 *
 ******************************************************************************/


/** <p>Niederreiter-Xing low discrepancy sequence with basis b=2 in dimension 
 *  at most 20. Uses the Gray code counter and bitwise operations for very fast
 *  point generation.</p>
 *  
 *  <p>Use \f$N=2^n-1, n=1,2,\dots\f$  points for QMC integration.
 *  At \f$N=2^n-1\f$ the Gray code counter G(k) is in sync with the integer sequence
 *  \f$k=1,2,\dots\f$ again, that is \f$\{G(1),\dots,G(N)\}=\{1,\dots,N\}.</p>
 *
 *  <p>Contains the 
 *  <a href="http://www.dismat.oeaw.ac.at/pirs/niedxing.html"> 
 *  generator matrices</a> for dimensions j=4...20 in the array 
 *  gMC decoded from the row encoded form of the download and reencoded
 *  column by column for use with the Gray code counter method in the
 *  NX point generation.</p>
 *  
 *  <p>The columns of the binary matrices are encoded as a decimal integers 
 *  with the bits in each column forming the digits of the encoding modulo two. 
 *  The significance of the bits increases bottom up (the top bit is the most
 *  significant). This ordering of the bits is necessary for the fast point
 *  generation with the Gray code counter.</p>
 *
 * @author Michael J Meyer
 */
class NX : public LowDiscrepancySequence {
    
    static const int M=30;               // we are using 30 bit integers
    static const Real N=1073741842;      // 2^30
    
    IntArray1D x_int;    // current vector of NX integers
    
    /** Array of generator matrices. These start at dimension
	 *  d=4, thus shift d->d-4 in generator matrix index. In each dimension
	 *  there is a sequence of binary generator matrices which are encoded 
	 *  column by column as decimal integers.
	 *
	 *  For any dimension dim set d=max(dim,4). Then genMats[d-4] is the encoding 
	 *  for the generator matrices C(j) for dimension dim and genMats[d-4][j] is 
	 *  the column encoding of the generator matrix C(j) for
     *  coordinate j in the current dimension: for 0<=k<M
     *  genMats[d-4][j][k] is the kth column col_k of C(j) encoded as a decimal 
     *  integer with binary digits the entries of col_k increasing in
     *  significance from the bottom up.
     */
    int*** genMats; 
	
	int d;  // max(dim,4), see above.
	
	
public:

    /** <p>NX low discrepancy sequence in dimension <code>dim</code>.</p>
     *
     * @param dim dimension of the sequence.
     */
    NX(int dim);
	
	~NX();
    void restart();
        
    /** The next nx point in the unit cube [0,1]^dim.
     */
    const RealArray1D& nextPoint();
    
}; // end NXT




/*******************************************************************************
 *
 *                     SOBOL SEQUENCE
 *
 ******************************************************************************/


/** <p>Generator for the Sobol sequence. Uses the Gray code counter and bitwise 
  operations for very fast point generation.</p>
  
  <p>Use \f$N=2^n-1, n=1,2,\dots\f$  points for QMC integration.
  At \f$N=2^n-1\f$ the Gray code counter G(k) is in sync with the integer sequence
  \f$k=1,2,\dots\f$ again, that is \f$\{G(1),\dots,G(N)\}=\{1,\dots,N\}.</p>

  <p><b>Dimension:</b> the Sobol generator is implemented only in dimension
  at most 300. The current implementation relies on primitive polynomials 
  and initialization numbers from the book [J]:
  <i>Monte Carlo Methods in Finance</i> by Peter Jaeckel, Wiley,
  ISNB 047149741X. The CD sold with the book contains millions of 
  primitive polynomials allowing you to extend the generator to millions of 
  dimensions.</p>
 
  <p>If the dimension is small low discrepancy sequences are significantly 
  better Monte Carlo integrators than uniform sequences while this advantage 
  seems to fade as the dimension increases at least if the number N of points
  is restricted to values that are realistic in applications.</p>

  <p>This would argue that we apply the Sobol sequence to a small number of
  important dimensions while driving the remaining dimensions with a uniform
  sequence. On the other hand [J] presents evidence that the Sobol sequence
  keeps up with the uniform sequence at any number N of points even in high 
  dimensions if the initialization numbers are chosen properly.</p>

  <p>In this regard it should be noted that even the best uniform random number 
  generator, the Mersenne Twister is only known to deliver an equidistributed
  sequence up to dimension 623. If the sequence is not equidistributed
  we do not know wether the Monte Carlo integral converges to the true integral
  as the number N of points inreases to infinity. Low discrepancy sequences on 
  the other hand are equidistributed in every dimension and so the Monte
  Carlo integral is guarenteed to converge to the true value of the integral.
  </p>

  <p>The reader is advised to consult [J] for a detailed description of 
  techniques to reduce effective dimensionality and much additional 
  source code related to Monte Carlo simulation. It is an excellent reference 
  on the topic.</p>
 */
class Sobol : public LowDiscrepancySequence {
    
    static const int bits=32;             // we are using 32 bit integers
    static const Real N=4294967296.0;     // 2^32 to big for int
    
    
    UnsignedLongArray2D v;        // v[k] - array of direction numbers for dimension k
    IntArray2D p;                 // p[k] - coefficient array of the k-th primitive polynomial
    IntArray1D g;                 // g[k] - degree of the k-th primitive polynomial
    UnsignedLongArray1D x_int;    // current vector of Sobol integers
	
	int index;      // index of current Sobol point
    
public:
    
// THE SOBOL POINTS

    void restart();
    

   /** The next Sobol point in the unit cube [0,1]^dim.
    */
    const RealArray1D& nextPoint();
    
    
	
// CONSTRUCTOR
	
	
    /** <p>A primitive polynomial p(x) modulo 2 is encoded by a pair of numbers 
     *  (d,n) as follows: d=degree(p), the leading and trailing coefficient of p
     *  are 1 and the intermediate coefficients are the bits of n in the binary
     *  representation of n: for example the polynomial</p>
     *
     * <center> 1+x+x^2+x^4+x^5 </center>
     *
     * <p> with coefficients (1)1101(1) is encoded as (5,n) with n=1101=13.
     *  In other words the least significant bit of n corresponds to the 
     *  second highest power of x etc.</p>
     *
     *  <p> The routine allocates the coefficient array p[k] of this 
     *  polynomial and writes the coefficients into the array with powers of x 
     *  decreasing left to right.</p>
     *
     * @param k polynomial is to be stored as p[k]
     * @param n,d encodings of polynomial
     */
     void read_prim_pol(int k, int n, int d);
	 
	 /** Test function, prints the v[j][k] initialization.
	  */
     void printInitialization();


     /** Constructor
      * @param dim dimension of the Sobol sequence.
      */
     Sobol(int dim);   
    
}; // end Sobol


 


MTGL_END_NAMESPACE(Martingale)

#endif

