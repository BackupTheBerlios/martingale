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

#ifndef martingale_testmatrix_h    
#define martingale_testmatrix_h


#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)


/** A collection of free standing short test programs.
 */
 


/*******************************************************************************    
    
     TEST OF UPPER TRIANGULAR MATRIX INVERSE MATRIX TO MATRIX MULTIPLY  
	
*******************************************************************************/



/** Allocates upper triangular matrix and tests the inverse.
 */
void testMatrixInverse();



/** Test of the various matrix products in dimension 3 by 3 against
 *  products computed by hand. Ridiculous, but more massive tests of matrix
 *  multiplication are embedded in the tests of the matrix exponential.
 */
void testMatrixMultiply();

/*******************************************************************************    
    
              TEST OF MATRIX EXPONENTIALS  
	
*******************************************************************************/


/** <p>Computes the exponentials <code>H=exp(A), K=exp(-A)</code> for a randomly 
 *  initialized matrix <code>A</code> then prints the products <code>HK,KH</code>
 *  These must be the identity matrix or close to it. The user chooses the matrix
 *  class: lower and upper triangular and square.</p>
 */
void reportMatrixElements(Real** data, int dim, int matrix_type);


/** Randomly intializes a square matrix A and computes exp(A), exp(-A).
 *  Then tests how close exp(A)*exp(-A) and exp(-A)*exp(A) are to the identity
 *  matrix.
 */
void testMatrixExponentials();


/** Randomly intializes a symmetric matrix A and computes exp(A) in two different ways:
 *  <p>
 *  (A) Using the general approach to matrix exponentials.<br>
 *  (B) Diagonalizing the matrix A and applying exp() to each eigenvalue while leaving
 *  the eigenvectors unchanged.
 *  <p>
 *  It is then tested if the two computations agree.
 */
void testSymmetricMatrixExponentials();





/** We compute the exponential exp(A) of the 2 by 2 matrix A={{a,b},{-b,a}}. 
 *  These matrices are isomorphic to the complex numbers via A <-> a+ib. 
 * We then compare our matrix exponential to the STL complex exponential.
 */
void testComplexExponential();

	
	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 