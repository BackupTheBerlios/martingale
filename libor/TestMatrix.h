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


#include <complex>
#include "TypedefsMacros.h"
#include "Matrices.h"
#include "Utils.h"
#include "Random.h"


MTGL_BEGIN_NAMESPACE(Martingale)

/** A collection of free standing short test programs.
 */
 



/*******************************************************************************    
    
     TEST OF UPPER TRIANGULAR MATRIX INVERSE MATRIX TO MATRIX MULTIPLY  
	
*******************************************************************************/



/** Allocates upper triangular matrix and tests the inverse.
 */
void testMatrixInverse()
{
	int n=6;
	UTRRealMatrix U(n);
	for(int i=0;i<n;i++)
	for(int j=i;j<n;j++) U(i,j)=1+Random::U01();
	
    UTRRealMatrix I(n); for(int i=0;i<n;i++) I(i,i)=1.0;
		
	UTRRealMatrix& Uinv=U.inverse();
	UTRRealMatrix Uinv1(Uinv);
	
	Uinv*=U;
	U*=Uinv1;
	
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.000000001; // zero denominator reset to epsilon
	
	Uinv.testEquals(I,precision,epsilon,"Uinv*U=I");
	cout << Uinv << endl;
	U.testEquals(I,precision,epsilon,"U*Uinv=I");
	cout << U << endl;
	
} // testMatrixInverse




/** Test of the various matrix products in dimension 3 by 3 against
 *  products computed by hand. Ridiculous, but more massive tests of matrix
 *  multiplication are embedded in the tests of the matrix exponential.
 */

void testMatrixMultiply()
{
    Real x[3][3]={{2,1,3},{0,2,-1},{1,-1,3}};
	RealMatrix X(x), X1=X;
		
	Real u[3][3]={{1,2,-1},{0,-1,3,},{0,0,2}};
	UTRRealMatrix U(u), U1=U;
		
    Real l[3][3]={{1,0,0},{2,-1,0},{-1,3,2}};
	LTRRealMatrix L(l), L1=L;
		
	// X^U hand computed
	Real xwu[3][3]={{1,8,6},{5,-5,-2},{-4,10,6}};
	RealMatrix XwedgeU(xwu);
		
     // X*U hand computed
	Real xsu[3][3]={{2,3,7},{0,-2,4},{1,3,2}};
	RealMatrix XstarU(xsu);
		
    // U^L hand computed
	Real uwl[3][3]={{1,0,3},{0,1,3},{0,0,4}};
	UTRRealMatrix UwedgeL(uwl);
		
    // L^U hand computed = (U^L)'
	Real lwu[3][3]={{1,0,0},{0,1,0},{3,3,4}};
	LTRRealMatrix LwedgeU(lwu);
		
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
	XwedgeU.testEquals(X1^=U,precision,epsilon,"X^U test");
	XstarU.testEquals(X*=U,precision,epsilon,"X*U test");
	UwedgeL.testEquals(U1^=L,precision,epsilon,"U^L test");
	UwedgeL.testEquals((L^=U).transpose(),precision,epsilon,"U^L=(L^U)' test");	

    
} // end testMatrixMultiply


/*******************************************************************************    
    
              TEST OF MATRIX EXPONENTIALS  
	
*******************************************************************************/


/** <p>Computes the exponentials <code>H=exp(A), K=exp(-A)</code> for a randomly 
 *  initialized matrix <code>A</code> then prints the products <code>HK,KH</code>
 *  These must be the identity matrix or close to it. Goes through all the matrix
 *  classes: lower and upper triangular and square.</p>
 *
 * <p>The exponential <code>exp(A)</code> is computed as 
 * <code>lim_m(1+A/m)^m, m=2^n, n->oo</code>. The user chooses the precision 
 * level <code>n</code>.
 */

void reportMatrixElements(Real** data, int dim, int matrix_type)
{
    Real min=100000000000000.0, 
	     max=-100000000000000.0;

    // diagonal smallest and largest off diagonal element
	switch(matrix_type){
		
		case 0:
	    for(int i=0;i<dim;i++)
	    for(int j=i+1;j<dim;j++){
		
            Real aij=data[i][j-i];
			if(aij>max)max=aij;
		    if(aij<min)min=aij;
	    }
		std::cout << endl << "Diagonal:" << endl;
    	for(int i=0;i<dim-1;i++) std::cout << data[i][0] << ", ";
	    std::cout << data[dim-1][0] << endl;
        break;
		
		case 1:
	    for(int i=0;i<dim;i++)
	    for(int j=0;j<i;j++){
		
            Real aij=data[i][j];
			if(aij>max)max=aij;
		    if(aij<min)min=aij;
	    }
		std::cout << endl << "Diagonal:" << endl;
    	for(int i=0;i<dim-1;i++) std::cout << data[i][i] << ", ";
    	std::cout << data[dim-1][dim-1] << endl;
		
        break;
		
	    case 2:
	    for(int i=0;i<dim;i++)
	    for(int j=0;j<dim;j++){
		
            Real aij=data[i][j];
			if((i!=j)&&(aij>max)) max=aij;
		    if((i!=j)&&(aij<min)) min=aij;
	    }
		std::cout << endl << "Diagonal:" << endl;
	    for(int i=0;i<dim-1;i++) std::cout << data[i][i] << ", ";
	    std::cout << data[dim-1][dim-1] << endl;
		
	} // end switch
    
	std::cout << "Smallest off diagonal element: " << min << endl;	
	std::cout << "Largest off diagonal element: " << max << endl;

} // end report




void testMatrixExponentials()
{
    int repeat=1;
	while(repeat==1){
	
	    std::cout << endl << endl
		          << "A matrix A will be initialized with random entries from [-r,r]." << endl
		          << "We then compute the products P=exp(A)exp(-A) and Q=exp(-A)exp(A)." << endl 
		          << endl
		          << "These should be equal to the identity matrix, to check this" << endl 
		          << "we print the norm, diagonal and largest and smallest off diagonal element " << endl 
		          << "of P and Q." << endl
		          << endl
		          << "The parameter a controls the matrix norm." << endl
		          << "If this norm is less than one the computation is extremely accurate." << endl
		          << "As the norm increases accuracy is lost." << endl
		          << "Eventually the computation blows up spectacularly even under long double precision." << endl
		          << "Square matrices blow up far more quickly than triangular ones." << endl
		          << "Note however that the errors have to be assessed relative to the norms of P and Q."
		          << endl;
		
		std::cout << "Enter r: ";
		Real r; cin>>r;
		
	    std::cout << "Which type of matrix? (required)" << endl
		          << "Upper triangular (0)" << endl
		          << "Lower triangular (1)" << endl
		          << "Square (2)" << endl
		          << "matrix type = ";
		int type=2; cin>>type;
		
	    std::cout << "Enter dimension: ";
	    int dim; std::cin>>dim; 
	
		if(type==0){

	        UTRRealMatrix A(dim), B(dim);
	        for(int i=0;i<dim;i++)
	        for(int j=i;j<dim;j++){
		
		        A(i,j)=r*Random::U01(); 
		        B(i,j)=-A(i,j);
	        }
		    UTRRealMatrix E(A.exp()), E1(E);
	        UTRRealMatrix F(B.exp());
		    std::cout << endl << endl << "Norm of A : " << A.norm()
			                  << endl << "Norm of exp(A): " << E.norm();
		    E*=F; std::cout << endl << "exp(A)exp(-A): ";  reportMatrixElements(E.getData(),dim,0);
	        F*=E1; std::cout << endl << "exp(-A)exp(A): "; reportMatrixElements(F.getData(),dim,0);

		} 
			
		if(type==1){

	        LTRRealMatrix A(dim), B(dim);
	        for(int i=0;i<dim;i++)
	        for(int j=0;j<=i;j++){
		
		        A(i,j)=r*Random::U01(); 
		        B(i,j)=-A(i,j);
	        }
		    LTRRealMatrix E(A.exp()), E1(E);
	        LTRRealMatrix F(B.exp());
		    std::cout << endl << endl << "Norm of A : " << A.norm() 
			                  << endl << "Norm of exp(A): " << E.norm();        
		    E*=F; std::cout << endl << "exp(A)exp(-A):";  reportMatrixElements(E.getData(),dim,1);
	        F*=E1; std::cout << endl << "exp(-A)exp(A):"; reportMatrixElements(F.getData(),dim,1);
		}
			
		if(type==2){

	        RealMatrix A(dim), B(dim);
	        for(int i=0;i<dim;i++)
	        for(int j=0;j<dim;j++){
		
		        A(i,j)=r*Random::U01(); 
		        B(i,j)=-A(i,j);
	        }
		    RealMatrix E(A.exp()), E1(E);
	        RealMatrix F(B.exp());
		    std::cout << endl << endl << "Norm of A : " << A.norm()
			                  << endl << "Norm of exp(A): " << E.norm();
		    E*=F; std::cout << endl << "exp(A)exp(-A):";   reportMatrixElements(E.getData(),dim,2);
	        F*=E1; std::cout << endl << "exp(-A)exp(A): "; reportMatrixElements(F.getData(),dim,2);
		}
	
	    std::cout << endl << endl << endl
	              << "Repeat (1/0)? repeat= ";
	    cin >> repeat;
    } // end main loop   
		
} // end testMatrixExponentials





/** We compute the exponential exp(A) of the 2 by 2 matrix A={{a,b},{-b,a}}. 
 *  These matrices are isomorphic to the complex numbers via A <-> a+ib. 
 * We then compare our matrix exponential to the STL complex exponential.
 */

void testComplexExponential()
{
     int repeat=1;
	 while(repeat==1){
		 
		 std::cout << "Computing the exponential exp(A) where A={{a,b},{-b,a}}." << endl
		           << "This matrix is the complex number z=a+ib." << endl
		           << "The matrix exponential will be compared to the STL exp(z)." << endl
		           << endl;
		 Real a,b;
		 std::cout << "Enter a: "; cin>>a;
		 std::cout << "Enter b: "; cin>>b;	 
		 Real Z[2][2]={{a,b},{-b,a}};
		 RealMatrix A(Z);
		 RealMatrix E(A.exp());
		 
		 std::complex<Real> z(a,b), e=exp(z);

		 std::cout << "\nMatrix A: " << A 
		           << "Corresponding complex number z: " << z 
		           << "\n\nMatrix exponential exp(A):" << E 
		           << "Complex exponential exp(z): " << e;
		 
		 std::cout << "\n\nRepeat? (0/1), repeat = ";
		 cin>>repeat;
	 } // end while
} // end testComplexExponential


	
	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 