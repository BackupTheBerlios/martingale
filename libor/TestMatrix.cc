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



#include "TypedefsMacros.h"
#include "TestMatrix.h"
#include "Matrix.h"
#include "Utils.h"
#include "Random.h"
#include <complex>


MTGL_BEGIN_NAMESPACE(Martingale)


void 
testMatrixInverse()
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





void 
testMatrixMultiply() {
	
Real sq1[5][5]=
{{2,-1,3,-2,4},
 {1,0,-2,-1,3},
 {-2,-1,0,-3,2},
 {3,0,-1,1,2},
 {2,3,-4,-2,1}};
RealMatrix SQ1(sq1), SQ11(sq1);					 

 
Real sq2[5][5]=
{{1,-1,2,-3,4},
 {-1,2,-2,-1,1},
 {-3,-1,0,-1,2},
 {3,2,-1,-2,2},
 {2,-1,3,-2,1}};
RealMatrix SQ2(sq2);							 
 
 
Real utr1[5][5]  = 
{{2,-1,3,-2,4},
 {0,3,-2,-1,3},				         
 {0,0,4,-1,2},					    
 {0,0,0,2,-3},
 {0,0,0,0,4}};
UTRRealMatrix UTR1(utr1), UTR11(utr1);								 

 
Real utr2[5][5]  = 
{{3,1,2,-3,-1},
 {0,1,-2,1,-3},
 {0,0,2,-3,1},
 {0,0,0,1,2},
 {0,0,0,0,-3}};
UTRRealMatrix UTR2(utr2);							 

 
Real ltr1[5][5]  = 
{{2,0,0,0,0},
 {-3,2,0,0,0},
 {3,-2,4,0,0},
 {-2,-1,5,3,0},
 {-3,1,-1,0,4}};
LTRRealMatrix LTR1(ltr1), LTR11(ltr1);								 
 
 
Real ltr2[5][5]  = 
{{3,0,0,0,0},
 {-2,1,0,0,0},
 {2,-1,4,0,0},
 {-1,1,4,-5,0},
 {-2,2,4,-1,-3}};
LTRRealMatrix LTR2(ltr2);					  


Real precision=0.001,       // maximum acceptable relative error in percent
     epsilon=0.00000000001; // zero denominator reset to epsilon
 

// tests of rectangular matrix products
cout << "\n\n\nTesting rectangular matrix products:\n";
 
// SQ1*SQ2
Real sq1_sq2[5][5]=
{{-4, -15, 20, -12, 13},
 {10, -4, 12, -5, 1},
 {-6, -8, 7, 9, -13},
 {13, -2, 11, -14, 14},
 {7, 3, 3, -3, 0}};
RealMatrix SQ1_SQ2(sq1_sq2);
 
SQ1_SQ2.testEquals(SQ1*=SQ2,precision,epsilon,"SQ1*=SQ2 test");


// SQ1*SQ2' : 
Real sq1_sq2t[5][5]=
{{31, -4, 5, 13, 22},
{12, 7, 4, 13, 1},
{16, 5, 14, 2, 5},
{6, 0, -6, 12, 3},
{1, 15, -5, 22, -6}};
RealMatrix SQ1_SQ2t(sq1_sq2t);

SQ1_SQ2t.testEquals(SQ11^=SQ2,precision,epsilon,"SQ1^=SQ2 test");


// tests of lower triangular matrix products
cout << "\n\n\nTesting lower triangular matrix products:\n";

// LTR1*LTR2 : 
Real ltr1_ltr2[5][5]=
{{6, 0, 0, 0, 0},
{-13, 2, 0, 0, 0},
{21, -6, 16, 0, 0},
{3, -3, 32, -15, 0},
{-21, 10, 12, -4, -12}};
LTRRealMatrix LTR1_LTR2(ltr1_ltr2);

LTR1_LTR2.testEquals(LTR1*=LTR2,precision,epsilon,"LTR1*=LTR2 test");

// LTR1*UTR2' : 
Real ltr1_utr2t[5][5]=
{{6, 0, 0, 0, 0},
{-7, 2, 0, 0, 0},
{15, -10, 8, 0, 0},
{-6, -8, 1, 3, 0},
{-14, -9, 2, 8, -12}};
LTRRealMatrix LTR1_UTR2t(ltr1_utr2t);

LTR1_UTR2t.testEquals(LTR11^=UTR2,precision,epsilon,"LTR1^=UTR2 test");


// tests of upper triangular matrix products
cout << "\n\n\nTesting upper triangular matrix products:\n";

// UTR1*UTR2 : 
Real utr1_utr2[5][5]=
{{6, 1, 12, -18, -12},
{0, 3, -10, 8, -22},
{0, 0, 8, -13, -4},
{0, 0, 0, 2, 13},
{0, 0, 0, 0, -12}};
UTRRealMatrix UTR1_UTR2(utr1_utr2);

UTR1_UTR2.testEquals(UTR1*=UTR2,precision,epsilon,"UTR1*=UTR2 test");



// UTR1*LTR2' : 
Real utr1_ltr2t[5][5]=
{{6, -5, 17, 19, -4},
{0, 3, -11, 0, -10},
{0, 0, 16, 21, 11},
{0, 0, 0, -10, 7},
{0, 0, 0, 0, -12}};
UTRRealMatrix UTR1_LTR2t(utr1_ltr2t);

UTR1_LTR2t.testEquals(UTR11^=LTR2,precision,epsilon,"UTR1^=LTR2 test");

} // end testMatrixMultiply



void 
reportMatrixElements(Real** data, int dim, int matrix_type)
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



void 
testMatrixExponentials()
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
		int type=3; cin>>type;
		
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




void 
testSymmetricMatrixExponentials()
{
    int repeat=1;
	while(repeat==1){
	
	    std::cout << endl << endl
		          << "A symmetric matrix A will be initialized with random entries from [-r,r]." << endl
		          << "We then compute the products exp(A) by diagonalization of A " 
		          << " and by series expansion." << endl
		          << "It is then tested if the two results agree." << endl;
		
		std::cout << "Enter r: ";
		Real r; cin>>r;
		
	    std::cout << "Enter dimension: ";
	    int dim; std::cin>>dim; 

	    RealMatrix A(dim); 
		UTRRealMatrix U(dim);
	    for(int i=0;i<dim;i++)
	    for(int j=i;j<dim;j++) A(j,i)=A(i,j)=U(i,j)=r*Random::U01(); 

		RealMatrix E(A.exp()), E1(U.matrixFunction(&exp));
	    
		Real precision=0.001,       // maximum acceptable relative error in percent
		     epsilon=0.00000000001; // zero denominator reset to epsilon
	
		E.testEquals(E1,precision,epsilon,"The matrix exponentials agree:");
	
	    std::cout << endl << endl << endl
	              << "Repeat (1/0)? repeat= ";
	    cin >> repeat;
    } // end main loop   
		
} // end testMatrixExponentials




void 
testComplexExponential()
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
