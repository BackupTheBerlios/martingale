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

#ifndef martingale_matrices_h
#define martingale_matrices_h
#define SUBSCRIPT_CHECK

#include <string>
#include <sstream>
#include <iostream>
#include "TypedefsMacros.h"
#include "tnt/tnt_array1d.h"
#include "tnt/tnt_array2d.h"
#include "jama/jama_eig.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "Array.h"
//#include "Utils.h"                 

/*
 * DataStructures.h
 *
 * Created on February 19, 2003, 6:00 PM
 */

MTGL_BEGIN_NAMESPACE(Martingale)

template<class S> class LTRMatrix;
template<class S> class UTRMatrix;
template<class S> class Matrix;



/***************************************************************************************

                      Useful global functions

***************************************************************************************/



/** The relative error in percent of an approximation to an exact value,
 *  positive if the approx > exact negative otherwise.
 *  Type S must support -(S&) and comparison ">".
 *
 * @param exact true value.
 * @approx approximation.
 * @param epsilon zero rounded up to epsilon to avoid division by zero.
 */
template<typename S>
S relativeError(S exact, S approx, S epsilon)
{
	S d=exact;
	if((-epsilon<d)&&(d<epsilon)) d=epsilon;
 
	return 100*(approx-exact)/d;
} 



// DIAGONALIZATION, RANK REDUCED FACTORIZATION, EIGEN ANALYSIS

/** JAMA object containing eigenvalues (sorted ascending) and ON matrix of associated 
 *  eigenvectors. Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 */
template<typename S>
JAMA::Eigenvalue<S>* eigenDecomposition(const UTRMatrix<S>& C)
{
	int dim=C.getDimension(),
	    b=C.getIndexBase();

	TNT::Array2D<S> cv(dim,dim);
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) cv[i][j]=cv[j][i]=C(i+b,j+b);
	return new JAMA::Eigenvalue<S>(cv);
}

/** <p>Let D be the symmetric matrix with upper half C. This function computes 
 *  the matrix R which best approximates D as a product RR' and has rank r. The matrix R
 *  is computed by diagonalizing D, setting all but the r largest eigenvalues
 *  equal to zero and taking the square root of the remaining eigenvalues.
 *  Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 *
 * <p>Maintains the row index base of C, columns indexed from zero.
 */
template<typename S>
Matrix<S>& rank_Reduced_Root(const UTRMatrix<S>& C, int r)
{
	int dim=C.getDimension(),
	    base=C.getIndexBase();
    JAMA::Eigenvalue<S>* eigval=eigenDecomposition(C);
	// eigenvalues sorted ascending
	TNT::Array1D<S> l(dim);
	eigval->getRealEigenvalues(l);
	// corresponding eigenvectors
	TNT::Array2D<S> U(dim,dim);
	eigval->getV(U);
	
	Matrix<S>& root=*(new Matrix<S>(dim,r,base,0));
	for(int i=0;i<dim;i++)
	for(int j=dim-1;j>=dim-r;j--) root(i+base,dim-j-1)=sqrt(l[j])*U[i][j];
			
	return root;
} // end rankReducedRoot
		
	
/** C is interpreted as the upper half of a multinormal covariance matrix.
 *  This method prints how much variability is captured by the 5 largest
 *  eigenvalues of C. Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 */
template<typename S>
static void factorAnalysis(const UTRMatrix<S>& C, int r)
{
	int dim=C.getDimension();
    JAMA::Eigenvalue<S>* eigval=eigenDecomposition(C);
	// eigenvalues sorted ascending
	TNT::Array1D<S> l(dim);
	eigval->getRealEigenvalues(l);
		
	S traceNorm=l[0];
	for(int i=1;i<dim;i++) traceNorm+=l[i];
			
	S sum=0.0; 
	cout << endl << endl
		 << "Variability captured by the i largest eigenvalues: ";
	for(int i=dim-1;i>=dim-r;i--){
			
		sum+=l[i];
		cout << "\ni=" << dim-i << ": " << 100.0*sum/traceNorm << "%";
	}
} // end factorAnalysis
	
	
/** Let D be the symmetric matrix with upper half C. D must be positive semidefinite.
 *  Computes the best approximate factorization \f$C\simeq RR'\f$, where R has rank r
 *  and returns the relative error in the trace norm
 *  \f[||A||=Tr(A'A)=\sum\nolimits_{ij}A^2_{ij}.\f]
 *  Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 *
 * @param message short description of test.
 */
template<typename S>
void factorizationTest(const UTRMatrix<S>& C, int r, string message="")
{		
	Matrix<S> D(C);	
	Real D_norm=D.norm();
	
	Matrix<S>& R=rank_Reduced_Root(C,r);
	
	R^=R;      // R*=R'
	R*=(-1.0);
	D+=R;      // D=D-RR'
	Real err=D.norm();
			
	cout << message << "\nerror:  "<< 100*err/D_norm << "%";
		
} // end factorizationTest
	



/***************************************************************************************

                      C-Style vectors

***************************************************************************************/


/** <p>Straightforward array based vector with elements of type <code>S</code>.
 *  If <code>S</code> is a user defined type the default constructor will 
 *  be called to determine the size of the vector elements. 
 *  </p>
 *
 * <p>Provides the overloaded operators <code>+=,/=</code> needed for use in
 * <code>RandomObjects</code> as componentwise operations on the coordinate type.
 * </p>
 */
template<class S>
class vector {

protected:

    int dim;         // dimension
	int b;           // index base: indices i=b,b+1,...,b+dim-1
    S* dptr;         // pointer to data

public:
	
// ACCESSORS
	
   /** Dimension. */
   int getDimension() const { return dim; }
   /** Sets new dimension, only downsizing, no new memory allocated. */
   void setDimension(int d) { dim=d; }
   /** Index base b: indices i=b,b+1,...,b+dim-1. */
   int getIndexBase() const { return b; }
   /** Set index base b: indices i=b,b+1,...,b+dim-1. */
   void setIndexBase(int base) { b=base; }
   S* getData() const { return dptr; }
      
   /** Subscripting, no bounds checking */
   S& operator[](int i)
   {
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::checkSubscript(i,b,dim,"Vector");
	   #endif	   
	   return dptr[i-b]; 
    }
     
   const S& operator[](int i) const { return dptr[i-b]; }
   
   
   // CONSTRUCTION - DESTRUCTION
   ~vector(){ delete[] dptr; }
   
   /** Constructor initializes with zeroes
    */
   explicit vector(int d, int base=0) : 
   dim(d), b(base), dptr(new S[d]) 
   { 
	   for(int i=0;i<dim;i++) dptr[i]=0; 
   }

   
   // shallow copy not good enough if we return from function
   vector(const vector& v) : 
   dim(v.getDimension()), b(v.getIndexBase()), dptr(new S[dim])
   {
	   S* vdptr=v.getData(); 
	   for(int i=0;i<dim;i++) dptr[i]=vdptr[i];
   }
   
   
  /** Initialize all  entries with the same value,
   *  index base is b=0.
   *
   * @param d dimension
   * @param v constant entry
   */
   vector(int n, S& v) : 
   dim(n), b(0), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=v;
   }

   
  /** Construct from data array, index base is b=0.
   *
   * @param n dimension
   * @param a data array
   */
   vector(int n, S* a) : dim(n), b(0), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=a[i];
   }
   
   
  /** Construct from data array, index base is b=0.
   * @param a data array
   */
   template<int n>
   vector(S a[n]) : dim(n), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=a[i];
   }

   
   /** Dimension of x must be less than or equal to the dimension
    *  of <code>this</code>, no new memory allocated.
    */
   vector& operator=(const vector& x)
   {
       dim=x.getDimension(); 
	   b=x.getIndexBase();
	   S* xdptr=x.getData();
	   for(int i=0;i<dim;i++) dptr[i]=xdptr[i];
	   return *this;
   } 
   
   
   
// ALGEBRAIC OPERATIONS NEEDED BY RANDOMOBJECT
   
   /** Switch components i,j. Indexing based on the index
    *  base of the vector.
    */
   vector& switchComponents(int i, int j)
   { 
       S tmp=dptr[i-b]; dptr[i-b]=dptr[j-b]; dptr[j-b]=tmp; 
   }
   
   
   /** dimension must be same (not checked), index bases need not be same,
    *  index base remains unchanged.
    */
   vector& operator +=(const vector& x)
   { 
	   S* xdptr=x.getData();
	   for(int i=0;i<dim;i++) dptr[i]+=xdptr[i];
	   return *this;
   }
   
   /** dimension must be same (not checked), index bases need not be same,
    *  index base remains unchanged.
    */
   vector& operator -=(const vector& x)
   { 
	   S* xdptr=x.getData();
	   for(int i=0;i<dim;i++) dptr[i]-=xdptr[i];
	   return *this;
   }
   
   /** division by integer (for use with class RandomObject)
    */
   vector& operator /=(int N)
   { 
	   Real f=1.0/N;
	   for(int i=0;i<dim;i++) dptr[i]*=f;
	   return *this;
   }
   
   
   /** multiplication by scalar
    */
   vector& operator *=(const S& lambda)
   { 
	   for(int i=0;i<dim;i++) dptr[i]*=lambda;
	   return *this;
   }
   

   /** left multiplication by lower triangular matrix
    */
   vector& operator *=(const LTRMatrix<S>& A)
   { 
	   S tmp[dim];
	   S** a=A.getData();
	   for(int i=0;i<dim;i++){
		   
		   S sum=0;
		   for(int j=0;j<=i;j++) sum+=a[i][j]*dptr[j];
		   tmp[i]=sum;
	   }
	   // copy temp into this
	   for(int i=0;i<dim;i++) dptr[i]=tmp[i];
	   return *this;
   }
   
 
   /** left multiplication by upper triangular matrix
    */
   vector& operator *=(const UTRMatrix<S>& A)
   { 
	   S tmp[dim];
	   S** a=A.getData();
	   for(int i=0;i<dim;i++){
		   
		   S sum=0;
		   for(int j=i;j<dim;j++) sum+=a[i][j-i]*dptr[j];
		   tmp[i]=sum;
	   }
	   // copy temp into this
	   for(int i=0;i<dim;i++) dptr[i]=tmp[i];
	   return *this;
   }
   

   /** left multiplication by rectangular matrix
    */
   vector& operator *=(const Matrix<S>& A)
   { 
	   int d=A.getnRows();  // new dimension
	   S* Ax=new S[d];      // new memory
	   S** a=A.getData();
	   for(int i=0;i<d;i++){
		   
		   S sum=0;
		   for(int j=0;j<dim;j++) sum+=a[i][j]*dptr[j];
		   Ax[i]=sum;
	   }

       // free the old memory and reset
	   delete[] dptr; dptr=Ax; dim=d;
	   return *this;
   }
   
   
   
// MEAN, EUCLIDEAN NORM
   
      /** Arithmetic mean of the components. */
   S mean() const
   {  
	   Real sum=0; 
	   for(int i=0;i<dim;i++) sum+=dptr[i];
	   return sum/dim;
   }
   
   /** The euclidean norm. Type S must support an absolute value <code>fabs(S& s)</code>
    *  mapping to Real.
    */
   Real norm() const
   {  
	   Real sum=0, absx_i; 
	   for(int i=0;i<dim;i++){ absx_i=fabs(dptr[i]); sum+=absx_i*absx_i; }
	   return sqrt(sum);
   }
   
   

// TEST FOR EQUALITY


/** Test for entry by entry equality.
 *  Type S must support absolute value and comparison ">".
 *  Equality is defined by an upper bound on the acceptable relative 
 *  error in percent.
 *  
 * @param v vector to compare this with.
 * @param precision upper bound on the acceptable relative error in percent.
 * @param epsilon zero denominator reset to epsilon, should be
 * smaller than relevant orders of magnitude.
 * @param test string printed to identify test.
 */
void testEquals(const vector<S>& v, S precision, S epsilon, string test) const
{
	int vdim=v.getDimension(), vbase=v.getIndexBase(); 
	if((vdim!=dim)||(vbase!=b)){
		
		std::cout << "\n"+test+": failed, dimensions or index range not equal.";
		return;
	}
		
	for(int i=0;i<dim;i++)
	{ 
		S err=relativeError(dptr[i],v[i+b],epsilon);
		if(err>precision){
			
		     std::cout << "\n"+test+": failed, "
		               << "component v["<<i<<"], relative error " << err;
		     return;
		}
	} // end for i
	
    std::cout << "\n"+test+": passed.";
} // end testEquals

	
	      
}; // end vector




/** print vector
 */
template<class S>
std::ostream& operator << (std::ostream& os, const vector<S>& v)
{
	int d=v.getDimension();
    S* vdptr=v.getData();
	os << endl << "vector of dimension " << d << ":" << endl;
    for(int i=0;i<d-1;i++) os << vdptr[i] << ", ";
    os << vdptr[d-1];
    return os << endl << endl;
} // end operator <<

	


/***************************************************************************************
 *
 *                             MATRIX CLASSES
 *
***************************************************************************************/

/*! \file Matrices.h
 * <p>Vectors, square lower triangular, square upper triangular and rectangular matrices 
 * stored in row major order. Barebones naive implementation with no checking of dimensional 
 * compatibility. All operations implemented as straightforward  
 * C-style loops with no optimizations. No attempts to avoid temporaries.
 * To check subscripting bounds #define SUBSCRIPT_CHECK.</p>
 *
 * <p>Matrix indices can be based on an arbitrary base b. This means that the indices
 * used with the subscripting operators are 
 * \f[i,j=b,b+1,...,b+dim-1.\f]
 * This feature allows the use of the natural indices from the application domain.
 * The subscripting operators take care of the necessary translations to the actual 
 * array indices. The index bases are irrelevant if the underlying data array is
 * accessed directly through the data pointers.</p>
 *
 * <p>Matrix multiplication implemented as an operator <code>*=(const Matrix& B)</code>
 * must decide wether <code>this</code> is multiplied with <code>B</code> on the right
 * or on the left. The choices made are the ones which most directly support the
 * applications to the simulation of the Libor dynamics. Use transposition 
 * (AB)'=B'A' to get the other cases.</p>
 *
 * <p>For matrices allocated in row major order the loops for the product AB' have
 * a better memory access pattern (row dot row) than the loops for the product
 * AB (row dot column). There is very little difference in low dimensions but for 
 * 1000 by 1000 matrices AB' is already three times faster than AB.</p>
 *
 * <p>Thus we provide both the operator *= multiplying <code>this</code> on the right 
 * B (slower) and the operator ^= which multiplies on the right with B'.</p>
 */



/***************************************************************************************
 *
 *                             SQUARE LOWER TRIANGULAR MATRICES
 *
***************************************************************************************/

/** <p>Square lower triangular matrix with indices with arbitrary base (the same in both 
 *  dimensions) and entries of type <code>S</code> stored in row major order. 
 *  If <code>S</code> is a user defined type 
 *  the default constructor will be called to determine the size of the matrix elements. 
 *  For the algebraic operations to work the type <code>S</code> should support 
 *  <code>+=(const S&), *(const S&, const S&), *=(Real f)</code>.</p>
 */
template<class S>
class LTRMatrix {

protected:

    int dim;        // dimension
	int base;       // indices i,j=base,...,base+dim-1
    S** dptr;       // pointer to data

 public:

// ACCESSORS 

/** dimension */
int getDimension() const { return dim; }
/** index base b: indices i,j=b,...,b+dim-1. */
int getIndexBase() const { return base; }
/** Set index base b: indices i,j=b,...,b+dim-1. */
void setIndexBase(int b) { base=b; }
/** pointer to data */
S** getData() const { return dptr; }


/** Subscripting, only \f$j\leq i\f$, no bounds checking.
 */
S& operator()(int i, int j)
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,base,base,dim,dim,"LTRMatrix");
	#endif	  
	return dptr[i-base][j-base]; 
}


/** Subscripting, only \f$j\leq i\f$, no bounds checking.
 */
const S& operator()(int i, int j) const 
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,base,base,dim,dim,"LTRMatrix");
	#endif	
	return dptr[i-base][j-base]; 
}


// EQUALITY CHECK


/** Test for entry by entry equality.
 *  Type S must support absolute value and comparison ">".
 *  Equality is defined by an upper bound on the acceptable relative 
 *  error in percent.
 *  
 * @param A matrix to compare this with.
 * @param precision upper bound on the acceptable relative error in percent.
 * @param epsilon zero denominator reset to epsilon, should be
 * smaller than relevant orders of magnitude.
 * @param test string printed to identify test.
 */
void testEquals(const LTRMatrix<S>& A, S precision, S epsilon, string test) const
{
	int d=A.getDimension(), baseA=A.getIndexBase();
	if((d!=dim)||(baseA!=base)){
		
		std::cout << "\n"+test+": failed, dimensions or index range not equal.";
		return;
	}
	
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) { 
		
		S err=relativeError(dptr[i][j],A(i+base,j+base),epsilon);
		if((err>precision)||(-err>precision)){
			
			std::cout << "\n"+test+": failed, "
			          << "component A["<<i<<"]["<<j<<"], relative error " << err;
		    return;
		}
	} // end for i

    std::cout << "\n"+test+": passed.";
} // end testEquals


	

// CONSTRUCTION - DESTRUCTION

/** Constructor for indices with base ib, all components intialized with zeroes.<br>
 *  @param d dimension of square matrix, default 1.
 *  @param b index base (i,j=b,...,b+d-1), default 0.
 */
LTRMatrix(int d=1, int b=0) : 
dim(d), base(b), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++){
		
		dptr[i]=new S[i+1];
		for(int j=0;j<=i;j++) dptr[i][j]=0;
	}
} // end constructor


/** Copy constructor
 */
LTRMatrix(const LTRMatrix& A) : 
dim(A.getDimension()), base(A.getIndexBase()), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++)dptr[i]=new S[i+1];
		
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) dptr[i][j]=A(i+base,j+base);
}

	
/** destructor */
~LTRMatrix()
{
	for(int i=0;i<dim;i++){ delete[] dptr[i]; }
    delete[] dptr;
}


/** Assignement, dimension of B must be less than or equal to 
 *  dimension of <code>this</code> (not checked), no new memory
 *  allocated.
 */
LTRMatrix& operator=(const LTRMatrix& B)
{
    dim=B.getDimension();
	base=B.getIndexBase();
    for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) dptr[i][j]=B(i+base,j+base);
		
	return *this;
}

 
/** Construct from data array, only lower triangular half is used.
 * @param a data array.
 * @param b index base (default 0).
 */
template<int n>
LTRMatrix(S a[n][n], int b=0) : 
dim(n), base(b), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++)dptr[i]=new S[i+1];
	
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) dptr[i][j]=a[i][j];
}


// NORMS

/** Square root of sum of absolute values of entries squared. 
 *  Type S must support absolute value <code>fabs(S& s)</code> mapping to Real.
 */
Real norm() const
{
	Real sum=0.0, absij;
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) { absij=fabs(dptr[i][j]); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm
	



// ALGEBRAIC OPERATIONS


/** Transpose. Must allocate new memory to preserve row major order.
 */
UTRMatrix<S>& transpose() const
{
	UTRMatrix<S>& at=*(new UTRMatrix<S>(dim,base));
    for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) at(i+base,j+base)=dptr[j][i];
	
	return at;
}
	

/** Matrix addition. Dimensions must match (not checked), 
 *  bases need not.
 */
LTRMatrix<S>& operator += (const LTRMatrix<S>& B)
{
	// direct data access for speed
	S** b=B.getData();
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) dptr[i][j]+=b[i][j];
		
	return *this;
} // end operator +=


/** Multiplication of matrix by a Real scalar
 */
LTRMatrix<S>& operator *= (Real f)
{
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) dptr[i][j]*=f;
		
	return *this;
} // end operator +=


/** RIGHT multiplication of <code>this</code> by B. 
 *  Dimensions must match (not checked), bases need not.
 */
// A product of lower triangular matrices is again lower triangular
LTRMatrix<S>& operator *= (const LTRMatrix<S>& B)
{
	LTRMatrix<S> temp(dim); 
	// direct data access for speed
	S** a=dptr;
	S** b=B.getData();
	S** tmp=temp.getData();
	
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++){
		
	    S sum=a[i][i]*b[i][j];  // k=i below
		for(int k=i+1;k<=j;k++) sum+=a[i][k]*b[k][j];
	    tmp[i][j]=sum;
	}
	
	// copy temp into this
    for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++)	a[i][j]=tmp[i][j];
		
	return *this;
} // end operator *=



/** <p>RIGHT multiplication of <code>this</code> by the transpose B'. 
 *  Because of row major allocation the matrix product AB' is faster
 *  than AB.</p>
 *
 *  <p>Dimensions must match (not checked), bases need not.</p>
 */
// A product of lower triangular matrices is again lower triangular
LTRMatrix<S>& operator ^= (const UTRMatrix<S>& B)
{
	LTRMatrix<S> temp(dim); 
	// direct data access for speed
	S** a=dptr;
	S** b=B.getData();
	S** tmp=temp.getData();
	
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++){
		
	    S sum=a[i][j]*b[j][0];
		for(int k=j+1;k<=i;k++) sum+=a[i][k]*b[j][k-j];
	    tmp[i][j]=sum;
	}
	
	// copy temp into this
    for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++)	a[i][j]=tmp[i][j];
		
	return *this;
} // end operator ^=



/** The upper triangular half of the product <code>c=aa'</code>, 
 *  where <code>a=this</code>.
 */
UTRMatrix<S>& aat() const
{
	S** a=dptr;
	UTRMatrix<S>& c=*(new UTRMatrix<S>(dim,base));
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++){
		
		S sum=0.0; 
		for(int k=0;k<=i;k++) sum+=a[i][k]*a[j][k];
		c(i+base,j+base)=sum;
	}
	return c;
} // end aat



// MATRIX FUNCTIONS: some analytic functions f(A), where A=this.

/** Matrix exponential <code>exp(A)</code>. Note <code>exp(A)</code> 
 *  is again lower triangular.
 */
LTRMatrix<S>& exp() const
{
    int n=2;
	Real fn=norm(), f=4.0, g=1.0;
	while((f<fn)){ f*=2.0; n++; }                   // f=2^n
	g=1.0/f; g/=64.0;                               // g=1/2^{n+6}

	// temporary matrices
	LTRMatrix<S> U(dim);            // U=(A/2^{n+6})=g*A 
	UTRMatrix<S> Ut(dim);           // U', we use ^=U ie. *=U'
	LTRMatrix<S> SUk(dim);          // sum of U^k/k!, k=0,1,..,6

    S** a=dptr;
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++){ 
		
		U(i,j)=Ut(j,i)=SUk(i,j)=g*a[i][j];  // U=A/2^{n+6}
		if(i==j) SUk(i,j)+=1.0;             // I+U
	}
	
	// Suk=I+U+U^2/2!+...+U^8/8! 
	for(int k=2;k<=8;k++){ U^=Ut; g=1.0/k; U*=g; SUk+=U; }
	
	// exp(A)=Suk^m, m=2^{n+6};
    LTRMatrix<S>& F=*(new LTRMatrix<S>(SUk));    
	for(int k=0;k<n+6;k++) F*=F;
	F.setIndexBase(base);		
	return F;

} // end exp


	
}; // end LTRMatrix


/** print lower triangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const LTRMatrix<S>& A)
{
	int dim=A.getDimension(), b=A.getIndexBase();
	os << endl << "Lower triangular matrix, dimension " << dim << ":" 
	   << endl << endl;
	for(int i=b;i<dim+b;i++){
		
		for(int j=b;j<i;j++) os << A(i,j) << ", ";
		os << A(i,i) << endl;
	}
    return os << endl << endl;
} // end operator <<





/***************************************************************************************
 *
 *                             SQUARE UPPER TRIANGULAR MATRICES
 *
***************************************************************************************/



/** <p>Square upper triangular matrix with indices with arbitrary base (the same in 
 *  both dimensions) and entries of type <code>S</code>. 
 *  Stored in row major order. If <code>S</code> is a user defined 
 *  type the default constructor will be called to determine the size of the matrix 
 *  elements. For the algebraic operations to work the type <code>S</code> should 
 *  support <code>+=(const S&), *(const S&, const S&), *=(Real f)</code>.
 *  </p>
 */
template<class S>
class UTRMatrix {

protected:

	int dim;       // dimension
	int base;      // indices i,j=base,...,base+dim-1
	S** dptr;      // pointer to data

public:
	
// ACCESSORS

/** nRows=nColumns */
int getDimension() const { return dim; }
/** indices i,j=base,...,base+dim-1 */
int getIndexBase() const { return base; }
/** Set index base b: indices i,j=b,...,b+dim-1. */
void setIndexBase(int b){ base=b; }
/** pointer to data array */
S** getData() const { return dptr; }


/** Subscripting, only \f$i\leq j\f$, no bounds checking.
 */
S& operator()(int i, int j)
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,base,base,dim,dim,"UTRMatrix");
	#endif	
	return dptr[i-base][j-i]; 
}


/** Subscripting, only \f$i\leq j\f$, no bounds checking.
 */
const S& operator()(int i, int j) const 
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,base,base,dim,dim,"UTRMatrix");
	#endif	
	return dptr[i-base][j-i]; 
}
	
	
// PRINTOUT
	
/** Prints the symmetric square matrix of which <code>this</code> is the upper half. 
 */
void print() const
{
	cout << endl << endl
	     << "Symmetric " << dim << " by " << dim << " matrix: "
	     << endl << endl;
	
	S** A=dptr;
	for(int i=0;i<dim;i++){
		
		for(int j=0;j<i;j++)     cout << A[j][i-j] << ", ";
		for(int j=i;j<dim-1;j++) cout << A[i][j-i] << ", ";
		cout << A[i][dim-i-1] << endl;
	}
} // end print
    
	
	
// EQUALITY CHECK

/** Test for entry by entry equality.
 *  Type S must support absolute value and comparison ">".
 *  Equality is defined by an upper bound on the acceptable relative 
 *  error in percent.
 *  
 * @param A matrix to compare this with.
 * @param precision upper bound on the acceptable relative error in percent.
 * @param epsilon zero denominator reset to epsilon, should be
 * smaller than relevant orders of magnitude.
 * @param test string printed to identify test.
 */
void testEquals(const UTRMatrix<S>& A, S precision, S epsilon, string test) const
{
	int d=A.getDimension(), baseA=A.getIndexBase();
	if((d!=dim)||(baseA!=base)){
		
		std::cout << "\n"+test+": failed, dimensions or index range not equal.";
		return;
	}
	
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) { 
		
		S err=relativeError(dptr[i][j-i],A(i+base,j+base),epsilon);
		if((err>precision)||(-err>precision)){
			
			std::cout << "\n"+test+": failed, "
			          << "component A["<<i<<"]["<<j<<"], relative error " << err;
		    return;
		}
	} // end for i

    std::cout << "\n"+test+": passed.";
} // end testEquals



// CONSTRUCTION - DESTRUCTION


/** Constructor, initializes all components with zeroes.<br>
 *  @param d dimension of square matrix 
 *  @param b index base: b<=i,j<b+dim
 */
UTRMatrix(int d=1, int b=0) : 
dim(d), base(b), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++){
		
	     dptr[i]=new S[dim-i];
	     for(int j=i;j<dim;j++)dptr[i][j-i]=0;
	}
} // end constructor


/** Copy constructor
 */
UTRMatrix(const UTRMatrix<S>& A) : 
dim(A.getDimension()), base(A.getIndexBase()), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++)dptr[i]=new S[dim-i];
		
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) dptr[i][j-i]=A(i+base,j+base);
}

	

/** Construct from data array, only upper triangular half is used.<br>
 * @param a data array
 * @param b index base (default 0).
 */
template<int n>
UTRMatrix(S a[n][n], int b=0) : 
dim(n), base(b), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++)dptr[i]=new S[dim-i];
	
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) dptr[i][j-i]=a[i][j];
}
	

/** destructor */
~UTRMatrix()
{
	for(int i=0;i<dim;i++){ delete[] dptr[i]; }
    delete[] dptr;
}


/** Assignement, dimension of B must be less than or equal to 
 *  the dimension of <code>this</code> (not checked), no new
 *  memory allocated.
 */
UTRMatrix<S>& operator =(const UTRMatrix<S>& B)
{
    dim=B.getDimension();
	base=B.getIndexBase();
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) dptr[i][j-i]=B(i+base,j+base);
	
	return *this;
}


// NORMS

/** Square root of sum of entries squared. Type S must support an absolute
 *  value <code>fabs(S& s)<S></code> mapping to Real.
 */
Real norm() const
{
	Real sum=0, absij;
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) { absij=fabs(dptr[i][j-i]); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm

	



// ALGEBRAIC OPERATIONS



/** Transpose. Must allocate new memory to preserve row major order.
 */
LTRMatrix<S>& transpose() const
{
	LTRMatrix<S>& at=*(new LTRMatrix<S>(dim,base));
    for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++) at(i+base,j+base)=dptr[j][i-j];
	
	return at;
}


/** The matrix inverse. */
UTRMatrix<S>& inverse() const
{
	for(int i=0;i<dim;i++){
	  if(dptr[i][0]==0.0){ 
		 cout << "UTRMatrix::inverse(): matrix is singular, quitting."; 
		 exit(0); 
	  }
	} // end for i
		
	UTRMatrix<S>& inverse = *(new UTRMatrix<S>(dim,base));
			
    S** A=dptr; S** I=inverse.getData();	
	// normalize diagonal
	for(int i=0;i<dim;i++) I[i][0]=1.0/A[i][0];       // A[i][0]=A_ii
	

	for(int j=dim-1;j>=0;j--)
	for(int i=j-1;i>=0;i--) 
	// subtract row_j*A_ij/A_ii from row_i
	for(int k=j;k<dim;k++) I[i][k-i]-=I[j][k-j]*A[i][j-i]/A[i][0];
		
	return inverse;
}
		
		

/** Matrix addition. Matrix dimensions must be equal (not checked),
 *  index bases need not be equal.
 */
UTRMatrix<S>& operator += (const UTRMatrix<S>& B)
{
	// direct data access for speed
	S** b=B.getData();
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) dptr[i][j-i]+=b[i][j-i];
		
	return *this;
} // end operator +=


/** Multiplication of matrix by a Real scalar
 */
UTRMatrix<S>& operator *= (Real f)
{
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) dptr[i][j-i]*=f;
		
	return *this;
} // end operator +=



/** RIGHT multiplication of <code>this</code> by B. 
 *  Matrix dimensions must be equal (not checked), index bases
 *  don't have to be equal.
 */
// A product of upper triangular matrices is again upper triangular
UTRMatrix<S>& operator *= (const UTRMatrix<S>& B)
{
	UTRMatrix<S> temp(dim); 
	// direct data access for speed
	S** a=dptr;
	S** b=B.getData();
	S** tmp=temp.getData();
	
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++){
		
		S sum=a[i][0]*b[i][j-i];  // k=i below
		for(int k=i+1;k<=j;k++) sum+=a[i][k-i]*b[k][j-k];
	    tmp[i][j-i]=sum;
	}
	
	// copy temp into this
    for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++)	a[i][j-i]=tmp[i][j-i];

	return *this;
}



/** RIGHT multiplication of <code>this</code> by the transpose B'. 
 *  Because of row major allocation the product AB' is faster than
 *  AB. Matrix dimensions must be equal (not checked), index bases
 *  don't have to be equal.
 */
// A product of upper triangular matrices is again upper triangular
UTRMatrix<S>& operator ^= (const LTRMatrix<S>& B)
{
	UTRMatrix<S> temp(dim); 
	// direct data access for speed
	S** a=dptr;
	S** b=B.getData();
	S** tmp=temp.getData();
	
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++){
		
		S sum=a[i][0]*b[j][i];  // k=i below
		for(int k=i+1;k<=j;k++) sum+=a[i][k-i]*b[j][k];
	    tmp[i][j-i]=sum;
	}
	
	// copy temp into this
    for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++)	a[i][j-i]=tmp[i][j-i];
		
	return *this;
} // end operator ^=


/** The upper triangular half of the product <code>c=aa'</code>, 
 *  where <code>a=this</code>.
 */
UTRMatrix<S>& aat() const
{
	S** a=dptr;
	UTRMatrix<S>& c=*(new UTRMatrix<S>(dim,base));
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++){
		
		S sum=0.0; 
        // note j>=i
		for(int k=j;k<dim;k++) sum+=a[i][k-i]*a[j][k-j];
		c(i+base,j+base)=sum;
	}
	return c;
} // end aat





// MATRIX FUNCTIONS: some analytic functions f(A), where A=this.

/** Matrix exponential <code>exp(A)</code>. Note <code>exp(A)</code> 
 *  is again lower triangular.
 */
UTRMatrix<S>& exp() const
{
    int n=2;
	Real fn=norm(), f=4.0, g=1.0;
	while((f<fn)){ f*=2.0; n++; }           // f=2^n
	g=1.0/f; g/=64.0;                       // g=1/2^{n+6}

	// temporary matrices
	UTRMatrix<S> U(dim);           // U=(A/2^{n+6})=g*A 
	LTRMatrix<S> Ut(dim);          // U', we use ^=U ie. *=U'
	UTRMatrix<S> SUk(dim);         // sum of U^k/k!, k=0,1,..,6

    S** a=dptr;
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++){ 
		
		U(i,j)=Ut(j,i)=SUk(i,j)=g*a[i][j-i];      // U,U'
		if(i==j) SUk(i,j)+=1.0;                   // I+U
	}
	
	// SUk=I+U+U^2/2!+...+U^8/8! 
	for(int k=2;k<=8;k++){ U^=Ut; g=1.0/k; U*=g; SUk+=U; }
	
	// exp(A)=SUk^m, m=2^{n+6};
    UTRMatrix<S>& F=*(new UTRMatrix<S>(SUk));     
	for(int k=0;k<n+6;k++) F*=F;
	F.setIndexBase(base);	
	return F;

} // end exp



/** Cholesky root L (lower triangular, LL'=A where A is the symmetric matrix with upper half <code>this</code>). 
 *  Terminates with error message if non positive definite detected. 
 *  Type S must support += and *, and sqrt(const S&). 0.0 must convert to type S.
 */
LTRMatrix<S>& ltrRoot() const
{
    LTRMatrix<S>& ltr=*(new LTRMatrix<S>(dim,base));
	S** L=ltr.getData();
  
    // computation of Lij reverse induction starting from i,j=0 ascending
	// get U_ij from relation row_i(U) dot row_j(U) = A_ij, j<=i, at that point all
	// U_rs with r<=i,s<=j (and not both equal) are already computed.
    for(int i=0;i<dim;i++)
    for(int j=0;j<=i;j++)                         // j\leq i
    { 
        // r_i(L).r_j(L)-L_ijL_jj
		S sum=0.0;  
        for(int k=0;k<j;k++) sum += L[i][k]*L[j][k];   
      
		S Aij=dptr[j][i-j], R=Aij-sum;                // R-L_ijL_jj=A_ij-r_i(L).r_j(L)=0
        if(j<i) L[i][j]=R/L[j][j];                    // L[j][j] already computed, compute L[i][j]
        else if(R>0) L[j][j]=sqrt(R);
        else{ cerr<<endl<<endl<<"Dimension "<<dim
                  <<", C["<<j<<"]["<<j<<"]-S="<<R<<endl
                  <<"ltrRoot(): Matrix not positive definite:" << endl
				  << *this << endl << "Terminating.";
               exit(1); 
        } // end else
    } //end for i
	
	return ltr;

} //end ltrRoot




/** Upper triangular root U of A, where A is the symmetric matrix of which
 *  <code>this</code>) is the upper half. Satisfies UU'=A. Differs from the 
 *  Cholesky root as it is upper triangular instead of lower triangular.
 *  Terminates with error message if non positive definite detected. 
 *  Type S must support += and *, and sqrt(const S&). 0.0 must convert to type S.
 */
UTRMatrix<S>& utrRoot() const
{
    UTRMatrix<S>& utr=*(new UTRMatrix<S>(dim,base));
	S** U=utr.getData();
	
	int n=dim-1;
  
	// computation of Lij reverse induction starting from i,j=n descending
	// get U_ij from relation row_i(U) dot row_j(U) = A_ij, j>=i, at that point all
	// U_rs with r>=i,s>=j (and not both equal) are already computed.
    for(int i=n;i>=0;i--)
    for(int j=n;j>=i;j--) {                      // j>=i
    
        // r_i(U).r_j(U)-U_ijU_jj
		S sum=0.0;  
        for(int k=j+1;k<dim;k++) sum += U[i][k-i]*U[j][k-j];   

   		S Aij=dptr[i][j-i], R=Aij-sum;           // R-U_ijU_jj=A_ij-r_i(U).r_j(U)=0    
        if(j>i) U[i][j-i]=R/U[j][0];             // U[j][0]=U[j][j-j]=U_jj              
        else if(R>0) U[j][0]=sqrt(R);
        else{ cerr<<endl<<endl<<"Dimension "<<dim
                  <<", A["<<j<<"]["<<j<<"]-S="<<R<<endl
                  <<"utrRoot(): Matrix not positive definite:" << endl
				  << *this << endl << "Terminating.";
               exit(1); 
        } // end else
    } //end for i
	
	return utr;

} //end utrRoot


/** <p>Let C be the symmetric matrix with upper half <code>this</code>. 
 *  This function computes the matrix R which best approximates D as a product RR' and has rank r. 
 *  The matrix R is computed by diagonalizing D, setting all but the r largest eigenvalues
 *  equal to zero and taking the square root of the remaining eigenvalues.
 *
 * <p>The row index base remains the same, columns of the root are indexed from zero.
 */
Matrix<S>& rankReducedRoot(int r) const { return rank_Reduced_Root(*this,r); }


/** Matrix is interpreted as the upper half of a multinormal covariance matrix.
 *  This method prints how much variability is captured by the r largest
 *  eigenvalues of C. Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 */
void analyseFactors(int r) const { factorAnalysis(*this,r); }


/** Let C be the symmetric matrix with upper half <code>this</code>. 
 *  C must be positive semidefinite. Computes the best approximate factorization 
 *  \f$C\simeq RR'\f$, where R has rank r and returns the relative error in the trace 
 *  norm
 *  \f[||A||=Tr(A'A)=\sum\nolimits_{ij}A^2_{ij}.\f]
 *  Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 *
 * @param message put test in context.
 */
void testFactorization(int r, string message="") const
{ factorizationTest(*this,r,message); }



}; // end UTRMatrix



/** print upper triangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const UTRMatrix<S>& U)
{
	os << endl << "Transposed matrix:";
    return os << U.transpose() << endl << endl;
} // end operator <<







/***************************************************************************************
 *
 *                             RECTANGULAR MATRICES
 *
***************************************************************************************/


/** <p>Recangular matrix with dimensions <code>rows,cols</code> and entries of type <code>S</code>
 *  stored in row major order. Indices can have arbitrary bases a,b. This means that we use indices 
 *  <code>
 *  <center>i=a,a+1,...,a+rows-1</center>
 *  <center>j=b,b+1,...,b+cols-1</center>
 *  </code>
 *  <p></p>
 *  in the subscripting operator.
 *  If <code>S</code> is a user defined type the default constructor 
 *  will be called to determine the size of the matrix elements. For the algebraic
 *  operations to work the type <code>S</code> should support 
 *  <code>+=(const S&), *(const S&, const S&), *=(Real f)</code>.</p>
 */
template<class S>
class Matrix {

protected:

    int rows,       // number of rows
	    cols,       // number of columns
	    a,          // row indices i=a,a+1,...,a+rows-1
	    b;          // column indices j=b,b+1,...,b+cols-1
    S** dptr;       // pointer to data

public:
	
// ACCESSORS

/** Number of rows */
int getnRows() const { return rows; }
/** Number of columns */
int getnCols() const { return cols; }

/** Row index base a, row indices i=a,a+1,...,a+rows-1 */
int getRowIndexBase() const { return a; }
/** Row index base set equal to r. */
void setRowIndexBase(int r){ a=r; }

/** Column index base b, column indices j=b,b+1,...,b+cols-1 */
int getColIndexBase() const { return b; }
/** Column index base set equal to r. */
void setColIndexBase(int r){ b=r; }

/** Row index base a, assumption is that both indices
 *  have the same base.
 */
int getIndexBase() const { return a; }
/** Bases for both row and column indices set equal to r.
 */
void setIndexBase(int r){ a=b=r; }

/** pointer to data */
S** getData() const { return dptr; }


/** Subscripting, no bounds checking
 */
S& operator()(int i, int j)
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,a,b,rows,cols,"Matrix");
	#endif		
	return dptr[i-a][j-b]; 
}


/** Subscripting, no bounds checking
 */
const S& operator()(int i, int j) const 
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,a,b,rows,cols,"Matrix");
	#endif	
	return dptr[i-a][j-b]; 
}


	
// EQUALITY CHECK


/** Test for entry by entry equality.
 *  Type S must support absolute value and comparison ">".
 *  Equality is defined by an upper bound on the acceptable relative 
 *  error in percent.
 *  
 * @param A matrix to compare this with.
 * @param precision upper bound on the acceptable relative error in percent.
 * @param epsilon zero denominator reset to epsilon, should be
 * smaller than relevant orders of magnitude.
 * @param test string printed to identify test.
 */
void testEquals(const Matrix<S>& A, S precision, S epsilon, string test) const
{
	int rows_A=A.getnRows(),   cols_A=A.getnCols(),
	    a_A=A.getRowIndexBase(), b_A=A.getColIndexBase();
	
	if((rows_A!=rows)||(cols_A!=cols)||(a_A!=a)||(b_A!=b)){
		
		std::cout << "\n"+test+": failed, dimensions or index range not equal.";
		return;
	}
	
	for(int i=0;i<a;i++)
	for(int j=0;j<b;j++) { 
		
		S err=relativeError(dptr[i][j],A(i+a,j+b),epsilon);
		if((err>precision)||(-err>precision)){
			
			std::cout << "\n"+test+": failed, "
			          << "component A["<<i<<"]["<<j<<"], relative error " << err;
		    return;
		}
	} // end for i

    std::cout << "\n"+test+": passed.";
} // end testEquals



// CONSTRUCTION - DESTRUCTION

/** Square matrix, all components intialized with zeroes.
 *  @param dim dimension (default 1).
 *  @param base index base, same for both indices (default 0).
 */
Matrix(int dim=1, int base=0) : 
rows(dim), cols(dim), a(base), b(base), dptr(new S*[dim])
{
	for(int i=0;i<dim;i++){
		
		dptr[i]=new S[dim];
		for(int j=0;j<dim;j++)dptr[i][j]=0;
	}
} // end constructor


/** Rectangular matrix, all components intialized with zeroes.
 *  @param nRows number of rows (default 1).
 *  @param nCols number of columns (default 1).
 *  @param row_base rwo index base (default 0).
 *  @param col_base column index base (default 0).
 */
Matrix(int nRows, int nCols, int row_base, int col_base) : 
rows(nRows), cols(nCols), a(row_base), b(col_base), dptr(new S*[nRows])
{
	for(int i=0;i<rows;i++){
	
		dptr[i]=new S[cols];
		for(int j=0;j<cols;j++)dptr[i][j]=0;	
	}
} // end constructor



/** Copy constructor
 */
Matrix(const Matrix& A) : 
rows(A.getnRows()), cols(A.getnCols()), 
a(A.getRowIndexBase()), b(A.getColIndexBase()), 
dptr(new S*[rows])
{
	for(int i=0;i<rows;i++){
		
		dptr[i]=new S[cols];
	    for(int j=0;j<cols;j++) dptr[i][j]=A(i+a,j+b);
	}
}


/** Construct from square data array, zero based indices only.
 * Cannot be extended to the nonsquare case as template parameters
 * apparently cannot be deduced.
 * @param A data array A[dim][dim] with entries of type S.
 * @param row_base row index base (default 0).
 * @param col_base column index base (default 0).
 */
template<int dim>
Matrix(S A[dim][dim], int row_base=0, int col_base=0) : 
rows(dim), cols(dim), a(row_base), b(col_base), dptr(new S*[rows])
{
	for(int i=0;i<rows;i++){
		
		dptr[i]=new S[cols];
	    for(int j=0;j<cols;j++) dptr[i][j]=A[i][j]; 
	}
}


/** Constructs a symmetric square matrix with upper half U
 * and same index base as U.
 * @param U upper half of matrix.
 */
Matrix(const UTRMatrix<S>& U) : 
rows(U.getDimension()), cols(U.getDimension()), 
a(U.getIndexBase()), b(U.getIndexBase()),
dptr(new S*[rows])
{
	S** u=U.getData();
	for(int i=0;i<rows;i++) dptr[i]=new S[cols];
	
	for(int i=0;i<rows;i++) 
	for(int j=i;j<cols;j++) dptr[i][j]=dptr[j][i]=u[i][j-i]; 
	
}


/** Assignement, dimensions of B must be less than or equal to 
 *  the dimensions of <code>this</code> (not checked), no new
 *  memory allocated.
 */
Matrix<S>& operator =(const Matrix<S>& B)
{
    rows=B.getnRows(); cols=B.getnCols();
	a=B.getRowIndexBase();    b=B.getColIndexBase();
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++) dptr[i][j]=B(i+a,j+b);
	
	return *this;
}
	

/** destructor */
~Matrix()
{
	for(int i=0;i<rows;i++){ delete[] dptr[i]; }
    delete[] dptr;
}


// NORMS

/** Square root of sum of entries squared. Type S must support an absolute value
 * <code>fabs(S& s)</code> mapping to Real.
 */
Real norm() const
{
	Real sum=0, absij;
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++) { absij=fabs(dptr[i][j]); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm


/** L2-norm of row i. Type S must have a function 
 * <code>fabs(S& s)</code> (absolute value mapping to Real).
 */
Real rowNorm(int i) const
{
	Real sum=0, absij;
	for(int j=0;j<cols;j++) { absij=fabs(dptr[i][j]); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm



/** L2-norm of column j. Type S must have a function 
 * <code>Real fabs(S s)</code> (absolute value mapping to Real).
 */
Real colNorm(int j)
{
	Real sum=0, absij;
	for(int i=0;i<rows;i++) { absij=fabs(dptr[i][j]); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm
	



// ALGEBRAIC OPERATIONS



/** Change into the transposed matrix
 */
void transposeSelf()
{
	// allocate new storage to preserve row major order
	S** A=dptr;
	S** At=new S*[cols];
	for(int j=0;j<cols;j++){
		
		  At[j]=new S[rows];
          for(int i=0;i<rows;i++)At[j][i]=A[i][j];
	}
		
	//free the old memory
    for(int i=0;i<rows;i++){ delete[] dptr[i]; }
    delete[] dptr;
	
	// switch row and column numbers
	int tmp=rows; rows=cols; cols=tmp;
	
	// reset
	dptr=At;
}



/** Scalar multiplication A -> fA
 */
Matrix& operator *= (Real f)
{
	// direct data access for speed
	S** A=dptr; 
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++) A[i][j]*=f;
		
	return *this;
} // end operator *=



/** Multiplication of row i by scalar f.
 *  Uses row indices relative to the base, not absolute indices.
 */
Matrix& scaleRow(int i, S f)
{
	// direct data access for speed
	S** A=dptr;
	for(int j=0;j<cols;j++) A[i-a][j]*=f;
		
	return *this;
} // end scaleRow


/** Multiplication of column j by scalar f.
 *  Uses column indices relative to the base, not absolute indices.
 */
Matrix& scaleCol(int j, S f)
{
	// direct data access for speed
	S** A=dptr;
	for(int i=0;i<rows;i++) A[i][j-b]*=f;
		
	return *this;
} // end scaleCol


/** Addition <code>this+=C</code>. 
 *  Matrices must be of equal dimensions, not checked.
 */
Matrix& operator += (Matrix& c)
{
	S** A=dptr; S** C=c.getData();
	for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++) A[i][j]+=C[i][j];
		
	return *this;
}



/** Right multiplication of <code>this</code> by C.
 *  Matrix dimensions must be compatible (not checked),
 *  index bases don't have to be equal, remain unchanged.
 */
Matrix& operator *=(const UTRMatrix<S>& c)
{
	// B is square cols by cols 
	// so multiplication with B leaves dimensions unchanged.
	Matrix<S> temp(rows,cols,0,0); 
	// direct data access for speed
	S** A=dptr;
	S** C=c.getData();
	S** tmp=temp.getData();
	
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++){
		
		S sum=A[i][0]*C[0][j];  // k=0 below
		for(int k=1;k<=j;k++)   sum+=A[i][k]*C[k][j-k];  
			
	    tmp[i][j]=sum;
	}
	
    for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++) A[i][j]=tmp[i][j];
		
	return *this;
} // end operator *=


/** RIGHT multiplication of <code>this</code> by transpose C'.
 *  Because of of row major allocation AC' is faster than AC. 
 *  Matrix dimensions must be compatible (not checked),
 *  index bases don't have to be equal, remain unchanged.
 */
Matrix& operator ^= (const UTRMatrix<S>& c)
{
	// B is square cols by cols 
	// so multiplication with B leaves dimensions unchanged.
	Matrix<S> temp(rows,cols,0,0); 
	// direct data access for speed
	S** A=dptr;
	S** C=c.getData();
	S** tmp=temp.getData();
	
	// recall that c is stored as its lower triangular transpose
	// and so c(i,j)=C[j][i]
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++){
		
		S sum=A[i][j]*C[j][0];   // k=j below
		for(int k=j+1;k<cols;k++) sum+=A[i][k]*C[j][k-j];  
			
	    tmp[i][j]=sum;
	}
	
    for(int i=0;i<rows;i++)
	for(int j=0;j<cols;j++) A[i][j]=tmp[i][j];
		
	return *this;
} // end operator ^=



/** Right multiplication of <code>this</code> by C.
 *  Matrix dimensions must be compatible (not checked),
 *  index bases don't have to be equal.
 */
Matrix& operator *=(const Matrix<S>& c)
{
    // we'll write B for the data of C, this*C is rows by cols_C
	int cols_C=c.getnCols();
	S** AC=new S*[rows];
	for(int i=0;i<rows;i++) AC[i]= new S[cols_C];
	// direct data access for speed
	S** A=dptr;
	S** C=c.getData();
	
	for(int i=0;i<rows;i++)
	for(int j=0;j<cols_C;j++){
		
		S sum=A[i][0]*C[0][j];  // k=0 below
		for(int k=1;k<cols;k++)   sum+=A[i][k]*C[k][j];  
			
	    AC[i][j]=sum;
	}
	
    // free the old memory
	for(int i=0;i<rows;i++) delete[] dptr[i]; delete[] dptr;
	// reset
	dptr=AC; cols=cols_C;
		
	return *this;
} // end operator *=



/** RIGHT multiplication of <code>this</code> by transpose C'.
 *  Because of of row major allocation AC' is faster than AC. 
 *  Matrix dimensions must be equal (not checked),
 *  index bases don't have to be equal.
 */
Matrix& operator ^= (const Matrix<S>& c)
{
	// this*C' is rows by rows_C
	int rows_C=c.getnRows();
	S** ACt=new S*[rows];
	for(int i=0;i<rows;i++) ACt[i]= new S[rows_C];
	// direct data access for speed
	S** A=dptr;
	S** C=c.getData();
	
    // AC' row dot row multip[lication
	for(int i=0;i<rows;i++)
	for(int j=0;j<rows_C;j++){
		
		S sum=A[i][0]*C[j][0];  // k=0
		for(int k=1;k<cols;k++) sum+=A[i][k]*C[j][k];  
			
	    ACt[i][j]=sum;
	}
	
    // free the old memory
	for(int i=0;i<rows;i++) delete[] dptr[i]; delete[] dptr;
	// reset
	dptr=ACt; cols=rows_C;
	
	return *this;
} // end operator ^=



// MATRIX FUNCTIONS: some analytic functions f(A), where A=this.

/** Matrix exponential <code>exp(A)</code>. 
 *  Matrix must be square, not checked. Same index bases as A.
 */
Matrix<S>& exp() const
{
	int dim=rows; 
	Real fn=norm(), f=4.0, g=1.0;
	
	int p=2;
	while((f<fn)){ f*=2.0; p++; }           // f=2^p
	g=1.0/f; g/=64.0;                       // g=1/2^{p+6}

	// temporary matrices, all index bases zero
	Matrix<S> U(dim);           // U=(A/2^{p+6})=(g*A) 
	Matrix<S> Ut(dim);          // U', we use ^=U ie. *=U'
	Matrix<S> SUk(dim);         // sum of U^k/k!, k=0,1,..,6

    S** A=dptr;
	for(int i=0;i<dim;i++)
	for(int j=0;j<dim;j++){ 
		
		U(i,j)=Ut(j,i)=SUk(i,j)=g*A[i][j];  // U,U'
		if(i==j) SUk(i,j)+=1.0;             // I+U
	}
	
	// Suk=I+U+U^2/2!+...+U^6/6! 
	for(int k=2;k<=8;k++){ U^=Ut; g=1.0/k; U*=g; SUk+=U; }
	
	// exp(A)=SUk^m, m=2^{p+6};
    Matrix<S>& F=*(new Matrix<S>(SUk));     
	for(int k=0;k<p+6;k++) F*=F;
		
    F.setRowIndexBase(a);
	F.setColIndexBase(b);
	return F;

} // end exp

		
	
}; // end Matrix




/** print rectangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const Matrix<S>& A)
{
	int rows=A.getnRows(), 
	    cols=A.getnCols();
	S** D=A.getData();
	os << endl << "Rectangular " << rows << " by " << cols << " matrix:"
	   << endl << endl;
	for(int i=0;i<rows;i++){
		
	    for(int j=0;j<cols-1;j++) os << D[i][j] << ", ";
		os << D[i][cols-1] << endl;
	}
    return os << endl << endl;
} // end operator <<





/***************************************************************************************
 *
 *            MATRIX ARRAYS
 *
***************************************************************************************/


/** Array of upper triangular matrices. Read only structure.
 *  Only the array of pointers to 
 *  the matrices is allocated. These pointers then have to be set to 
 *  actual matrices as needed.
 */
class UTRMatrixSequence {
	
	int n;                         // number of matrices
	// array of pointers to const UTRMatrix: *matrix[t] is read only
	Array1D<const UTRMatrix<Real>*> matrices;    
	
public:
	
	/** Constructor.
	 *
	 * @param m number of matrices.
	 */
	UTRMatrixSequence(int m) : 
	n(m), matrices(m) 
    {  }  
	
    /** Destructor. Deletes the matries pointed to. */
	~UTRMatrixSequence(){ for(int t=0;t<n;t++) delete matrices[t]; }

	
	/** Const reference to matrix[t].
	 */
    const UTRMatrix<Real>& getMatrix(int t) const { return *(matrices[t]); }
	
    /** Set pointer to matrix[t].
	 */
    void setMatrix(int t, const UTRMatrix<Real>& matrix) { matrices[t]=&matrix; }
	
}; // end UTRMatrixSequence




/** Array of rectangular matrices. Only the array of pointers to 
 *  the matrices is allocated. These pointers then have to be set to 
 *  actual matrices as needed.
 */
class MatrixSequence {
	
	int n;                                    // number of matrices
	Array1D<const Matrix<Real>*> matrices;    // matrices[t]: pointer to matrix(t)
	
public:
	
	/** Constructor.
	 *
	 * @param m number of matrices.
	 */
	MatrixSequence(int m) : 
	n(m), matrices(m) 
    {  }  
	
    /** Destructor. Deletes the matrices pointed to*/
	~MatrixSequence(){ for(int t=0;t<n;t++) delete matrices[t]; }

	
	/** Reference to matrix[t].
	 */
    const Matrix<Real>& getMatrix(int t) const { return *(matrices[t]); }
	
    /** Set pointer to matrix[t].
	 */
    void setMatrix(int t, const Matrix<Real>& matrix) { matrices[t]=&matrix; }
	
}; // endMatrixSequence




MTGL_END_NAMESPACE(Martingale)

#endif

