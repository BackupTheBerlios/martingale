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

#ifndef martingale_matrix_h
#define martingale_matrix_h


#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>                      // min,max
#include <cmath>
#include "TypedefsMacros.h"
//#include "tnt/tnt_array1d.h"
//#include "tnt/tnt_array2d.h"
#include "jama/jama_eig.h"
//#include "tnt_array1d_utils.h"
//#include "tnt_array2d_utils.h"
#include "Array.h"
//#include "Utils.h"                 

/*
 * DataStructures.h
 *
 * Created on February 19, 2003, 6:00 PM
 */

MTGL_BEGIN_NAMESPACE(Martingale)



/*! \file Matrix.h
 * <p>Real vectors, square lower triangular, square upper triangular, symmetric and rectangular
 * matrices stored in row major order. Barebones naive implementation with no checking of 
 * dimensional compatibility. All operations implemented as straightforward  
 * C-style loops with no optimizations. 
 * To check subscripting bounds #define SUBSCRIPT_CHECK.</p>
 *
 * <p>Matrix indices can be based on an arbitrary bases a,b. This means that the indices
 * used with the subscripting operators are 
 * \f[i=a,a+1,...,a+rows-1;\quad j=b,b+1,...,b+cols-1,\f]
 * where rows and cols are the number of rows and columns respectively.
 * This feature allows the use of the natural indices from the application domain
 * and will be referred to as "natural indexation" in each application. 
 * The subscripting operators take care of the necessary translations to the actual 
 * array indices. The index bases are irrelevant if the underlying data array is
 * accessed directly through the data pointers.</p>
 *
 * <p>For matrices allocated in row major order the loops for the product AB' have
 * a better memory access pattern (row dot row) than the loops for the product
 * AB (row dot column). There is very little difference in low dimensions but for 
 * 1000 by 1000 matrices AB' is already three times faster than AB.</p>
 *
 * <p>Thus we provide both the operator *= multiplying <code>this</code> on the right 
 * B (slower) and the operator ^= which multiplies on the right with B'.
 */



/***************************************************************************************

                      Useful global functions

***************************************************************************************/



/** The relative error in percent of an approximation to an exact value,
 *  positive if the approx > exact negative otherwise.
 *  Type S must support -(S&) and comparison ">".
 *
 * @param S real scalar type (float, double, long double).
 * @param exact true value.
 * @param approx approximation.
 * @param epsilon zero rounded up to epsilon to avoid division by zero.
 */
template<typename S>
S relativeError(S exact, S approx, S epsilon)
{
	S d=exact;
	if((-epsilon<d)&&(d<epsilon)) d=epsilon;
 
	return 100*(approx-exact)/d;
} 




/***************************************************************************************

                      C-Style vectors

***************************************************************************************/

// forward declarations
template<typename S,typename MatrixType> class Matrix;


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
class Vector {

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
 
      
   /** Subscripting relative to index base. */
   S& operator[](int i)
   {
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::checkSubscript(i,b,dim,"Vector");
	   #endif	   
	   return dptr[i-b]; 
    }
   
   /** Subscripting relative to index base. */
   const S& operator[](int i) const 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::checkSubscript(i,b,dim,"Vector");
	   #endif	  
	   return dptr[i-b]; 
   }
   
   
   // CONSTRUCTION - DESTRUCTION
   ~Vector(){ delete[] dptr; }
   
   /** Constructor initializes with zeroes
    */
   explicit Vector(int d, int base=0) : 
   dim(d), b(base), dptr(new S[d]) 
   { 
	   for(int i=0;i<dim;i++) dptr[i]=0; 
   }

   
   // shallow copy not good enough if we return from function
   Vector(const Vector& v) : 
   dim(v.getDimension()), b(v.getIndexBase()), dptr(new S[dim])
   {
	   S* vdptr=v.getData(); 
	   for(int i=0;i<dim;i++) dptr[i]=vdptr[i];
   }
   
   
  /** Initialize all  entries with the same value,
   *  index base is b=0.
   *
   * @param n dimension
   * @param v constant entry
   */
   Vector(int n, S& v) : 
   dim(n), b(0), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=v;
   }

   
  /** Construct from data array, index base is b=0.
   *
   * @param n dimension
   * @param a data array
   */
   Vector(int n, S* a) : dim(n), b(0), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=a[i];
   }
   
   
  /** Construct from data array, index base is b=0.
   * @param a data array
   */
   template<int n>
   Vector(S a[n]) : dim(n), dptr(new S[n])
   {
	    for(int i=0;i<dim;i++) dptr[i]=a[i];
   }

   
   /** Dimension of x must be less than or equal to the dimension
    *  of <code>this</code>, no new memory allocated.
    */
   Vector& operator=(const Vector& x)
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
   Vector& switchComponents(int i, int j)
   { 
       S tmp=dptr[i-b]; dptr[i-b]=dptr[j-b]; dptr[j-b]=tmp; 
   }
   
   
   /** dimension must be same (not checked), index bases need not be same,
    *  index base remains unchanged.
    */
   Vector& operator +=(const Vector& x)
   { 
	   S* xdptr=x.getData();
	   for(int i=0;i<dim;i++) dptr[i]+=xdptr[i];
	   return *this;
   }
   
   /** dimension must be same (not checked), index bases need not be same,
    *  index base remains unchanged.
    */
   Vector& operator -=(const Vector& x)
   { 
	   S* xdptr=x.getData();
	   for(int i=0;i<dim;i++) dptr[i]-=xdptr[i];
	   return *this;
   }
   
   /** division by integer (for use with class RandomObject)
    */
   Vector& operator /=(int N)
   { 
	   Real f=1.0/N;
	   for(int i=0;i<dim;i++) dptr[i]*=f;
	   return *this;
   }
   
   
   /** multiplication by scalar
    */
   Vector& operator *=(const S& lambda)
   { 
	   for(int i=0;i<dim;i++) dptr[i]*=lambda;
	   return *this;
   }

   
   /** Left multiplication by matrix A
    */
   template<typename MatrixBaseType>
   Vector& operator *=(const Matrix<S,MatrixBaseType>& A)
   { 
	   int d=A.rows();        // new dimension
	   S* Ax=new Real[d];         // new memory
	   for(int i=0;i<d;i++){
		   
		   S sum=0.0;
		   for(int j=A.rowBegin(i);j<A.rowEnd(i);j++) sum+=A.entry(i,j)*dptr[j];
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
	   S sum=0; 
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
void testEquals(const Vector<S>& v, S precision, S epsilon, string test) const
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


// PRINTING

std::ostream&  printSelf(std::ostream& os) const
{
	os << "\n\nVector of dimension " << dim << ":" << endl;
    for(int i=0;i<dim-1;i++) os << dptr[i] << ", ";
    os << dptr[dim-1];
    return os << endl;
} // end operator <<
	
	      
}; // end Vector


/** print Vector
 */
template<class S>
std::ostream& operator << (std::ostream& os, const Vector<S>& v)
{
	return v.printSelf(os);
} 


typedef Vector<Real> RealVector;
	


/***************************************************************************************
 *
 *                             MATRIX CLASSES
 *
***************************************************************************************/

/*! \file Matrix.h
 * <p>Real vectors, square lower triangular, square upper triangular and rectangular
 * matrices stored in row major order. Naive implementation. C-style loops with no optimizations. 
 * To check subscripting bounds #define SUBSCRIPT_CHECK.</p>
 *
 * <p>Matrix indices can be based on an arbitrary bases a,b. This means that the indices
 * used with the subscripting operators are 
 * \f[i=a,a+1,...,a+rows-1;\quad j=b,b+1,...,b+cols-1,\f]
 * where rows and cols are the number of rows and columns respectively.
 * This feature allows the use of the natural indices from the application domain
 * and will be referred to as "natural indexation" in each application. 
 *
 * <p>For matrices allocated in row major order the loops for the product AB' have
 * a better memory access pattern (row dot row) than the loops for the product
 * AB (row dot column). There is very little difference in low dimensions but for 
 * 1000 by 1000 matrices AB' is already three times faster than AB.</p>
 *
 * <p>Thus we provide both the operator *= multiplying <code>this</code> on the right 
 * B (slower) and the operator ^= which multiplies on the right with B'.
 * These operations are only implemented between matrices of the same type. 
 * Global functions implement binary operations between matrices of different type.
 *
 * <p><b>Symmetric matrices></b> The matrix classes which are implemented are closed 
 * under all algebraic operations. This is not true for symmetric matrices. The product 
 * of two symmetric matrices A,B is symmetric if and only if AB=BA. For this reason 
 * symmetric matrices are not implemented as a special class. 
 *
 * <p>If agebraic operations are needed symmetric matrices have to be treated as 
 * rectangular with redundant storage of entries. Upper triangular matrices stand in for
 * symmetric matrices in many cases such as such as eigenvector analysis, factor analysis 
 * and matrix functions. By convention the upper triangular matrix then represents the
 * symmetric matrix of which it is the upper half. In these cases the necessary algebra 
 * is carried out by the TNT classes.
 */

// IMPLEMENTATION NOTES: we use the Barton-Nachman trick.
// The type of matrix (triangular,symmetric,...) is determined by the base class
// from which the matrix class derives and this base class is a template parameter
// for the derived class template.
//
// Matrices are stored in row major order and usually traversed row by row.
// For each i we need to know where nonzero entries begin and end in this row.
// The functions rowBegin(i) and rowEnd(i) provide this information.
// The function rowEnd(i) returns the index one past the last nonzero entry
// in row i. The functions colBegin(i) and colEnd(i) do the same for
// the columns. The information returned by these functions is type specific
// (upper triangular, lower triangular,...) and not specific to each matrix 
// object and allows optimizations depending on the matrix type (matrix
// products).



/***************************************************************************************
 *
 *                      MATRIX BASE TYPES
 *
***************************************************************************************/

// forward declaration
template<typename S> class LowerTriangular;


// IMPLEMENTATION NOTES: the matrix base classes define the behaviour for the 
// various special matrix classes. They are put at the base of the derived class
// Matrix using the Barton-Nachman trick (base class is a template parameter for
// the derived template class). 
// 
// Consequently many functions will seem to have unnecessary general parameter
// signatures (the most general applying to any matrix type). This is necessary
// as the derived class Matrix has to be written in the most general case and
// has to call functions without any assumptions of the type of the base class.


/** Storage allocation and subscripting for square upper triangular matrices.
 *  @param S the type of matrix entries.
 */
template<typename S>
class UpperTriangular {
	
protected:
	
	int dim;      // number of rows = number of cols	    
	S** dptr;     // pointer to data
	
public:
	
	typedef LowerTriangular<S> TransposeType;
	
/** Number of rows. */	
int rows() const { return dim; }
/** Number of columns. */	
int cols() const { return dim; }
/** Empty function. Action needed only in nonsquare case.*/
void setCols(int) {   }

/** Memory allocation. 
 *  @param nRows number of rows.
 *  @param nCols number of columns.
 */
S** allocate(int nRows, int nCols)
{
	if(nRows!=nCols)
	{ cout<<"\n\nUpper triangular matrix must be square. Terminating."; exit(1); }
	S** data=new S*[nRows];
	for(int i=0;i<nRows;i++){
		
        data[i]=new S[nRows-i];
        for(int j=i;j<nRows;j++)data[i][j-i]=0;
    }
	return data;
}


/** Number of rows must equal number of columns.
 *  @param nRows number of rows.
 *  @param nCols number of columns. 
 */	
UpperTriangular(int nRows, int nCols) : dim(nRows), dptr(allocate(nRows,nCols)) {  } 
	
~UpperTriangular(){  for(int i=0;i<dim;i++) delete[] dptr[i]; delete[] dptr; }
	
	

/** Subscripting from zero index base: \f$i\leq j\f$. */
S& entry(int i, int j){ return dptr[i][j-i]; }
/** Subscripting from zero index base: \f$i\leq j\f$. */
const S& entry(int i, int j) const { return dptr[i][j-i]; }

		
/** Beginning of nonzero entries in row i relative to zero index base.*/
int rowBegin(int i) const { return i; }

/** End of nonzero entries in row i relative to zero index base.*/
int rowEnd(int i) const { return dim; }

/** Beginning of nonzero entries in column j relative to zero index base.*/
int colBegin(int j) const { return 0; }

/** End of nonzero entries in column j relative to zero index base.*/
int colEnd(int j) const { return j+1; }

	
}; // end UpperTriangular




/** Storage allocation and entries for square lower triangular matrices.
 *  @param S the type of matrix entries.
 */
template<typename S>
class LowerTriangular{
		
protected:
	
	int dim;      // number of rows = number of cols	    
	S** dptr;     // pointer to data
	
public:
	
	typedef UpperTriangular<S> TransposeType;
	
/** Number of rows. */	
int rows() const { return dim; }
/** Number of columns. */	
int cols() const { return dim; }
/** Empty function. Action needed only in nonsquare case.*/
void setCols(int) {   }


/** Memory allocation. 
 *  @param nRows number of rows.
 *  @param nCols number of columns.
 */
S** allocate(int nRows, int nCols)
{
	if(nRows!=nCols)
	{ cout<<"\n\nLower triangular matrix must be square. Terminating."; exit(1); }
	S** data=new S*[nRows];
	for(int i=0;i<nRows;i++){
		
        data[i]=new S[i+1];
        for(int j=0;j<=i;j++)data[i][j]=0;
   }
   return data;
}

		
/** Number of rows must equal number of columns.
 *  @param nRows number of rows.
 *  @param nCols number of columns. 
 */	
LowerTriangular(int nRows, int nCols) : dim(nRows), dptr(allocate(nRows,nCols)) {  } 

~LowerTriangular(){  for(int i=0;i<dim;i++) delete[] dptr[i]; delete[] dptr; }
	
	
/** Subscripting from zero index base: \f$j\leq i\f$. */
S& entry(int i, int j){ return dptr[i][j]; }
/** Subscripting from zero index base: \f$j\leq i\f$. */
const S& entry(int i, int j) const { return dptr[i][j]; }
	
		
/** Beginning of nonzero entries in row i relative to zero index base.*/
int rowBegin(int i) const { return 0; }

/** End of nonzero entries in row i relative to zero index base.*/
int rowEnd(int i) const { return i+1; }

/** Beginning of nonzero entries in column j relative to zero index base.*/
int colBegin(int j) const { return j; }

/** End of nonzero entries in column j relative to zero index base.*/
int colEnd(int j) const { return dim; }

		
}; // end LowerTriangular




/** Storage allocation and entries for rectangular matrices.
 *  @param S the type of matrix entries.
 */
template<typename S>
class Rectangular {
	
protected:
	
	int rows_,        // number of rows 
	    cols_;        // number of columns
	S** dptr;         // pointer to data

	
public:
	
	typedef Rectangular<S> TransposeType;
	
/** Number of rows. */	
int rows() const { return rows_; }
/** Number of columns. */	
int cols() const { return cols_; }
/** Set field for number of columns. No action on memory.*/
void setCols(int nCols) { cols_=nCols;  }


/** Memory allocation. 
 *  @param nRows number of rows.
 *  @param nCols number of columns.
 */
S** allocate(int nRows, int nCols)
{
	S** data=new S*[nRows];
	for(int i=0;i<nRows;i++){
		
        data[i]=new S[nCols];
        for(int j=0;j<nCols;j++)data[i][j]=0;
   }
   return data;
}
	
		
/** 
 *  @param nRows number of rows.
 *  @param nCols number of columns. 
 */	
Rectangular(int nRows, int nCols) : 
rows_(nRows), cols_(nCols), dptr(allocate(nRows,nCols))
{    } 

~Rectangular(){  for(int i=0;i<rows_;i++) delete[] dptr[i]; delete[] dptr; }
	
	
/** Subscripting from zero index base. */
S& entry(int i, int j){ return dptr[i][j]; }
/** Subscripting from zero index base. */
const S& entry(int i, int j) const { return dptr[i][j]; }

	
/** Beginning of nonzero entries in row i relative to zero index base.*/
int rowBegin(int i) const { return 0; }

/** End of nonzero entries in row i relative to zero index base.*/
int rowEnd(int i) const { return cols_; }

/** Beginning of nonzero entries in column j relative to zero index base.*/
int colBegin(int j) const { return 0; }

/** End of nonzero entries in column j relative to zero index base.*/
int colEnd(int j) const { return rows_; }

	
}; // end Rectangular



/***************************************************************************************
 *
 *                      MATRIX PRODUCT RETURN TYPE
 *
***************************************************************************************/


/** In general the return type is the most general: "Rectangular" */
template<typename S, typename MatrixType1, typename MatrixType2>
struct ProductType { typedef Rectangular<S> type; };

/** If both factors have the same type then the product has this type also.
 *  So we specialize:
 */
template<typename S, typename MatrixType>
struct ProductType<S,MatrixType,MatrixType> { typedef MatrixType type; };
	



/***************************************************************************************
 *
 *                             MATRIX
 *
***************************************************************************************/


/** <p>Square lower triangular, upper triangular, symmetric or general rectangular matrix 
 *  with indices based on arbitrary bases a,b stored in row major order. 
 *  This means that the indices used with the subscripting operators are 
 *  \f[i=a,a+1,...,a+rows-1;\quad j=b,b+1,...,b+cols-1,\f]
 *  where rows and cols are the number of rows and columns respectively.
 *
 * @param S type of matrix entries.
 * @param MatrixBaseType: UpperTriangular<S>, LowerTriangular<S>, Symmetric<S> or Rectangular<S>.
 */
template<typename S, typename MatrixBaseType>
class Matrix : public MatrixBaseType {

protected:

	int a;           // row index base
	int b;           // column index base
	
	static S cache[SMALL][SMALL];          // workspace for small matrix optimization.


public:
	 
	/** The type of the matrix: UpperTriangular<S>, LowerTriangular<S>,...*/
	typedef typename MatrixBaseType::TransposeType TransposeType;


// FREE MEMORY
	
void deallocate(){ for(int i=0;i<rows();i++) delete[] dptr[i]; delete[] dptr; }
		
// ACCESSORS

/** Pointer to data. */
S** getData() const { return dptr; }

/** Row index base a, row indices i=a,a+1,...,a+rows-1 */
int getRowIndexBase() const { return a; }
/** Row index base set equal to r. */
void setRowIndexBase(int r){ a=r; }

/** Column index base b, column indices j=b,b+1,...,b+cols-1 */
int getColIndexBase() const { return b; }
/** Column index base set equal to r. */
void setColIndexBase(int c){ b=c; }


/** Subscripting relative to the index bases a,b, that is, we use indices
 * \f[i=a,a+1,\dots,a+rows-1;\quad j=b,b+1,\dots,b+cols-1.\f]
 */
S& operator()(int i, int j)
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,a,b,rows(),cols(),"Matrix");
	#endif		
	return entry(i-a,j-b);
}


/** Subscripting relative to the index bases a,b, that is, we use indices
 * \f[i=a,a+1,\dots,a+rows-1;\quad j=b,b+1,\dots,b+cols-1.\f]
 */
const S& operator()(int i, int j) const 
{ 
	#ifdef SUBSCRIPT_CHECK
	  SubscriptCheck::checkSubscript(i,j,a,b,rows(),cols(),"Matrix");
	#endif		
	return entry(i-a,j-b);
}	
	

// EQUALITY CHECK

/** Test for entry by entry equality. Matrix types must be equal.
 *  Type S must support absolute value and comparison ">".
 *  Equality is defined by an upper bound on the acceptable relative 
 *  error in percent.
 *  
 * @param B matrix to compare this with.
 * @param precision upper bound on the acceptable relative error in percent.
 * @param epsilon zero denominator reset to epsilon.
 * @param test string printed to identify test.
 */
void testEquals(const Matrix& B, S precision, S epsilon, string test) const
{
	int a_B=B.getRowIndexBase(), b_B=B.getColIndexBase();
	if((B.rows()!=rows())||(B.cols()!=cols())){
		
		std::cout << "\nMatrix::testEquals(): "+test+": failed, dimensions not equal.";
		return;
	}
	
	if((a_B!=a)||(b_B!=b)){
		
		std::cout << "\nMatrix::testEquals(): "+test+": failed, index bases not equal.";
		return;
	}
	
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) { 
		
		    S err=relativeError(entry(i,j),B.entry(i,j),epsilon);
		    if((err>precision)||(-err>precision)){
			
			    std::cout << "\n"+test+": failed, "
			              << "component B["<<i<<"]["<<j<<"], relative error " << err;
		        return;
		    } // endif
	} // end for i

    std::cout << "\n"+test+": passed.";
} // end testEquals


// CONSTRUCTION - DESTRUCTION


/** Constructor for square matrix with both row and column index base = a, 
 *  all components intialized with zeroes.
 *  @param dim square matrix: rows=cols=dim.
 *  @param c index base (i,j=c,...,c+dim-1), default 0.
 */
Matrix(int dim, int c=0) : MatrixBaseType(dim,dim), a(c), b(c) {    }


/** Rectangular matrix, all components intialized with zeroes.
 *  @param nRows number of rows.
 *  @param nCols number of columns.
 *  @param row_base row index base.
 *  @param col_base column index base.
 */
Matrix(int nRows, int nCols, int row_base, int col_base) : 
MatrixBaseType(nRows,nCols),
a(row_base), b(col_base) 
{    }



/** Copy constructor, matrix types must be equal.
 */
Matrix(const Matrix& A) : MatrixBaseType(A.rows(),A.cols()),
a(A.getRowIndexBase()), b(A.getColIndexBase()) 
{
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)=A.entry(i,j);
}


/** Construct from square data array, zero based indices only.
 * Cannot be extended to the nonsquare case as template parameters
 * cannot be deduced. All matrix types, only the relevant array entries are
 * used.
 * @param A data array A[dim][dim] with entries of type S.
 * @param row_base row index base (default 0).
 * @param col_base column index base (default 0).
 */
template<int dim>
Matrix(S A[dim][dim], int row_base=0, int col_base=0) : 
MatrixBaseType(dim,dim),
a(row_base), b(col_base) 
{
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)=A[i][j]; 
}


/** Assignement, matrix types and sizes must be equal.
 */
Matrix& operator =(const Matrix& B)
{
    if(this=&B) return *this;
	if((rows()!=B.rows())||(cols()!=B.cols())){
		
		cout << "\n\nMatrix assignement: matrix sizes different. Terminating.";
		exit(1);
	}

	a=B.getRowIndexBase(); b=B.getColIndexBase();
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)=B.entry(i,j);
	
	return *this;
}



// NORMS

/** Square root of sum of absolute values of entries squared. 
 *  Type S must support absolute value <code>abs(S& s)</code> mapping to Real.
 */
Real norm() const
{
	Real sum=0.0, absij;
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) { absij=std::abs(entry(i,j)); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm



/** L2-norm of row i. Type S must have a function 
 * <code>abs(S& s)</code> (absolute value mapping to Real).
 */
Real rowNorm(int i) const
{
	Real sum=0, absij;
	for(int j=rowBegin(i);j<rowEnd(i);j++) { absij=std::abs(entry(i,j)); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm



/** L2-norm of column j. Type S must have a function 
 * <code>Real fabs(S s)</code> (absolute value mapping to Real).
 */
Real colNorm(int j)
{
	Real sum=0, absij;
	for(int i=colBegin(j);i<colEnd(j);i++) { absij=std::abs(entry(i,j)); sum+=absij*absij; }
	
	return sqrt(sum);
} // end norm
	


/** Multiplication of row i by scalar f.
 *  Uses row index i relative to the base.
 */
Matrix& scaleRow(int i, Real f)
{
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i-a,j)*=f; 
	return *this;
} // end scaleRow


/** Multiplication of column j by scalar f.
 *  Uses column index j relative to the base.
 */
Matrix& scaleCol(int j, Real f)
{
	for(int i=colBegin(j);i<colEnd(j);i++) entry(i,j-b)*=f;		
	return *this;
} // end scaleRow


	

// ALGEBRAIC OPERATIONS

/** The quadratic form \f$(Cx,x)\f$ where C is the symmetric matrix of 
 *  which this is the upper half.
 */
S quadraticForm(const Vector<S>& x) const
{
	int dim_x=x.getDimension(), 
	    b_x=x.getIndexBase();
	if((rows()!=dim_x)||(cols()!=dim_x)){
		
		cout << "\n\nMatrix::quadraticForm(): dimensional mismatch. Terminating.";
		exit(1);
	}
	S d=0, u=0, x_i,x_j;
	// diagonal
	for(int i=0;i<rows();i++) { x_i=x[i+b_x]; d+=x_i*x_i*entry(i,i); }
	// above the diagonal
	for(int i=0;i<rows();i++)
	for(int j=i+1;j<cols();j++) { x_i=x[i+b_x]; x_j=x[j+b_x]; u+=x_i*x_j*entry(i,j); }

	return d+2*u;
}



/** Transpose. Must allocate new memory to preserve row major order.
 */
Matrix<S,TransposeType>& transpose() const
{
	Matrix<S,TransposeType>* 
	trnsps=new Matrix<S,TransposeType>(cols(),rows(),b,a);
    for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) trnsps->entry(j,i)=entry(i,j);
	
	return *trnsps;
}
	

/** Matrix addition. Dimensions and types must match (not checked), 
 *  index bases need not.
 */
Matrix& operator += (const Matrix& B)
{
	if((rows()!=B.rows())||(cols()!=B.cols())){
		
		cout << "\n\nMatrix +=: matrix sizes different. Terminating.";
		exit(1);
	}
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)+=B.entry(i,j);
		
	return *this;
} // end operator +=


/** Multiplication of matrix by a Real scalar
 */
Matrix& operator *= (Real f)
{
	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)*=f;
		
	return *this;
} // end operator +=


/** Multiplication on the right by the matrix B of the same type.
 *  Dimensional compatibility not checked.
 */
// product will have the same type as this.
Matrix& operator *= (const Matrix& B)
{
	if(cols()!=B.rows()) {
		
		cout << "\n\nMatrix *=: matrix sizes not compatible. Terminating.";
		exit(1);
	}
	
	// if both matrices are small and square use the workspace, no memory allocation.
	// matrix size does not change and so the product matrix can use the old memory.
	if((rows()<SMALL)&&(rows()==cols())&&(cols()==B.cols())){
	
	    for(int i=0;i<rows();i++)
	    for(int j=rowBegin(i);j<rowEnd(i);j++){
		
            S sum=0.0;
		    int u=std::max(rowBegin(i),B.colBegin(j)),
		        v=std::min(rowEnd(i),B.colEnd(j));
		    for(int k=u;k<v;k++) sum+=entry(i,k)*B.entry(k,j);
		    cache[i][j]=sum;		
	    }
	    // copy the entries from the workspace
	    for(int i=0;i<rows();i++)
	    for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)=cache[i][j];
		return *this;
	} // end if
	
	// otherwise allocate new memory for the product
	Matrix<S,MatrixBaseType>* product = new Matrix<S,MatrixBaseType>(rows(),B.cols(),a,b);
    for(int i=0;i<rows();i++)
	for(int j=product->rowBegin(i);j<product->rowEnd(i);j++){
		
        S sum=0.0;
		int u=std::max(rowBegin(i),B.colBegin(j)),
		    v=std::min(rowEnd(i),B.colEnd(j));
		for(int k=u;k<v;k++) sum+=entry(i,k)*B.entry(k,j);
		product->entry(i,j)=sum;
	}
	deallocate();                // free old memory
	dptr=product->getData();
	setCols(B.cols());
	return *this;
		
} // end operator *=



/** <p>RIGHT multiplication of <code>this</code> by the transpose B'. 
 *  Because of row major allocation the matrix product AB' is faster
 *  than AB. Dimensions and types of A,B' must match (not checked), 
 *  bases need not.</p>
 */
// product will have the same type as this.
Matrix& operator ^= (const Matrix<S,TransposeType>& B)
{
	 if(cols()!=B.cols()){
		   
		   cout << "\n\nMatrix ^=: dimensions not compatible. Terminating.";
	       exit(1);
	 }

	// if both matrices are small and square use the workspace, no memory allocation.
	// matrix size does not change and so the product matrix can use the memory of dptr
	if((rows()<SMALL)&&(rows()==cols())&&(rows()==B.rows())){
	
	    for(int i=0;i<rows();i++)
	    for(int j=rowBegin(i);j<rowEnd(i);j++){
		
            Real sum=0.0;
		    int u=std::max(rowBegin(i),B.rowBegin(j)),
		        v=std::min(rowEnd(i),B.rowEnd(j));
		    for(int k=u;k<v;k++) sum+=entry(i,k)*B.entry(j,k);
		    cache[i][j]=sum;		
	    }
	    // copy the entries from the workspace
	    for(int i=0;i<rows();i++)
	    for(int j=rowBegin(i);j<rowEnd(i);j++) entry(i,j)=cache[i][j];
		return *this;
	} // end if
	
	// otherwise allocate new memory for the product,
	Matrix<S,MatrixBaseType>* product = new Matrix<S,MatrixBaseType>(rows(),B.rows(),a,b);
    for(int i=0;i<rows();i++)
	for(int j=product->rowBegin(i);j<product->rowEnd(i);j++){
		
        S sum=0.0;
		int u=std::max(rowBegin(i),B.rowBegin(j)),
		    v=std::min(rowEnd(i),B.rowEnd(j));
		for(int k=u;k<v;k++) sum+=entry(i,k)*B.entry(j,k);
		product->entry(i,j)=sum;
	}
	deallocate();                     // free the old memory
	dptr=product->getData();
	setCols(B.rows());
	return *this;
		
} // end operator ^=




// MATRIX FUNCTIONS: some analytic functions f(A), where A=this.


/** The upper half of the product <code>C=RR'</code>, where <code>R=this</code>.
 *  Zero based indexation.
 */
Matrix< S,UpperTriangular<S> >& aat() const
{
    Matrix< S,UpperTriangular<S> >* 
	c = new Matrix< S,UpperTriangular<S> >(rows(),rows(),0,0); // aat
	for(int i=0;i<rows();i++)
	for(int j=i;j<rows();j++){
		
		S sum=0.0;                                      // row_i\cdot row_j
        int u=std::max(rowBegin(i),rowBegin(j)),
		    v=std::min(rowEnd(i),rowEnd(j));
		for(int k=u;k<v;k++) sum+=entry(i,k)*entry(j,k);
		c->entry(i,j)=sum;
	}
	return *c;
} // end aat



/** Matrix exponential <code>exp(A)</code>. 
 */
Matrix& exp() const
{
    if(rows()!=cols()){
		
		cout << "\n\nMatrix::exp(): matrix not square. Terminating.";
		exit(1);
	}
	
	int n=2;
	Real fn=norm(), f=4.0, g=1.0;
	while((f<fn)){ f*=2.0; n++; }                   // f=2^n
	g=1.0/f; g/=64.0;                               // g=1/2^{n+6}

	// temporary matrices
	Matrix<S,MatrixBaseType> U(rows());             // U=(A/2^{n+6})=g*A 
	Matrix<S,TransposeType>  Ut(rows());            // U', we use ^=U ie. *=U'
	Matrix<S,MatrixBaseType> SUk(rows());           // sum of U^k/k!, k=0,1,..,6

	for(int i=0;i<rows();i++)
	for(int j=rowBegin(i);j<rowEnd(i);j++){ 
		
		U(i,j)=Ut(j,i)=SUk(i,j)=g*entry(i,j);       // U=A/2^{n+6}
		if(i==j) SUk(i,j)+=1.0;                     // I+U
	}
	
	// Suk=I+U+U^2/2!+...+U^8/8! 
	for(int k=2;k<=8;k++){ U^=Ut; g=1.0/k; U*=g; SUk+=U; }
	
	// exp(A)=Suk^m, m=2^{n+6};
    Matrix<S,MatrixBaseType>& F=*(new Matrix<S,MatrixBaseType>(SUk));    
	for(int k=0;k<n+6;k++) F*=F;
	F.setRowIndexBase(a);	
	F.setColIndexBase(b);	
	return F;

} // end exp


// PRINTING

std::ostream& printSelf(std::ostream& os) const
{
	os << "\n\nMatrix, rows: " << rows() << endl << endl;
	for(int i=0;i<rows();i++){
		
		int u=rowBegin(i), v=rowEnd(i); 
		for(int j=0;j<u;j++) os << 0.0 << ", ";
		for(int j=u;j<v-1;j++) os << entry(i,j) << ", ";
		os << entry(i,v-1) << endl;
	}
    return os << endl << endl;
} 


}; // end Matrix




// cache must be initialized outside the class:
template<typename S,typename MatrixBaseType>
S Matrix<S,MatrixBaseType>::cache[SMALL][SMALL] = { };


/** print matrix
 */
template<typename S, typename MatrixBaseType>
std::ostream& operator << (std::ostream& os, const Matrix<S,MatrixBaseType>& A) 
{
    return A.printSelf(os);
} // end operator <<





// SOME  MATRIX TYPES

// since there are no templated typedefs:
template<typename S>
class UTRMatrix : public Matrix< S,UpperTriangular<S> >{
	
public:

/** @param S type of entries.
 *  @param dim number of rows = number of columns.
 *  @param a index base (for both rows and columns).
 */
UTRMatrix(int dim, int a=0) : Matrix< S,UpperTriangular<S> >(dim,dim,a,a) {  }
	
};
	
	
// upper triangular and symmetric real types are more elaborate below.
typedef Matrix< Real,Rectangular<Real> > RealMatrix;
typedef Matrix< Real,LowerTriangular<Real> > LTRRealMatrix;



/***************************************************************************************
 *
 *                             UPPER TRIANGULAR MATRICES
 *
***************************************************************************************/



// DIAGONALIZATION, RANK REDUCED FACTORIZATION, EIGEN ANALYSIS

// forward declarations
class UTRRealMatrix;


/** JAMA object containing eigenvalues (sorted ascending) and ON matrix of associated 
 *  eigenvectors of the symmetric matrix of which C is the upper half.
 */
JAMA::Eigenvalue<Real>* eigenDecomposition(const UTRRealMatrix& C);


/** <p>Let D be the symmetric matrix with upper half C. This function computes 
 *  the matrix R which best approximates D as a product RR' and has rank r. The matrix R
 *  is computed by diagonalizing D, setting all but the r largest eigenvalues
 *  equal to zero and taking the square root of the remaining eigenvalues.
 *  Template parameter S is the type of the matrix entries and
 *  must be a scalar type (float, double, long double).
 *
 * <p>Maintains the row index base of C, columns indexed from zero.
 */
RealMatrix& rank_Reduced_Root(const UTRRealMatrix& C, int r);
		


/** C is interpreted as a multinormal covariance matrix.
 *  This method prints how much variability is captured by the 5 largest
 *  eigenvalues of C. 
 */
void factorAnalysis(const UTRRealMatrix& C, int r);
	


/** Let D be the symmetric matrix with upper half C. D must be positive semidefinite.
 *  Computes the best approximate factorization \f$C\simeq RR'\f$, where R has rank r
 *  and returns the relative error in the trace norm
 *  \f[||A||=Tr(A'A)=\sum\nolimits_{ij}A^2_{ij}.\f]
 *
 * @param message short description of test.
 */
void factorizationTest(const UTRRealMatrix& C, int r, string message="");



/** <p>Real upper triangular matrix. Index base for rows and columns must be equal.
 *  To make the subscripting operator as fast as possible it is limited to
 *  access elements on or above the diagonal (where the marix is known to be zero).
 *  Out of bounds error below the diagonal 
 *
 *  <p>Subscripting speed is crucial. For this reason upper triangular matrices
 *  often stand in for symmetric matrices which have a slower subscripting operator.
 *  The upper triangular matrix then represents the symmetric matrix of which it is 
 *  the upper half. This expalains why we implement operations for upper triangular
 *  matrics which make sense for symmetric matrices (eigen value analysis, 
 *  pseudo square roots,...).
 */
class UTRRealMatrix : public Matrix< Real,UpperTriangular<Real> > {

public:
	
/** Constructor, initializes all components with zeroes.<br>
 *  @param d dimension of square matrix 
 *  @param b index base: b<=i,j<b+d
 */
UTRRealMatrix(int d, int b=0) : Matrix< Real,UpperTriangular<Real> >(d,d,b,b) {  }

	
/** Construct from data array, only upper triangular half is used.<br>
 * @param A data array
 * @param b index base (default 0).
 */
template<int n>
UTRRealMatrix(Real A[n][n], int b=0) : Matrix< Real,UpperTriangular<Real> >(A,b) {   }
	
	
/** The symmetric matrix of which this is the upper half. 
 *  Return by value is deliberate.
 */
RealMatrix& symmetricCompletion() const;


/** The matrix inverse. This is upper triangular of the same size.
 */
UTRRealMatrix& inverse() const;
	

/** Cholesky root L (lower triangular, LL'=A where A is the symmetric matrix  
 *  with upper half <code>this</code>). 
 *  Terminates with error message if this is not positive definite. 
 */
LTRRealMatrix& ltrRoot() const;


/** Upper triangular root U of A, where A is the symmetric matrix of which
 *  <code>this</code>) is the upper half. Satisfies UU'=A. Differs from the 
 *  Cholesky root as it is upper triangular instead of lower triangular.
 *  Terminates with error message if this is not positive definite. 
 */
UTRRealMatrix& utrRoot() const;


/** <p>Let C be the symmetric matrix with upper half <code>this</code>. 
 *  This function computes the matrix R which best approximates D as a 
 *  product RR' and has rank r. The matrix R is computed by 
 *  diagonalizing D, setting all but the r largest eigenvalues
 *  equal to zero and taking the square root of the remaining eigenvalues.
 *
 * <p>The row index base remains the same, columns of the root are indexed from zero.
 */
RealMatrix& rankReducedRoot(int r) const;


/** Matrix exponential. 
 */
UTRRealMatrix& exp() const;


/** Computes the function f(A) of the symmetric matrix A of which <code>this</code> 
 *  is the upper half. The matrix funtion f(A) is computed as follows: diagonalize A as
 *  \f[A=U D U',\quad\hbox{where}\quad D=diag(\lambda_j)\f]
 *  is the diagonal matrix with the eigenvalues of A on the diagonal the columns of U 
 *  are associated eigenvectors. Set \f$f(D)=diag(f(\lambda_j))\f$ (f is applied to 
 *  each eigenvalue). Then
 *  \f[f(A)=Uf(D)U'.\f]
 */ 
RealMatrix& matrixFunction(Real (*f)(Real));


/** Matrix is interpreted as the upper half of a multinormal covariance matrix.
 *  This method prints how much variability is captured by the r largest
 *  eigenvalues of C. 
 */
void analyseFactors(int r) const;


/** Let C be the symmetric matrix with upper half <code>this</code>. 
 *  C must be positive semidefinite. Computes the best approximate factorization 
 *  \f$C\simeq RR'\f$, where R has rank r and returns the relative error in the trace 
 *  norm
 *  \f[||A||=Tr(A'A)=\sum\nolimits_{ij}A^2_{ij}.\f]
 *
 * @param message put test in context.
 */
void testFactorization(int r, string message="") const;


}; // end UTRRealMatrix




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
	Array1D<const UTRRealMatrix*> matrices;    
	
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
    const UTRRealMatrix& getMatrix(int t) const { return *(matrices[t]); }
	
    /** Set pointer to matrix[t].
	 */
    void setMatrix(int t, const UTRRealMatrix& matrix) { matrices[t]=&matrix; }
	
}; // end UTRMatrixSequence




/** Array of rectangular matrices. Only the array of pointers to 
 *  the matrices is allocated. These pointers then have to be set to 
 *  actual matrices as needed.
 */
class MatrixSequence {
	
	int n;                                    // number of matrices
	Array1D<const RealMatrix*> matrices;      // matrices[t]: pointer to matrix(t)
	
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
    const RealMatrix& getMatrix(int t) const { return *(matrices[t]); }
	
    /** Set pointer to matrix[t].
	 */
    void setMatrix(int t, const RealMatrix& matrix) { matrices[t]=&matrix; }
	
}; // endMatrixSequence




MTGL_END_NAMESPACE(Martingale)

#endif

