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

#ifndef martingale_array_h    
#define martingale_array_h
// #define SUBSCRIPT_CHECK              // enable to check for subscript out of range error                     

#include "TypedefsMacros.h"
#include <cstdlib>                       // for exit()
#include <string>

MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file Array.h
 *  <p> Barebones array classes in dimension 1,2,3.
 *  Provides subscripting from arbitrary index base to enable 
 *  indexing with indices which are natural for the intended application.
 *  <a name="index-base"></a>
 *  If \f$n_j\f$ denotes the size of the array in dimension j then indexing from base
 *  \f$b_j\f$ means that the subscripting operator uses indices
 *  \f[k=b_j,b_j+1,\dots,b_j+n_j-1\f]
 *  in dimension j. The array elements are of type S and the default constructor for 
 *  S is used to determine the size of these elements. On construction the array elements 
 *  are initialized with the integer 0. Consequently 0 must be convertible to type S. This 
 *  is useful if the array is populated with pointers where 0 signifies the null pointer.
 *
 * <p>To turn on subscript checking for out of range error #define SUBSCRIPT_CHECK.
 *  Note that this can slow down code using subscripting in innermost nested loops
 *  by a factor of ten or more.
 *
 * @param S type of array elements.
 */
 


/***************************************************************************************
 *
 *           SUBSCRIPT CHECKING
 *
 **************************************************************************************/

// linker errors if this is a namespace instead of a struct. Why????
struct SubscriptCheck {

/** Checks if \f$i\in[b_1,b_1+n_1)\f$. */
static void checkSubscript(int i, int b1, int n1, string str)
{
	if((i<b1)||(i>b1+n1-1)){
		
	    cout << "\n\nSubscript out of range: " << str 
	         << "\ni = " << i << " not in [" << b1 << ", " << b1+n1-1 << "]"
		     << "\n index base = " << b1 << ", dimension = " << n1
	         << "\nTerminating.";
	    exit(0);
	}
} // end checkSubscript


/** Checks if \f$i\in[b_1,b_1+n_1)\f$ and \f$j\in[b_2,b_2+n_2)\f$. */
static void checkSubscript(int i, int j, int b1, int b2, int n1, int n2, string str)
{
	if((i<b1)||(i>b1+n1-1)){ 
		
	    cout << "\n\nSubscript out of range: " << str 
	         << "\ni = " << i << " not in [" << b1 << ", " << b1+n1-1 << "]"
		     << "\n index base = " << b1 << ", rows = " << n1
		     << "\nTerminating.";
		exit(0);
	}
		
	
	if((j<b2)||(j>b2+n2-1)){
		
	    cout << "\n\n Subscript out of range: " << str 
	         << "\nj = " << j << " not in [" << b2 << ", " << b2+n2-1 << "]"
		     << "\n index base = " << b2 << ", cols = " << n2
	         << "\nTerminating.";
	    exit(0);
	}
} // end checkSubscript



/** Checks if \f$i\in[b_1,b_1+n_1)\f$, \f$j\in[b_2,b_2+n_2)\f$ and
 *  \f$k\in[b_3,b_3+n_3)\f$.
 */
static void checkSubscript
(int i, int j, int k, int b1, int b2, int b3, int n1, int n2, int n3, string str)
{
	if((i<b1)||(i>b1+n1-1)){ 
		
	    cout << "\n\nSubscript out of range: " << str 
	         << "\ni = " << i << " not in [" << b1 << ", " << b1+n1-1 << "]"
		     << "\nTerminating.";
		exit(0);
	}
		
	
	if((j<b2)||(j>b2+n2-1)){
		
	    cout << "\n\n Subscript out of range: " << str 
	         << "\nj = " << j << " not in [" << b2 << ", " << b2+n2-1 << "]"
	         << "\nTerminating.";
	    exit(0);
	}
	
	
	if((k<b3)||(k>b3+n3-1)){
		
	    cout << "\n\n Subscript out of range: " << str 
	         << "\nj = " << k << " not in [" << b3 << ", " << b3+n3-1 << "]"
	         << "\nTerminating.";
	    exit(0);
	}
} // end checkSubscript


}; // end SubscriptCheck




/**********************************************************************************
 *
 *         1 DIMENSIONAL   ARRAY
 *
 *********************************************************************************/


/** Lightweight array in dimension 1. 
 *
 * @param S type of array elements.
 */
template<typename S>
class Array1D {
	
protected:
	
	/** Index base. */
	int b;
	/** Number of elements. */
	int n;         
				
	/** Data array. */
	S* dptr;

	
public:
	

// ACCESSORS

	/** <a href="index-base">Index base</a>. 
	 */
    int getIndexBase() const { return b; }
	
	void setIndexBase(int base){ b=base; }
		
    /** Number of array elements. */
	int getDimension() const { return n; }
	
	/** Resets dimension, NO EFFECT ON MEMORY. */
	void setDimension(int q) { n=q; }
	
	/** Pointer to data array */
	S* getData() const { return dptr; }

	
// CONSTRUCTOR
	
	/** @param n_ number of array elements.
	 *  @param b_ <a href="index-base">index base</a>.
	 */
	explicit Array1D(int n_, int b_=0) : 
	b(b_), n(n_) 
	{  
         dptr=new S[n];
		 for(int i=0;i<n;i++) dptr[i]=0;
	} 
	
	
	Array1D(const Array1D& x) : 
	b(x.getIndexBase()), 
	n(x.getDimension()) 
	{  
         dptr=new S[n];
		 S* xdptr=x.getData();
		 for(int i=0;i<n;i++) dptr[i]=xdptr[i];
	} 

	
	~Array1D(){ delete[] dptr; }
	
   /** Subscripting.*/
   const S& operator[](int i) const 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,b,n,"Array1D");
	   #endif	
	   return dptr[i-b]; 
   }
   
   /** Subscripting.*/
   S& operator[](int i)  
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,b,n,"Array1D");
	   #endif	
	   return dptr[i-b];
   }
   
   /** Copies the elements of b which must be of the same
    *  length as a (not checked).
    */
   Array1D<S>& operator=(const Array1D<S>& B)
   {
	   S* dptrB=B.getData();
	   for(int j=0;j<n;j++) dptr[j]=dptrB[j];
	   b=B.getIndexBase();  
	   
	   return *this;
   }
   
// EUCLLIDEAN NORM, SCALING
   
    /** Euclidean norm. */     
    Real norm() const
    {
        Real u=0;
        for(int i=0;i<n;i++) u+=dptr[i]*dptr[i];
        return sqrt(u);
    }
	
	/** Multiplication by a scalar, returns scaled array. */
	Array1D& scale(S f)
    {
        for(int i=0;i<n;i++) dptr[i]*=f;
        return *this;
    }  
	
	/** Dot product with the array X of the same length.
	 *  Lengths are not checked.
	 */
	S dotProduct(Array1D<S>& X) const
    {
		Real* a=dptr;
		Real* x=X.getData();
		Real sum=0.0;
		for(int j=0;j<n;j++) sum+=a[j]*x[j];
			
		return sum;
    }
	
	std::ostream& printSelf(std::ostream& os) const
    {
	     os << endl << "Array1D of dimension " << n << ":" << endl;
         for(int i=0;i<n-1;i++) os << dptr[i] << ", ";
         return os << dptr[n-1] << endl << endl;
    } 
			 
		   
}; // end Array1D



/**********************************************************************************
 *
 *         2 DIMENSIONAL   ARRAY
 *
 *********************************************************************************/


/** Lightweight rectangular array in dimension 2. 
 *
 * @param S type of array elements.
 */
template<typename S>
class Array2D {
	
protected:
	
	/** Index bases in dimension 1,2,3 */
	int b1,b2;
	/** Number of elements in dimension 1,2,3. */
	int n1,n2;         
				
	/** Data array. */
	S** dptr;

	
public:
	

// ACCESSORS

	/** <a href="index-base">Index base</a> in dimension j. 
	 */
    int getIndexBase(int j) const 
    {
		switch(j){
			
			case 1 : return b1;
			case 2 : return b2;
		}
	} // end getIndexBase
	
    /** Number \f$n_j\f$ of array elements in dimension j. */
	int getSize(int j) const 
    {
		switch(j){
			
			case 1 : return n1;
			case 2 : return n2;
		}
	} // end getIndexBase

	
// CONSTRUCTOR
	
	/** @param n_1 number of array elements in dimension 1.
	 *  @param b_1 <a href="index-base">index base</a> in dimension 1.
	 */
	Array2D(int n_1, int n_2, int b_1=0, int b_2=0) :
	b1(b_1), b2(b_2),
	n1(n_1), n2(n_2) 
	{  
         dptr=new S*[n_1];
		 for(int i=0;i<n1;i++){
			 
			 dptr[i]=new S[n2];
			 for(int j=0;j<n2;j++)dptr[i][j]=0;
		 }
	} // end constructor

	
	~Array2D()
    {
		 for(int i=0;i<n1;i++) delete[] dptr[i];
		 delete[] dptr;
	}
	
   /** Subscripting with indices based on the index bases.*/
   const S& operator()(int i, int j) const 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,j,b1,b2,n1,n2,"Array2D");
	   #endif	
	   return dptr[i-b1][j-b2]; 
   }
   
   /** Subscripting with indices based on the index bases.*/
   S& operator()(int i, int j) 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,j,b1,b2,n1,n2,"Array2D");
	   #endif	
	   return dptr[i-b1][j-b2]; 
   }
			 
   /** Pointer to row i. This allows subscripting as X[i][j]
    *  with absolute (zero based) indices i,j.
    */
   S* operator[](int i)  
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,b1,n1,"Array2D, row i");
	   #endif	
	   return dptr[i-b1];
   }
   
   std::ostream& printSelf(std::ostream& os) const
   {
	    os << endl << "Rectangular " << n1 << " by " << n2 << " array:"
	       << endl << endl;
	    for(int i=0;i<n1;i++){
		
	        for(int j=0;j<n2-1;j++) os << dptr[i][j] << ", ";
		    os << dptr[i][n2-1] << endl;
	    }
        return os << endl << endl;
	}
				

}; // end Array2D




	     


/**********************************************************************************
 *
 *         3 DIMENSIONAL   ARRAY
 *
 *********************************************************************************/


/** Lightweight rectangular array in dimension 3. 
 *
 * @param S type of array elements.
 */
template<typename S>
class Array3D {
	
protected:
	
	/** Index bases in dimension 1,2,3 */
	int b1,b2,b3;
	/** Number of elements in dimension 1,2,3. */
	int n1,n2,n3;         
				
	/** Data array. */
	S*** dptr;

	
public:
	

// ACCESSORS

	/** <a href="index-base">Index base</a> in dimension dim. 
	 */
    int getIndexBase(int j) const 
    {
		switch(j){
			
			case 1 : return b1;
			case 2 : return b2;
			case 3 : return b3;
		}
	} // end getIndexBase
	
    /** Number \f$n_j\f$ of array elements in dimension j. */
	int getSize(int j) const 
    {
		switch(j){
			
			case 1 : return n1;
			case 2 : return n2;
			case 3 : return n3;
		}
	} // end getIndexBase

	
// CONSTRUCTOR
	
	/** @param n_1 number of array elements in dimension 1.
	 *  @param b_1 <a href="index-base">index base</a> in dimension 1.
	 */
	Array3D(int n_1, int n_2, int n_3, int b_1=0, int b_2=0, int b_3=0) : 
	b1(b_1), b2(b_2), b3(b_3),
	n1(n_1), n2(n_2), n3(n_3)
	{  
         dptr=new S**[n_1];
		 for(int i=0;i<n1;i++){
			 
			 dptr[i]=new S*[n2];
			 for(int j=0;j<n2;j++){
				 
				 dptr[i][j]=new S[n3];
				 for(int k=0;k<n3;k++) dptr[i][j][k]=0;
			}
		}
	} // end constructor
	
	~Array3D()
    {
		 for(int i=0;i<n1;i++){
			 
			 for(int j=0;j<n2;j++) delete[] dptr[i][j];
			 delete[] dptr[i];
		}
		delete[] dptr;
	}
	
   /** Subscripting.
    */
   const S& operator()(int i, int j, int k) const 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,j,k,b1,b2,b3,n1,n2,n3,"Array3D");
	   #endif	
	   return dptr[i-b1][j-b2][k-b3]; 
   }
   
   /** Subscripting.
    */
   S& operator()(int i, int j, int k) 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      SubscriptCheck::
	      checkSubscript(i,j,k,b1,b2,b3,n1,n2,n3,"Array3D");
	   #endif	
	   return dptr[i-b1][j-b2][k-b3]; 
   }
			 
		   
				

}; // end Array3D




typedef Array1D<Real> RealArray1D;
typedef Array2D<Real> RealArray2D;
typedef Array3D<Real> RealArray3D;

typedef Array1D<int> IntArray1D;
typedef Array2D<int> IntArray2D;
typedef Array3D<int> IntArray3D;

typedef Array1D<long> LongArray1D;
typedef Array2D<long> LongArray2D;

typedef Array1D<unsigned long> UnsignedLongArray1D;
typedef Array2D<unsigned long> UnsignedLongArray2D;






/**********************************************************************************
 *
 *         LIBOR  ARRAY
 *
 *********************************************************************************/


/** <p>Array in dimension 2 to contain Libors under the assumption that nSteps time 
 *  steps are taken in each accrual interval. For each accrual interval \f$[T_t,T_{t+1}]\f$
 *  the array has block of nSteps rows of length n-t-1 (the number of Libors still alive
 *  in this time interval), one row for each time step. Each such row has column index
 *  base t+1.
 *
 * <p>In other words: during time step s we are in row r[s][] with t=s/nStep and elements
 * in this row are subscripted as r(s,j), j=t+1,...,n-1 corresponding to the indices of the
 * surviving Libors.
 *
 * @param S type of array elements.
 */
template<typename S>
class LiborArray2D {
	
protected:
	
	/** Dimension of the Libor process. */
	int n;  
	
	/** Number of time steps of the underlying Libor process in 
	 *  each Libor accrual interval.
	 */
	int nSteps;   
				
	/** Data array. */
	S** dptr;

	
public:
	

	
// CONSTRUCTOR
	
	/** @param dim dimension of underlying Libor process (number n of accrual intervals).
	 *  @param steps number of time step in each accrual interval.
	 */
	LiborArray2D(int dim, int steps) : 
	n(dim), nSteps(steps)
	{  
         dptr=new S*[n*nSteps];
		 for(int t=0;t<n-1;t++)
		 for(int u=0;u<nSteps;u++){
			 
		    int s=t*nSteps+u;             // number of time step	 
		    dptr[s]=new S[n-t-1];
			for(int j=0;j<n-t-1;j++)dptr[s][j]=0;
		 }
	} // end constructor

	
	~LiborArray2D()
    {
		 for(int t=0;t<n-1;t++)
		 for(int u=0;u<nSteps;u++) delete[] dptr[t*nSteps+u];
		 delete[] dptr;
	}
	
   /** Subscripting.
    */
   const S& operator()(int s, int j) const 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      checkSubscripts(s,j);
	   #endif	
	   return dptr[s][j-1-s/nSteps]; 
   }
   
   /** Subscripting.
    */
   S& operator()(int s, int j) 
   { 
	   #ifdef SUBSCRIPT_CHECK
	      checkSubscripts(s,j);
	   #endif	
	   return dptr[s][j-1-s/nSteps]; 
   }
   
private:
   
   // time step s lands in the accrual interval (T_{t-1},T_t].
   int get_t(int s)
   {
	   if(s%nSteps==0) return s/nSteps;
	   return s/nSteps+1;
   }
   
   void checkSubscripts(int s, int j)
   {
	   if((s<0)||(s>n*nSteps-1)){
		   
		   cout << "\n\nLiborArray2D: time step s = " << s 
		        << " not in [0," << n*nSteps-1 << "]"
		        << "\nTerminating.";
		   exit(0);
	   }
	   
	   int t=get_t(s);
	   if((j<t)||(j>n-1)){
		   
		   cout << "\n\nLiborArray2D: time step s = " << s << ", t = " << t
		   
		        << "\nLibor index j = " << j << " not in ["<<t<<","<<n-1<<"]"
		        << "\nTerminating.";
		   exit(0);
	   }
   } // end checkSubscript
		   
				

}; // end LiborArray2D




	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
