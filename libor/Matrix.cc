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

#include "Matrix.h"
#include "tnt/tnt_array1d.h"
#include "tnt/tnt_array2d.h"
#include "jama/jama_eig.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "Array.h"


MTGL_BEGIN_NAMESPACE(Martingale)



/***************************************************************************************

                      Global functions

***************************************************************************************/



JAMA::Eigenvalue<Real>* 
eigenDecomposition(const UTRRealMatrix& C)
{
	int dim=C.rows();

	TNT::Array2D<Real> cv(dim,dim);
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) cv[i][j]=cv[j][i]=C.entry(i,j);
	return new JAMA::Eigenvalue<Real>(cv);
}


RealMatrix& 
rank_Reduced_Root(const UTRRealMatrix& C, int r)
{
	int dim=C.rows(),
	    a=C.getRowIndexBase();
    JAMA::Eigenvalue<Real>* eigval=eigenDecomposition(C);
	// eigenvalues sorted ascending
	TNT::Array1D<Real> l(dim);
	eigval->getRealEigenvalues(l);
	// corresponding eigenvectors
	TNT::Array2D<Real> U(dim,dim);
	eigval->getV(U);
	
	RealMatrix* root=new RealMatrix(dim,r,a,0);
	for(int i=0;i<dim;i++)
	for(int j=dim-1;j>=dim-r;j--) root->entry(i,dim-j-1)=sqrt(l[j])*U[i][j];
			
	return *root;
} // end rankReducedRoot
		
	

void 
factorAnalysis(const UTRRealMatrix& C, int r)
{
	int dim=C.rows();
    JAMA::Eigenvalue<Real>* eigval=eigenDecomposition(C);
	// eigenvalues sorted ascending
	TNT::Array1D<Real> l(dim);
	eigval->getRealEigenvalues(l);
		
	Real traceNorm=l[0];
	for(int i=1;i<dim;i++) traceNorm+=l[i];
			
	Real sum=0.0; 
	cout << endl << endl
		 << "Variability captured by the i largest eigenvalues: ";
	for(int i=dim-1;i>=dim-r;i--){
			
		sum+=l[i];
		cout << "\ni=" << dim-i << ": " << 100.0*sum/traceNorm << "%";
	}
} // end factorAnalysis
	
	

void 
factorizationTest(const UTRRealMatrix& C, int r, string message="")
{		
	int dim=C.rows();
	RealMatrix D(dim);	
	for(int i=0;i<dim;i++)
	for(int j=i;j<dim;j++) D(i,j)=D(j,i)=C.entry(i,j);
	Real D_norm=D.norm();
	
	RealMatrix& R=rank_Reduced_Root(C,r);
	
	R^=R;      // R*=R'
	R*=(-1.0);
	D+=R;      // D=D-RR'
	Real err=D.norm();
			
	cout << message << "\nerror:  "<< 100*err/D_norm << "%";
		
} // end factorizationTest
	



/***************************************************************************************
 *
 *                             UPPER TRIANGULAR MATRICES
 *
***************************************************************************************/

RealMatrix&
UTRRealMatrix::
symmetricCompletion() const
{
     RealMatrix* D = new RealMatrix(rows(),cols(),a,b);
	 for(int i=0;i<rows();i++)
	 for(int j=i;j<cols();j++) D->entry(j,i)=D->entry(i,j)=entry(i,j);
		 
	 return *D;
 }
	
	
		
UTRRealMatrix& 
UTRRealMatrix::
inverse() const
{
	int dim=rows();
	for(int i=0;i<dim;i++){
	  if(dptr[i][0]==0.0){ 
		 cout << "UTRMatrix::inverse(): matrix is singular, quitting."; 
		 exit(0); 
	  }
	} // end for i
    UTRRealMatrix* inverse = new UTRRealMatrix(dim,a);
			
    Real** A=dptr; 
	Real** I=inverse->getData();	
	// normalize diagonal
	for(int i=0;i<dim;i++) I[i][0]=1.0/A[i][0];       // A[i][0]=A_ii
	

	for(int j=dim-1;j>=0;j--)
	for(int i=j-1;i>=0;i--) 
	// subtract row_j*A_ij/A_ii from row_i
	for(int k=j;k<dim;k++) I[i][k-i]-=I[j][k-j]*A[i][j-i]/A[i][0];
		
	return *inverse;
}
	


LTRRealMatrix& 
UTRRealMatrix::
ltrRoot() const
{
    LTRRealMatrix* ltr=new LTRRealMatrix(rows(),a);
	Real** L=ltr->getData();
  
    // computation of Lij reverse induction starting from i,j=0 ascending
	// get U_ij from relation row_i(U) dot row_j(U) = A_ij, j<=i, at that point all
	// U_rs with r<=i,s<=j (and not both equal) are already computed.
    for(int i=0;i<rows();i++)
    for(int j=0;j<=i;j++)                         // j\leq i
    { 
        // r_i(L).r_j(L)-L_ijL_jj
		Real sum=0.0;  
        for(int k=0;k<j;k++) sum += L[i][k]*L[j][k];   
      
		Real Aij=dptr[j][i-j], R=Aij-sum;                // R-L_ijL_jj=A_ij-r_i(L).r_j(L)=0
        if(j<i) L[i][j]=R/L[j][j];                    // L[j][j] already computed, compute L[i][j]
        else if(R>0) L[j][j]=sqrt(R);
        else{ 
			
			 cerr << "\n\nUTRMatrix::ltrRoot(): matrix not positive definite:" 
			     << "\nDimension "<< dim <<", A["<<j<<"]["<<j<<"]-S="<< R << endl  
			     << "\n\nMatrix:\n\n" << *this << endl << "Terminating.";
               exit(1); 
        } // end else
    } //end for i
	
	return *ltr;

} //end ltrRoot



UTRRealMatrix& 
UTRRealMatrix::
utrRoot() const
{
    int dim=rows();
	UTRRealMatrix* utr=new UTRRealMatrix(dim,a);
	Real** U=utr->getData();
	
	int n=dim-1;
  
	// computation of Lij reverse induction starting from i,j=n descending
	// get U_ij from relation row_i(U) dot row_j(U) = A_ij, j>=i, at that point all
	// U_rs with r>=i,s>=j (and not both equal) are already computed.
    for(int i=n;i>=0;i--)
    for(int j=n;j>=i;j--) {                      // j>=i
    
        // r_i(U).r_j(U)-U_ijU_jj
		Real sum=0.0;  
        for(int k=j+1;k<dim;k++) sum += U[i][k-i]*U[j][k-j];   

   		Real Aij=dptr[i][j-i], R=Aij-sum;        // R-U_ijU_jj=A_ij-r_i(U).r_j(U)=0    
        if(j>i) U[i][j-i]=R/U[j][0];             // U[j][0]=U[j][j-j]=U_jj              
        else if(R>0) U[j][0]=sqrt(R);
        else{ 
			
			cerr << "\n\nUTRMatrix::utrRoot(): matrix not positive definite:" 
			     << "\nDimension "<< dim <<", A["<<j<<"]["<<j<<"]-S="<< R << endl  
			     << "\n\nMatrix:\n\n" << *this << endl << "Terminating.";
               exit(1); 
        } // end else
    } //end for i
	
	return *utr;

} //end utrRoot



RealMatrix& 
UTRRealMatrix::
rankReducedRoot(int r) const { return rank_Reduced_Root(*this,r); }



UTRRealMatrix& 
UTRRealMatrix::
exp() const 
{
	UTRRealMatrix* E = new UTRRealMatrix(rows(),a);
	Matrix< Real,UpperTriangular<Real> >& 
	exp_ =  Matrix< Real,UpperTriangular<Real> >::exp();
	for(int i=0;i<rows();i++)
	for(int j=i;j<rows();j++) E->entry(i,j)=exp_.entry(i,j);
	delete &exp_;
	return *E;
}



void 
UTRRealMatrix::
analyseFactors(int r) const { factorAnalysis(*this,r); }



void 
UTRRealMatrix::
testFactorization(int r, string message="") const { factorizationTest(*this,r,message); }



RealMatrix& 
UTRRealMatrix::
matrixFunction(Real (*f)(Real))
{
	int dim=rows();
    JAMA::Eigenvalue<Real>* eigval=eigenDecomposition(*this);
	// eigenvalues sorted ascending
	TNT::Array1D<Real> l(dim);
	eigval->getRealEigenvalues(l);
	// corresponding eigenvectors
	TNT::Array2D<Real> V(dim,dim);
	eigval->getV(V);
	
	RealMatrix& f_A = *(new RealMatrix(dim,dim,a,b));
	RealMatrix& U  = *(new RealMatrix(dim,dim,a,b));     // matrix V
	for(int i=0;i<dim;i++)
	for(int j=0;j<dim;j++){ 
		
		f_A.entry(i,j)=f(l[j])*V[i][j];     // U*f(D), with D=diag(l_j), f(D)=diag(f(l_j))
		U.entry(i,j)=V[i][j]; 
	}
	
	f_A^=U;
	delete &U;		
	return f_A;
} 




MTGL_END_NAMESPACE(Martingale)

