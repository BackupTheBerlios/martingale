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

#ifndef martingale_tnt_h    
#define martingale_tnt_h

#include "TypedefsMacros.h"
#include "tnt/tnt_array1d.h"
#include "tnt/tnt_array2d.h"
#include "jama/jama_eig.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "Utils.h"                      



MTGL_BEGIN_NAMESPACE(Martingale)

/** <p>Static methods to diagonalize a symmetric matrix C, associated eigenanalysis
 *  and approximate factorizations C=RR' with R of low rank. Eigenanalysis of 
 *  Libor factor loadings.
 *
 * <p>Redundant, has now been integrated into Matrices.h.
 */
class TnT {
	

	
public:
	
	
	static JAMA::Eigenvalue<Real>* eigenDecomposition(UTRMatrix<Real>& C)
    {
		int dim=C.getDimension();
		Real** c=C.getData();
		TNT::Array2D<Real> cv(dim,dim);
		for(int i=0;i<dim;i++)
		for(int j=i;j<dim;j++) cv[i][j]=cv[j][i]=c[i][j-i];
		return new JAMA::Eigenvalue<Real>(cv);
	}

	/** Let D be the symmetric matrix with upper half C. This function computes 
	 *  the matrix R which best approximates D as a product RR' and has rank r. The matrix R
	 *  is computed by diagonalizing D, setting all but the r largest eigenvalues
	 *  equal to zero and taking the square root of the remaining eigenvalues.
	 */
	static Matrix<Real>& rankReducedRoot(UTRMatrix<Real>& C, int r)
    {
		int dim=C.getDimension();
        JAMA::Eigenvalue<Real>* eigval=eigenDecomposition(C);
		// eigenvalues sorted ascending
		TNT::Array1D<Real> l(dim);
		eigval->getRealEigenvalues(l);
		// corresponding eigenvectors
		TNT::Array2D<Real> U(dim,dim);
		eigval->getV(U);
	
		Matrix<Real>& root=*(new Matrix<Real>(dim,r,0,0));
		for(int i=0;i<dim;i++)
		for(int j=dim-1;j>=dim-r;j--) root(i,dim-j-1)=sqrt(l[j])*U[i][j];
			
		return root;
	} // end rankReducedRoot

	/** Let D be the symmetric matrix with upper half C. This function computes 
	 *  the a matrix R which approximates D as a the product RR' and has rank r. The matrix R
	 *  is computed by setting all but the r last columns of the upper triangular root of
	 *  D equal to zero. If D is the covariance matrix of a multinormal vector Y then this
	 *  corresponds to conditioning Y on its last two components.
	 */
	static Matrix<Real>& conditioningRoot(UTRMatrix<Real>& C, int r)
    {
		Real**	rt=C.utrRoot().getData();
		int dim=C.getDimension();
		Matrix<Real>& root=*(new Matrix<Real>(dim,r,0,0));
		for(int j=dim-r;j<dim;j++) 
		for(int i=0;i<=j;i++)	root(i,j-dim+r)=rt[i][j-i];
			
		return root;
	} // end conditioningRoot

	
	static Real error(UTRMatrix<Real>& C,Matrix<Real>& P)
    {
		int dim=C.getDimension();
		Real** c=C.getData();
		Real** p=P.getData();
		
		Real S=0, delta;
		for(int i=0;i<dim;i++)
		for(int j=i;j<dim;j++) { delta=c[i][j-i]-p[i][j]; S+=delta*delta; }
			
		return sqrt(S);
	} // end error
		
	/** C is interpreted as the upper half of a multinormal covariance matrix.
	 *  This method prints how much variability is captured by the 5 largest
	 *  eigenvalues of C.
	 */
	static void factorAnalysis(UTRMatrix<Real>& C)
    {
		int dim=C.getDimension();
        JAMA::Eigenvalue<Real>* eigval=eigenDecomposition(C);
		// eigenvalues sorted ascending
		TNT::Array1D<Real> l(dim);
		eigval->getRealEigenvalues(l);
		
		Real traceNorm=l[0];
		for(int i=1;i<dim;i++) traceNorm+=l[i];
			
		Real sum=0.0; 
		cout << endl << endl
		     << "Variability captured by the i largest eigenvalues: ";
		for(int i=dim-1;i>=dim-5;i--){
			
			sum+=l[i];
			cout << "\ni=" << dim-i << ": " << 100.0*sum/traceNorm << "%";
		}
	} // end factorAnalysis
	
	
    /** Examines how much variability is 
	 *  captured by the 5 largest eigenvalues and 
	 *  tests the two approximate rank r factorizations of C.
	 */
	static void factorizationTest(UTRMatrix<Real>& C, int r)
    {		
		factorAnalysis(C);
		
		Matrix<Real>& R0=rankReducedRoot(C,r);
		Matrix<Real>& R1=conditioningRoot(C,r);
		
		R0^=R0; R1^=R1;  // R*=R'
			
		cout << "\n\n\n\nTrace norm of the covariation matrix C: " << C.norm()
		     << "\nError of the rank " << r << " factorization (trace norm):"
		     << "\n\nfactorization with diagonalization, error: " << error(C,R0)
		     << "\nfactorization with conditioning on last " << r << " variables, error: " << error(C,R1);
		
	} // end factorizationTest
	
	

	
}; // end TnT
	
	
	
	



MTGL_END_NAMESPACE(Martingale)

#endif
 
