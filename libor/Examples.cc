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


#include "Examples.h"
#include "TypedefsMacros.h"
#include "Matrix.h"
#include "Utils.h"
#include "Random.h"
#include "RandomObject.h"
#include "LiborMarketModel.h"
#include "LiborFactorLoading.h"
#include "PredictorCorrectorLMM.h"
#include "FastPredictorCorrectorLMM.h"
#include "StochasticProcesses.h"


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Examples)


void 
timeSTN()
{
    Timer watch;
	watch.start();
    for(int i=0;i<100000000;i++) Random::sTN(); 
    watch.stop();
	watch.report("100,000,000 STN()s:");
} 


void 
timeMatrixMultiply(int dim, int N)
{
    RealMatrix A(dim); 
	RealMatrix B(dim);
	Timer watch;
	
	std::cout << "\n\nTiming square matrix times square matrix\n"
	          << "in dimension " << dim << ", for " << N << " products."
	          << endl << endl 
	          << "Optimal memory access pattern: " << endl;
    watch.start();
	// loop over the multiplications
	for(int n=0;n<N;n++){
		
		for(int i=0;i<dim;i++){
		
			for(int j=0;j<dim;j++) A(i,j)=n+i+j; 
			for(int j=i;j<dim;j++) B(i,j)=n-i+j; 
		}
			
		A^=B;
	}
    watch.stop();
    watch.report("time: ");
	
    std::cout << "\n\nStandard memory access pattern: " << endl;
    watch.start();
	// loop over the multiplications
	for(int n=0;n<N;n++){
		
		for(int i=0;i<dim;i++){
		
			for(int j=0;j<dim;j++) A(i,j)=n+i+j; 
			for(int j=i;j<dim;j++) B(i,j)=n-i+j; 
		}
			
		A*=B;
	}
    watch.stop();
    watch.report("time: ");
} 


void 
matrixExponentials()
{
    std::cout << endl << "Enter matrix dimension: ";
	int dim; cin >> dim;
    std::cout << endl << "Enter step size dt (roughly 0.1): ";
	Real dt; cin >> dt;
	// make string conversion available
	Int_t Dim(dim);  Real_t Dt(dt);
	
	Timer watch;
	watch.start();
 
	// We use the product ^=Q ie. *=Q' where Q=I+Adt which is faster than
	// the product *=Q. Consequently we have to work with the lower triangular 
	// transpose of A instead of a itself.
	LTRRealMatrix A(dim), Q(dim), R(dim);       
	for(int i=0;i<dim;i++)
	for(int j=0;j<=i;j++){ 
		
		A(i,j)=Random::U01();
		Q(i,j)=dt*A(i,j); R(i,j)=-Q(i,j);
		if(i==j){ Q(i,j)+=1; R(i,j)+=1; }
	}
	
	int n = (int)(1.0/dt);
	UTRRealMatrix H(dim), K(dim);
	
	// intitialization of H(t), K(t), t=0
    for(int i=0;i<dim;i++){ H(i,i)=K(i,i)=1; }
		
	// forward computation of H(t), K(t) to t=1	
	for(int k=0;k<n;k++){ H^=Q; K^=R; }
				
	std::cout << "exp(A) is the following matrix:" << H;
	//std::cout << "exp(-A) is the following matrix:" << H;
	
	H*=K;
	std::cout << "exp(A)exp(-A) is the following matrix:" << H;
		
	watch.stop();
    watch.report("Matrix exponentials, dimension "+Dim.toString()+", stepsize "+Dt.toString());


} // end matrixExponentials



             
void 
liborPathTiming(int n, int N)
{
    Timer watch;
	for(int lmmType=0;lmmType<4;lmmType++){
			
        int volType=VolSurface::JR,
		    corrType=Correlations::CS;
		LiborMarketModel* lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
        watch.start();
		for(int i=0;i<N;i++) lmm->newPath();
		watch.stop();
		cout << '\n' << *(lmm->getType());
		watch.report(" ");
	}		
} // end liborPathTiming
	
	
void 
brownianMotionInBall(int dim, int T, Real dt, int N)
{
	Region< RealVector >* D=new Ball(dim,1.0);
	VectorProcess* X=new VectorBrownianMotion(dim,T,dt);
	StoppingTime* tau=new FirstExitTime<RealVector,Real>(X,D);
		
	Real sum=0;
	for(int k=0;k<N;k++){
			
		Real t=X->newPathSegment(tau);
		sum+=t;
	}
	cout << endl << endl
	     << "Mean time to hit boundary: " << dt*sum/N;
}
	



MTGL_END_NAMESPACE(Examples)
MTGL_END_NAMESPACE(Martingale)
