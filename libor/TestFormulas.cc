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


#include "TestFormulas.h"
#include "Random.h"
#include "QuasiMonteCarlo.h"
#include "VolatilityAndCorrelation.h"
#include "Utils.h"
#include <cmath>

using std::exp;



MTGL_BEGIN_NAMESPACE(Martingale)


void 
testExponentialIntegralFormulas()
{
	 printStars();
	 cout << "\n\nTesting some analytic formulas for expectations with Monte Carlo.";
	
	int N=10000000;
	Real analytic, montecarlo, sum;
	Real K,a,b,c,alpha,beta,N1,N2,D;
	
	// FORMULA C.12
	cout << "\n\nFormula C.12:";
	a=1.3, b=1.2, alpha=0.9, beta=1.6;
	
	// Monte  carlo expectation
	sum=0;
	for(int i=0;i<N;i++){
					   
			Real X=Random::sTN(), 
		         A=std::exp(a*X+b), 
		         B=Random::N(alpha*X+beta);
			sum+=A*B;
	}
	montecarlo=sum/N;
	
	// analytic expectation
	N1=beta+a*alpha;
	D=sqrt(1+alpha*alpha);
	analytic=std::exp(b+0.5*a*a)*Random::N(N1/D);

	cout << "\nAnalytic: " << analytic
	     << "\nMonte Carlo: " << montecarlo;	
	

	// FORMULA C.13
	cout << "\n\nFormula C.13:";
	K=0.8; a=1.3;
	
	// Monte  carlo expectation
	sum=0;
	for(int i=0;i<N;i++){
					   
			Real Y=Random::sTN(), A=std::exp(a*Y);
			sum+=(A>K)? A-K : 0.0;
	}
	montecarlo=sum/N;
	
	// analytic expectation
	N1=a*a-std::log(K);
	N2=-std::log(K);
	D=a;
	analytic=std::exp(0.5*a*a)*Random::N(N1/D)-K*Random::N(N2/D);

	cout << "\nAnalytic: " << analytic
	     << "\nMonte Carlo: " << montecarlo;	
	

	// FORMULA C.16
	cout << "\n\nFormula C.16:";
	a=0.3; b=0.4; c=0.5;
	Real y=1.2, z=2.3;
	
	// Monte  carlo expectation
	sum=0;
	for(int i=0;i<N;i++){
					   
			Real z1=Random::sTN(),
			     z2=Random::sTN(),
			Y=y+a*z1+b*z2,
			Z=z+c*z2,
			A=std::exp(Y), B=K*std::exp(Z);
			
			sum+=(A>B)? A-B : 0.0;
	}
	montecarlo=sum/N;
	
	// covariation matrix
	Real C11=a*a+b*b, C12=b*c, C22=c*c;
	cout << "\nCorrelation: " << C12/sqrt(C11*C22);
	// analytic expectation
	N1=C11-C12-std::log(K)+y-z;
	N2=C12-C22-std::log(K)+y-z;
	D=sqrt(C11+C22-2*C12);
	analytic=std::exp(y+0.5*C11)*Random::N(N1/D)-K*std::exp(z+0.5*C22)*Random::N(N2/D);

	cout << "\nAnalytic: " << analytic
	     << "\nMonte Carlo: " << montecarlo;
	
} // end testExponentialIntegralFormulas
	


void 
testVolSurfaceIntegrals(int N, Real precision)
{
	 printStars();
	 cout << "\n\nTesting the analytic VolSurface integrals "
	      << "with Monte Carlo.";
	VolSurface* vol;
	
	vol=VolSurface::sample(VolSurface::JR);
	vol->testVolSurfaceIntegrals(N,precision);
	
	vol=VolSurface::sample(VolSurface::M);
    vol->testVolSurfaceIntegrals(N,precision);
}
	
	
	     
	
MTGL_END_NAMESPACE(Martingale)
