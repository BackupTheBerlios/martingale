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

#include "StochasticGenerator.h"
#include "QuasiMonteCarlo.h"
#include "Matrices.h"
#include <iostream>
using namespace Martingale;


/*********************************************************************************
 *
 *         Stochastic Generator (Driver)       
 *
 *********************************************************************************/

std::ostream& 
StochasticGenerator::
printSelf(std::ostream& os) const 
{
	   return os << "\nUnspecified stochastic generator.";
}




// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, const StochasticGenerator& sg){ return sg.printSelf(os); }



/*********************************************************************************
 *
 *         Libor Drivers   
 *
 *********************************************************************************/

   
void 
MonteCarloLiborDriver::
newWienerIncrements(int t, int T, UTRRealMatrix& Z)
{
	 for(int s=t;s<T;s++)
	 for(int k=s+1;k<n;k++)Z(s,k)=Random::sTN();
}


std::ostream& 
MonteCarloLiborDriver::
printSelf(std::ostream& os) const { return os << "Mersenne Twister."; }




void 
SobolLiborDriver::
newWienerIncrements(int t, int T, UTRRealMatrix& Z)
{		
	// number of normal deviates needed per path
	int d = (T-t)*(2*n-(1+t+T))/2;
	// check if Sobol generator is intialized in the right dimension
	if(lds==0) lds=new Sobol(d);
	else if ((lds->getDimension())!=d){ delete lds; lds=new Sobol(d); }

	const RealArray1D& x=lds->nextQuasiNormalVector();
	int coordinate=0;
    for(int s=t;s<T;s++)
    for(int k=s+1;k<n;k++){ Z(s,k)=x[coordinate]; coordinate++; }
     
} // end newWienerIncrements
	
	
void 
SobolLiborDriver::
restart(){ if(lds) lds->restart(); }
	

std::ostream& 
SobolLiborDriver::
printSelf(std::ostream& os) const { return os << "Sobol sequence."; }



/*********************************************************************************
 *
 *         VectorProcess Drivers       
 *
 *********************************************************************************/

void 
MonteCarloVectorDriver::
newWienerIncrements(int t, int s, Real** Z)
{
	 for(int u=t;u<s;u++)
     for(int k=0;k<n;k++) Z[u][k]=Random::sTN();
}
	
void 
MonteCarloVectorDriver::
newWienerIncrements(int t, int s, RealMatrix& Z)
{   newWienerIncrements(t,s,Z.getData());  }
	
	
std::ostream& 
MonteCarloVectorDriver::
printSelf(std::ostream& os) const { return os << "Mersenne Twister."; }
	


void 
SobolVectorDriver::
newWienerIncrements(int t, int s, Real** Z)
{		
    // Use s Sobol vector of full dimension n*T otherwise the effective
	// dimension of the simulation is underestimated if a path is computed 
	// as a series of path segments.			
	if(lds==0) lds=new Sobol(n*T);
    const RealArray1D& x=lds->nextQuasiNormalVector();
    int coordinate=0;
	for(int u=t;u<s;u++)
    for(int k=0;k<n;k++){ Z[u][k]=x[coordinate]; coordinate++; }
} 
	
    
void 
SobolVectorDriver::
newWienerIncrements(int t, int s, RealMatrix& Z)
{   
	newWienerIncrements(t,s,Z.getData());  
}
	

void 
SobolVectorDriver::
restart(){ if(lds) lds->restart(); }
	

std::ostream& 
SobolVectorDriver::
printSelf(std::ostream& os) const { return os << "Sobol sequence."; }


/*********************************************************************************
 *
 *         ScalarProcess Drivers    
 *
 *********************************************************************************/


void 
MonteCarloScalarDriver::
newWienerIncrements(int t, int s, Real* Z)
{
	 for(int u=t;u<s;u++) Z[u]=Random::sTN();
}
	

std::ostream& 
MonteCarloScalarDriver::
printSelf(std::ostream& os) const { return os << "Mersenne Twister."; }



void 
SobolScalarDriver::
newWienerIncrements(int t, int s, Real* Z)
{		
	// Use s Sobol vector of full dimension T otherwise the effective
    // dimension of the simulation is underestimated if a path is computed 
    // as a series of path segments.
	if(lds==0) lds=new Sobol(T);
	const RealArray1D& x=lds->nextQuasiNormalVector();
    for(int u=t;u<s;u++) Z[u]=x[u]; 
} 
	

void 
SobolScalarDriver::
restart(){ if(lds) lds->restart(); }
	

std::ostream& 
SobolScalarDriver::
printSelf(std::ostream& os) const { return os << "Sobol sequence."; }



