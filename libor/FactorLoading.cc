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

#include <string>
#include <iostream>
#include "Matrices.h"
using namespace Martingale;


/*******************************************************************************
 *
 *                     FactorLoading
 * 
 ******************************************************************************/

 

// LOG-COVARIATION MATRICES 

const UTRRealMatrix&
FactorLoading::
covariationMatrix(int p,int q, Real t, Real T) const
{
    int size=q-p;       // matrix size
	UTRRealMatrix& cvm=*(new UTRRealMatrix(size,p));

    for(int i=p;i<q;i++)
    for(int j=i;j<q;j++)
    cvm(i,j)=integral_sgi_sgj_rhoij(i,j,t,T);

    return cvm;
}
   
   
const UTRRealMatrix&
FactorLoading::
covariationMatrix(Real t, Real T) const
{
	 UTRRealMatrix& cvm=*(new UTRRealMatrix(n));

     for(int i=0;i<n;i++)
     for(int j=i;j<n;j++)
     cvm(i,j)=integral_sgi_sgj_rhoij(i,j,t,T);

     return cvm;
}
   

const UTRRealMatrix& 
FactorLoading::
nu(Real t, int p) const
{
	UTRRealMatrix c(n-p,p);
	for(int j=p;j<n;j++) 
	for(int k=j;k<n;k++) c(j,k)=sigma(j,t)*sigma(k,t)*rho(j,k);
		   
	return c.utrRoot();
}


/*******************************************************************************
 *
 *            Constant volatility-correlation factor loading
 *
 ******************************************************************************/


// CORRELATIONS, VOLATILITIES, LOG-COVARIATION INTEGRALS

Real 
ConstantFactorLoading::
rho(int i, int j) const 
{ 
	if(i<=j) return corr(i,j); return corr(j,i); 
} 
	 
   
Real 
ConstantFactorLoading::	
sigma(int i, Real t) const { return sg[i]; }
   

Real 
ConstantFactorLoading::
integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const
{ 
	return sg[i]*sg[j]*(T-t); 
}
   

const UTRRealMatrix 
ConstantFactorLoading::
correlationMatrixRoot()
{ 
	return corr.utrRoot(); 
}
   
   
const UTRRealMatrix 
ConstantFactorLoading::
correlationMatrixRankReducedRoot(int r)
{ 
	return corr.rankReducedRoot(r); 
}
   
   
   
std::ostream& 
ConstantFactorLoading::
printSelf(std::ostream& os) const
{
	return
	os << endl << endl
	   << "Factor loading with constant volatilities:" << endl << vols << endl
       << "and constant correlations:" << endl << corr;
}



