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

#include "LiborMarketModel.h"
#include "DriftlessLMM.h"
#include "LowFactorDriftlessLMM.h"
#include "PredictorCorrectorLMM.h"
#include "FastPredictorCorrectorLMM.h"
#include "Array.h"
#include "QuasiMonteCarlo.h"
#include "LiborFactorLoading.h"
#include <iostream>
#include <string>

using std::ostream;
using std::string;


MTGL_BEGIN_NAMESPACE(Martingale)



/*******************************************************************************
 *
 *                     LiborMarketModelType
 * 
 ******************************************************************************/
 

// global insertion
ostream& 
LiborMarketModelType::
printSelf(ostream& os) const
{
	string simulation;
	switch(type){
		
		case LiborMarketModelType::DL   :   simulation="driftless"; break;
		case LiborMarketModelType::LFDL :   simulation="low factor driftless"; break;
		case LiborMarketModelType::PC   :   simulation="predictor-corrector"; break;
		default                         :   simulation="fast predictor-corrector";
	}
	
	return
	os << "\n\nLibor market model:"
	   << "\nSimulation: " << simulation 
	   << *flType;
}



/*******************************************************************************
 *
 *                  LIBOR MARKET MODEL
 *
 ******************************************************************************/



string 
LiborMarketModel::
lmmType(int type)
{
	switch(type){
			 
		 case PC   : return "PC";
         case FPC  : return "FPC";
		 case DL   : return "DL";
	     case LFDL : return "LFDL";
	}
	return "Unknown_LMM_type";
}

   
// CONSTRUCTORS

LiborMarketModel::
LiborMarketModel(LiborFactorLoading* fl, int lmmType) :
type(fl->getType(), lmmType),
n(fl->getDimension()),
delta(fl->getDeltas()),
T(fl->getTenorStructure()),
l(fl->getInitialLibors()),
x(fl->getInitialXLibors()),
factorLoading(fl)
{   } 


LiborMarketModel* 
LiborMarketModel::
sample(int n, int lmmType, int volType, int corrType)
{	
    LiborMarketModel* lmm=0;
	switch(lmmType){

		case PC :
        lmm=PredictorCorrectorLMM::sample(n,volType,corrType); break;
		
		case FPC :
        lmm=FastPredictorCorrectorLMM::sample(n,volType,corrType); break;
		
	    case LFDL :
        lmm=LowFactorDriftlessLMM::sample(n,volType,corrType); break;
		
		default :
        lmm=DriftlessLMM::sample(n,volType,corrType); 
	}
	
	return lmm;
}


// PATH GNERATION


void 
LiborMarketModel::	
newPath()
{
     newWienerIncrements(n-1);
     for(int t=0;t<n-1;t++)timeStep(t,1);
}


void 
LiborMarketModel::	
newPath(int t, int p)
{
     newWienerIncrements(t);
     for(int s=0;s<t;s++)timeStep(s,p);
}


	

// LIBOR VOLATILITIES
	
Real 
LiborMarketModel::	
vol(int i) const 
{ 
	return factorLoading->annualVol(i); 
}
 
	
// FORWARD TRANSPORTING FACTORS

Real
LiborMarketModel::	
L(int j, int t) const { return XL(j,t)/delta[j]; }

// H_0(0)
Real 
LiborMarketModel::
H0() const
{
     Real f=1.0;
     for(int k=0;k<n;k++)f*=(1+x[k]);
     return f;
}
	 
	
// H_i(0)
Real 
LiborMarketModel::
H_i0(int i) const
{
     Real f=1.0;
     if(i==n) return f;
		
     for(int k=i;k<n;k++)f*=(1+x[k]);
     return f;
}
     
     
// H_i(T_t)=B_i(T_t)/B_n(T_t)
Real 
LiborMarketModel::
H_it(int i, int t) 
{ 
     if(i==n)return 1.0;
	 const RealVector& X_t=XLvect(t,i);         // X_t[j]=X_j(T_t), j>=i
	 Real f=1.0;
     for(int j=i;j<n;j++)f*=(1.0+X_t[j]);
     return f;
} 
	 

// H_t(T_t)=1/B_n(T_t)
Real 
LiborMarketModel::
H_ii(int t) { return H_it(t,t); }

 
	 
// ZERO COUPON BONDS 

// B_i(0)
Real 
LiborMarketModel::
B0(int i) const
{ 
    Real f=1.0;
	// accumulate 1 from time t=0 to time t=T_i                    
    for(int j=0;j<i;j++)f*=(1.0+x[j]); 
    return 1.0/f;                                 // B_i(0)=1/f  
}   
 

// B_i(T_t)
Real 
LiborMarketModel::
B(int i, int t) 
{ 
    const RealVector& X_t=XLvect(t,t);         // X_t[k]=X_k(T_t), j>=i
	Real f=1.0; 
    // accumulate 1 forward from time t to time i
    for(int k=t;k<i;k++)f*=(1.0+X_t[k]); 

    return 1.0/f;                                 // B_i(T_t)=1/f  
}
	
	 
// FORWARD SWAP RATES                      
             
// S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)
Real 
LiborMarketModel::
swRate(int p, int q, int t) 
{ 
    Real num,        // numerator in the swap rate definition
         denom;      // denominator
        
    num=B(p,t)-B(q,t);
    denom=0;
    for(int j=p;j<q;j++)denom+=delta[j]*B(j+1,t);
       
    return num/denom;
} //end Swap_Rate 
     

// S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)
Real 
LiborMarketModel::
swapRate(int p, int q, int t) 
{ 
    const RealVector& X_t=XLvect(t,p);         // X_t[k]=X_k(T_t), k>=p
	Real f=1.0+X_t[q-1], 
	     S=delta[q-1];
    for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+X_t[k]); }
 
    return (f-1.0)/S;
} //end Swap_Rate


// S_{pq}(0)
Real 
LiborMarketModel::
swapRate(int p, int q) const
{ 
     Real f=1.0+x[q-1], S=delta[q-1];
     for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+x[k]); }
 
     return (f-1.0)/S;
} //end Swap_Rate

 

// ANNUITY NUMERAIRE
                  
// B_{pq}(T_t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(T_t)
Real 
LiborMarketModel::
b_pq(int p, int q, int t) 
{ 
     Real S=0.0;
     for(int k=p;k<q;k++)S+=delta[k]*B(k+1,t);
     return S;
} //end B_pq
     
     
// B_{pq}(T_t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(T_t)    
Real 
LiborMarketModel::
B_pq(int p, int q, int t) 
{
     const RealVector& X_t=XLvect(t,p);         // X_t[k]=X_k(T_i), k>=p
     Real S=0.0, F=B(q,t);
     for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+X_t[k]); }

     return S;
} //end B_pq

     
// B_{pq}(0)
Real 
LiborMarketModel::
B_pq(int p, int q) const
{ 
     Real S=0.0, F=B0(q);
     for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+x[k]); }
     return S;
} //end B_pq
	 
	 
	 
// H_{p,q}(T_t)=B_{pq}(T_t)/B_n(T_t)  
Real 
LiborMarketModel::
H_pq(int p, int q, int t) 
{ 
    return B_pq(p,q,t)*H_ii(t);
} 

     
// H_pq(0)
Real 
LiborMarketModel::
H_pq(int p, int q) const
{ 
     return B_pq(p,q)*H0();
} 
	 
	 
Real 
LiborMarketModel::	 
capletAggregateVolatility(int i) const
{ 
    std::cout << "LiborMarketModel#capletAggregateVolatility(): " 
              << "not implemented in this generality, aborting.";
    exit(1);
    return 0.0;
} 

             
Real 
LiborMarketModel::
swaptionAggregateVolatility(int p, int q, int t) const
{ 
    std::cout << "LiborMarketModel#swaptionAggregateVolatility(): " 
         << "not implemented in this generality, aborting.";
    exit(1);
    return 0.0;
} 


Real 
LiborMarketModel::
bondAggregateVolatility(Bond* B, int t) const
{ 
     cerr << "LiborMarketModel#bondAggregateVolatility(): " 
          << "not implemented in this generality, aborting.";
	 exit(1);
	 return 0.0;
} 





	
MTGL_END_NAMESPACE(Martingale)


