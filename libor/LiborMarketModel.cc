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
#include <string>

using namespace Martingale;

/*******************************************************************************
 *
 *                  LIBOR MARKET MODEL
 *
 ******************************************************************************/

   
// CONSTRUCTION


    LiborMarketModel::
	LiborMarketModel(LiborFactorLoading* fl) :
	n(fl->getDimension()),
	delta(fl->getDeltas()),
	T(fl->getTenorStructure()),
	l(fl->getInitialLibors()),
	x(fl->getInitialXLibors()),
	factorLoading(fl)
	{   } 
	
       
 
	
// FORWARD TRANSPORTING FACTORS

	 // H_0(0)
     Real LiborMarketModel::H0() const
     {
         Real f=1.0;
	     for(int k=0;k<n;k++)f*=(1+x[k]);
	     return f;
     }
	 
	 // H_i(0)
     Real LiborMarketModel::H_i0(int i) const
     {
         Real f=1.0;
	     if(i==n) return f;
		
	     for(int k=i;k<n;k++)f*=(1+x[k]);
	     return f;
     }
     
     
	 // H_i(T_t)=B_i(T_t)/B_n(T_t)
     Real LiborMarketModel::H_it(int i, int t) 
     { 
		 if(i==0)return 1.0;

		 const vector<Real>& X_t=XLvect(t,i);         // X_t[j]=X_j(T_t), j>=i
		 Real f=1.0;
         for(int j=i;j<n;j++)f*=(1.0+X_t[j]);
         return f;
     } 
	 


	 // H_t(T_t)=1/B_n(T_t)
     Real LiborMarketModel::H_ii(int t) { return H_it(t,t); }

 
	 
              
// ZERO COUPON BONDS 


     // B_i(0)
     Real LiborMarketModel::B0(int i) const
     { 
         Real f=1.0;
         // accumulate 1 from time t=0 to time t=T_i                    
         for(int j=0;j<i;j++)f*=(1.0+x[j]); 
         return 1.0/f;                                 // B_i(0)=1/f  
     }   
 

	// B_i(T_t)
    Real LiborMarketModel::B(int i, int t) 
    {     
        const vector<Real>& X_t=XLvect(t,i);         // X_t[k]=X_k(T_t), j>=i
		Real f=1.0; 
        // accumulate 1 forward from time t to time i
        for(int k=t;k<i;k++)f*=(1.0+X_t[k]); 
        return 1.0/f;                                 // B_i(T_t)=1/f  
    }
	

     

	 
// FORWARD SWAP RATES                      
             

     // S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)
     Real LiborMarketModel::swRate(int p, int q, int t) 
     { 
        Real num,        // numerator in the swap rate definition
             denom;      // denominator
        
        num=B(p,t)-B(q,t);
        denom=0;
        for(int j=p;j<q;j++)denom+=delta[j]*B(j+1,t);
       
        return num/denom;
     } //end Swap_Rate 
     

     
     // S_{pq}(T_t)=(B_p(T_t)-B_q(T_t))/B_{p,q}(T_t)
     Real LiborMarketModel::swapRate(int p, int q, int t) 
     { 
        const vector<Real>& X_t=XLvect(t,p);         // X_t[k]=X_k(T_t), k>=p
		Real f=1.0+X_t[q-1], 
		     S=delta[q-1];
        for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+X_t[k]); }
 
        return (f-1.0)/S;
     } //end Swap_Rate


     // S_{pq}(0)
     Real LiborMarketModel::swapRate(int p, int q) const
     { 
        Real f=1.0+x[q-1], S=delta[q-1];
        for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+x[k]); }
 
        return (f-1.0)/S;
     } //end Swap_Rate

 
	 
	 
// ANNUITY NUMERAIRE
                  
     
     // B_{pq}(T_t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(T_t)
     Real LiborMarketModel::b_pq(int p, int q, int t) 
     { 
         Real S=0.0;
         for(int k=p;k<q;k++)S+=delta[k]*B(k+1,t);
         return S;
     } //end B_pq
     
     
     // B_{pq}(T_t)=\sum\nolimits_{k=p}^{q-1}\delta_kB_{k+1}(T_t)    
	 Real LiborMarketModel::B_pq(int p, int q, int t) 
     { 
         const vector<Real>& X_t=XLvect(t,p);         // X_t[k]=X_k(T_i), k>=p
		 Real S=0.0, F=B(q,t);
         for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+X_t[k]); }

         return S;
     } //end B_pq

     
     // B_{pq}(0)
     Real LiborMarketModel::B_pq(int p, int q) const
     { 
         Real S=0.0, F=B0(q);
         for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+x[k]); }
         return S;
     } //end B_pq
	 
	 
	 
     // H_{p,q}(T_t)=B_{pq}(T_t)/B_n(T_t)  
	 Real LiborMarketModel::H_pq(int p, int q, int t) 
     { 
         return B_pq(p,q,t)*H_ii(t);
     } //end H_pq

     
     // H_pq(0)
     Real LiborMarketModel::H_pq(int p, int q) const
     { 
         return B_pq(p,q)*H0();
     } //end B_pq

             




/*******************************************************************************
 *
 *                          BONDS
 *
 ******************************************************************************/



// CONSTRUCTOR


    Bond::
    Bond(int k, int m, const RealArray1D& d, LiborMarketModel* lmm) :
    LMM(lmm),
	n(lmm->getDimension()), p(k), q(m), 
	c(d),
	b(q-p,p)
    {   
		c.setIndexBase(p);
		b[p]=c[p];
		for(int j=p;j<q-1;j++) b[j+1]=b[j]+c[j+1];
	}
	
	
    // zero coupon bond maturing at T_i
	Bond::
    Bond(int i, LiborMarketModel* lmm) :
    LMM(lmm),
	n(lmm->getDimension()), p(i), q(i+1), 
	c(1,p),
	b(1,p)
    {   
		c[p]=b[p]=1.0;
	}


    
// PPRICE
	
     //Forward price F(T_t)=\sum_{j=p}^{q-1}c_jH_j(T_t)
     Real Bond::
	 forwardPrice(int t) const
     {
         Real sum=0;
	   	 for(int j=p;j<q;j++){ sum+=c[j]*(LMM->H_it(j,t)); }
         return sum;	
     } //end analyticForwardPrice   

	 
     // Forward price F(0) at time 0
     Real Bond::
	 forwardPrice() const { return forwardPrice(0); }

	 
	 // Cash price B(T_t)=\sum_{j=p}^{q-1}c_jB_j(T_t)
     Real Bond::
	 cashPrice(int t) const { return forwardPrice(t)/(LMM->H_ii(t)); }

	 
     // Cash price B(0) at time 0
     Real Bond::
	 cashPrice() const { return forwardPrice(0)/(LMM->H0()); }
	 
	 
// PRINTING 

    std::ostream& Bond::printSelf(std::ostream& os) const
    {
	    return 
		os << "Bond along [T_"<<p<<",T_"<<q<<"], n="<<n<< endl
	    	<< "Coupons: " << c;
    }
	



