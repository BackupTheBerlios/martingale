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


/*
 * Derivatives.h
 *
 * Created on March 21, 2003, 12:45 PM
 */
 


#include "Option.h"
#include "LmmLattice.h"
//#include "BasketLattice.h"
#include "FinMath.h"
#include "Random.h"
#include "Bond.h"
#include "DriftlessLMM.h"
#include "LowFactorDriftlessLMM.h"
#include "LiborMarketModel.h"          // Bond
#include <algorithm>                   // max
#include <stdlib.h>                    // exit



MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::endl;



/*******************************************************************************
 *
 *                     OPTION
 *
 ******************************************************************************/
 
 
const RealArray1D
Option::
nextControlledForwardPayoff()
{
    RealArray1D Z(2); Z[0]=-1.0; Z[1]=0.0;
	return Z;
}


Real
Option::
analyticForwardPrice() const
{
	cout << "\n\nGeneric Option::analyticForwardPrice(): "
	     << "not implemented in this generality."
	     << "\nReturning value -1.0.";
	return -1.0;
}


Real
Option::
monteCarloForwardPrice(int N)
{
    return Pricing::monteCarloForwardPrice(this,N);
}


Real
Option::
controlledMonteCarloForwardPrice(int N)
{
    return Pricing::controlledMonteCarloForwardPrice(this,N);
}


Real
Option::
correlationWithControlVariate(int N)
{
    return Pricing::correlationWithControlVariate(this,N);
}


Real
Option::
latticeForwardPrice()
{
	cout << "\n\nGeneric Option::latticeForwardPrice(): "
	     << "not implemented in this generality."
	     << "\nReturning value -1.0.";
	return -1.0;
}

ostream& 
Option::
printSelf(ostream& os) const
{
	return cout << "\nGeneric Option.";
}



/*******************************************************************************
 *
 *                          LIBOR DERIVATIVE
 *
 ******************************************************************************/


LiborDerivative::
LiborDerivative
(LiborMarketModel* lmm, int s, bool an, bool lt, bool mc, bool cv) : 
// T=lmm->getTenorStructure()[s] is the continuous time to expiry
Option(lmm->getTenorStructure()[s],an,lt,mc,cv),
LMM(lmm), t(s)
{  
	// lattice pricing assumes constant vols and driftless LMM
	int volType = lmm->getType()->flType->volType;
	int lmmType = lmm->getType()->type;
	if((volType!=VolSurface::CONST)||
	   (lmmType!=LiborMarketModel::DL))
	hasLattice_=false;
}


int 
LiborDerivative::
effectiveDimension() const { return LMM->effectiveDimension(0,t); }
	
	
ostream& 
LiborDerivative::
printSelf(ostream& os) const { return os << "Generic Libor Derivative"; }


LmmLattice* 
LiborDerivative::
getDefaultLattice()
{
	LiborFactorLoading* fl = LMM->getFactorLoading();
	// time steps per accrual interval: 
	int steps=1;                               
	// lattice built up to time T_t, total time steps: t*nSteps
	while(t*steps+t<LATTICE_MAX_STEPS) steps++;
	bool verbose = false;
	int r=2;   // number of factors
	return new LmmLattice(r,fl,t,steps,verbose);
}


Real 
LiborDerivative::
forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s)
{
	cout << "\n\nLiborDerivative::forwardPayoff(StandardBrownianNode*, LmmLattice*, int):"
	     << "\nNot implemented in this generality. Aborting.";
	exit(1);
	return -1.0;
}


Real 
LiborDerivative::
latticeForwardPrice(LmmLattice* theLattice)
{
	return Pricing::latticeForwardPrice(theLattice,this);
}


void 
LiborDerivative::
testPrice(int N)
{ 
	 // all prices forward prices at time T_n
     Real aPrice,               // analytic price
          mcPrice,              // Monte carlo price
          cvmcPrice,            // Monte carlo price with control variate
	      lPrice,               // Price in the default lattice
	      corr,                 // correlation with control variate
	      epsilon=0.00000001;   // replacement for zero denominator
         
	  // message
	  cout << "\n\n\n\n\n" << *this << *LMM
           << "\nEffective dimension of the simulation = " << effectiveDimension();
	  if(hasMonteCarlo_)
	  cout << "\nNumber of samples (possibly paths): " << N;
	  cout << "\n\n\n\nFORWARD PRICES: " << endl << endl;  
		
	  // PRICES: 
	
	  // analytic
	  if(hasAnalytic_){
      
		  aPrice=analyticForwardPrice();
	      cout << "Analytic: " << aPrice << endl << endl;
	  }
	
	  // lattice 
	  if(hasLattice_){
		  
	      LmmLattice* theLattice = getDefaultLattice();
	      lPrice = latticeForwardPrice(theLattice);
	      cout << "Default lattice: " << lPrice;
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,lPrice,epsilon) << "%" ;
		  cout << endl << endl;
		  
		  theLattice->rescaleVols();
		  lPrice = latticeForwardPrice(theLattice);
	      cout << "Default lattice (volatilities rescaled): " << lPrice;
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,lPrice,epsilon) << "%";
		  cout << endl << endl;
		  delete theLattice;
      }
	  
	  
	  if(hasMonteCarlo_){
		  
          mcPrice=monteCarloForwardPrice(N);
	      cout << "Monte Carlo: " << mcPrice;
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,mcPrice,epsilon) << "%"; 
		  cout << endl << endl;
	  }
	  
	  
	  if(hasControlVariate_){
		  
		  corr = correlationWithControlVariate(1000);
	      cvmcPrice=controlledMonteCarloForwardPrice(N);
          cout << "Monte Carlo with control variate: " << cvmcPrice
		       << "\nCorrelation with control variate: " << 100.0*corr << "%";
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,cvmcPrice,epsilon) << "%";
		  cout << endl << endl;
	  }
		   
	  		   
      // QMC dynamics
	  LMM->switchToQMC();
	  
	  if(hasMonteCarlo_){
		  
          mcPrice=monteCarloForwardPrice(N);
          cout << "Quasi Monte Carlo: " << mcPrice;
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,mcPrice,epsilon) << "%"; 
		  cout << endl << endl;
	  }
		   
	  LMM->restartSobolGenerator();
	  
	  if(hasControlVariate_){
		  
	      cvmcPrice=controlledMonteCarloForwardPrice(N);		   
		  cout << "Quasi Monte Carlo with control variate: " << cvmcPrice;
		  if(hasAnalytic_)
		  cout << "\nRelative error: " << relativeError(aPrice,cvmcPrice,epsilon) << "%";
		  cout << endl << endl;
	  }

} // end testPrice
	


/*******************************************************************************
 *
 *                          CAPLETS
 *
 ******************************************************************************/


	

Caplet::
Caplet(int k, Real strike, LiborMarketModel* lmm) :
 // Libors L_j, j>=k needed until time t=min(k+1,n-1)
LiborDerivative
(lmm,min(k+1,lmm->getDimension()-1),true,true,true,true),              
i(k), 
kappa(strike),
delta_i((lmm->getDeltas())[i])
{   } 


Caplet* 
Caplet::
sample(int n, int lmmType, int volType, int corrType)
{
	LiborMarketModel* lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
	int i=n/3;
	Real kappa=(lmm->getInitialLibors())[i];
	return new Caplet(i,kappa,lmm);
}
		
	 

Real 
Caplet::
nextForwardPayoff()
{                         
	// Libors X_j, j>=i until time T_t, t=min(i+1,n-1)
	LMM->newPath(t,i);    
	
	Real X_iT_i,h,f;
    X_iT_i=LMM->XL(i,i);                     // Libor X_i(T_i)
    h=max(X_iT_i-delta_i*kappa,0.0);         // payoff at time T_{i+1}
	// H_{i+1}(T_{i+1}) accrual factor T_{i+1}->T_n at time T_{i+1}
    f=LMM->H_ii(i+1);                        
    return f*h;	
}

	 

Real
Caplet::
controlVariateMean(){ return (LMM->X_i0(i))*(LMM->H_i0(i+1)); }
 
	
const RealArray1D
Caplet::	
nextControlledForwardPayoff() 
{
	// Libors X_j, j>=i until time T_t, t=min(i+1,n-1)
	LMM->newPath(t,i);     
	RealArray1D Z(2);
	
	Real X_iT_i,h,f;
    X_iT_i=LMM->XL(i,i);                     // Libor X_i(T_i)
    h=max(X_iT_i-delta_i*kappa,0.0);         // payoff at time T_{i+1}
	// H_{i+1}(T_{i+1}) accrual factor T_{i+1}->T_n at time T_{i+1}
    f=LMM->H_ii(i+1);                        
    Z[0] = f*h;		                         // forward payoff
    Z[1] = X_iT_i*(LMM->H_it(i+1,i));        // control variate 
	
	return Z;
} 


Real 
Caplet::
forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s)
{ 
	return theLattice->forwardCapletPayoff(i,kappa,node,s); 
}
    
       
Real 
Caplet::
analyticForwardPrice() const
{
    Real Li0=LMM->L_i0(i),                                           // L_i(0)
         cSigma=LMM->capletAggregateVolatility(i),                   // aggregate volatility to expiry    
	     Nplus=FinMath::N(FinMath::d_plus(Li0,kappa,cSigma)),
	     Nminus=FinMath::N(FinMath::d_minus(Li0,kappa,cSigma)),
	     f=LMM->H_i0(i+1);                                           // H_{i+1}(0)=B_{i+1}(0)/B_n(0)
 
    return delta_i*(Li0*Nplus-kappa*Nminus)*f;  
} 
   

ostream& 
Caplet::
printSelf(ostream& os) const
{ return os << "\nCaplet Cplt([T_"<<i<<",T_"<<i+1<<"],"<<kappa<<")"; }



/*******************************************************************************
 *
 *                          SWAPTIONS
 *
 ******************************************************************************/


Swaption::
Swaption(int p_, int q_, int t_, Real strike, LiborMarketModel* lmm) :
// Libors L_j, j>=t needed until time t  (forward transporting)
LiborDerivative(lmm,t_,true,true,true,true),                                          
p(p_), q(q_), t(t_), kappa(strike)
{   } 


Swaption* 
Swaption::
sample(int p, int q, int lmmType, int volType, int corrType)
{
	LiborMarketModel* lmm=LiborMarketModel::sample(q,lmmType,volType,corrType);
	int t=p;
	Real kappa=lmm->swapRate(p,q);
	return new Swaption(p,q,t,kappa,lmm);
} 
	 
    

Real 
Swaption::
nextForwardPayoff()
{
	// Libors X_j, j>=t until time T_t
	LMM->newPath(t,t);     	
	
	Real S_pqT,H_pqT;
	S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
	if(S_pqT<kappa) return 0.0;
    H_pqT=LMM->H_pq(p,q,t);
    return H_pqT*(S_pqT-kappa);                   // payoff at time T_t
}

	 
Real 
Swaption::
controlVariateMean() 
{ 
	Real fp=LMM->H_i0(p),        // H_p(0)
		 fq=LMM->H_i0(q);        // H_q(0)
    return fp-fq;
} 
 
	 
const RealArray1D
Swaption::
nextControlledForwardPayoff() 
{ 
	// Libors X_j, j>=t until time T_t
	LMM->newPath(t,t); 
	RealArray1D Z(2);
	
	Real S_pqT,H_pqT;
	S_pqT=LMM->swapRate(p,q,t);         // swaprate S_{p,q}(T_t)
	if(S_pqT<kappa) Z[0]=0.0;
	else{
		
       H_pqT=LMM->H_pq(p,q,t);
       Z[0]=H_pqT*(S_pqT-kappa);           // payoff at time T_t
	}
	
    Real fp,fq;
	fp=LMM->H_it(p,t),                  // H_p(T_m)
	fq=LMM->H_it(q,t);                  // H_q(T_m)
	Z[1]=fp-fq;                         // control variate
	
	return Z;
} 


Real 
Swaption::
forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s)
{ 
	return theLattice->forwardSwaptionPayoff(p,q,kappa,node,s); 
}


Real 
Swaption::
analyticForwardPrice() const
{
    Real S_pq=LMM->swapRate(p,q),                                    // S_{p,q}(0)
         swpnSigma=LMM->swaptionAggregateVolatility(p,q,t),          // aggregate volatility to T_t    
	     Nplus=FinMath::N(FinMath::d_plus(S_pq,kappa,swpnSigma)),
	     Nminus=FinMath::N(FinMath::d_minus(S_pq,kappa,swpnSigma)),
	     H_pq=LMM->H_pq(p,q);                                        // forward B_{p,q}(0)

    return H_pq*(S_pq*Nplus-kappa*Nminus);  
} 

 

ostream& 
Swaption::
printSelf(ostream& os) const
{  return os << "\nSwaption Swpn(T_"<<t<<",[T_"<<p<<",T_"<<q<<"],"<<kappa<<")"; }

	

/*******************************************************************************
 *
 *                   CALL ON BONDS
 *
 ******************************************************************************/



BondCall::
BondCall(Bond* D, Real strike, int s) : 
// Libors L_j, j>=s needed until time T_s (forward transporting)
LiborDerivative(D->getLMM(),s,true,true,true,true),     
B(D),
p(D->get_p()), q(D->get_q()), 
K(strike), t(s)
{   }
	
	
BondCall* 
BondCall::
sample(int p, int q, int lmmType, int volType, int corrType)
{
	LiborMarketModel* 
	lmm=LiborMarketModel::sample(q,lmmType,volType,corrType);
	RealArray1D c(q-p,p);
	for(int j=p;j<q;j++) c[j]=0.5+Random::U01();
	Bond* bond=new Bond(p,q,c,lmm);
	Real K=bond->cashPrice();
	
	int t=p;     // exercise at T_t
	return new BondCall(bond,K,t);
}
	 

BondCall* 
BondCall::
sampleCallOnZeroCouponBond(int p, int lmmType, int volType, int corrType)
{
	LiborMarketModel* 
	lmm=LiborMarketModel::sample(p+3,lmmType,volType,corrType);
	int t=p-1;
	Bond* B=new Bond(p,lmm);
	Real K=B->cashPrice();
		
	return new BondCall(B,K,t);
}

     
Real 
BondCall::
nextForwardPayoff() 
{ 
	// Libors X_j, j>=t until time T_t
	LMM->newPath(t,t); 	
	Real F=B->forwardPrice(t),
		 Ht=LMM->H_it(t,t);                // H_t(T_t)=1/B_n(T_t)
	return max(F-K*Ht,0.0);
}

	 
Real 
BondCall::
controlVariateMean() { return B->forwardPrice(); }


const RealArray1D
BondCall::
nextControlledForwardPayoff() 
{
	// Libors X_j, j>=t until time T_t
	LMM->newPath(t,t); 
	RealArray1D Z(2);
	
	Real F=B->forwardPrice(t),
		 Ht=LMM->H_it(t,t);                // H_t(T_t)=1/B_n(T_t)
	Z[0] = max(F-K*Ht,0.0);                // payoff
	Z[1] = B->forwardPrice(t);             // control variate
	
	return Z;
} 


Real 
BondCall::
forwardPayoff(StandardBrownianNode* node, LmmLattice* theLattice, int s)
{ 
	return theLattice->forwardBondCallPayoff(this,node,s); 
}



Real 
BondCall::
analyticForwardPrice() const
{
	Real F=B->forwardPrice(),
	     H_t0=LMM->H_i0(t),                                          // H_t(0)
	     Q=F/H_t0,
         bondSigma=LMM->bondAggregateVolatility(B,t),                // aggregate volatility to T_t    
	     Nplus=FinMath::N(FinMath::d_plus(Q,K,bondSigma)),
	     Nminus=FinMath::N(FinMath::d_minus(Q,K,bondSigma));
	                                             
    return F*Nplus-K*H_t0*Nminus;  
} 
   
   
  
   
ostream& 
BondCall::
printSelf(ostream& os) const
{  
   return
   os << "Call on bond B with strike K =: " << K << endl << *B;
}



MTGL_END_NAMESPACE(Martingale)