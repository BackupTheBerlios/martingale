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



MTGL_BEGIN_NAMESPACE(Martingale)


using std::ostream;
using std::cout;
using std::endl;



/*******************************************************************************
 *
 *                     OPTION
 *
 ******************************************************************************/


PathGenerator* 
Option::
getPathGenerator()
{
	cout << "\n\nGeneric Option::getPathGenerator(): "
	     << "not implemented in this generality."
	     << "\nAborting.";
	exit(0);
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
monteCarloForwardPrice(int nPaths)
{
    return Pricing::monteCarloForwardPrice(this,nPaths);
}


Real
Option::
controlledMonteCarloForwardPrice(int nPaths)
{
    return Pricing::controlledMonteCarloForwardPrice(this,nPaths);
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


void 
LiborPathGenerator::
newPath(){ LMM->newPath(t,i); }


LiborDerivative::
LiborDerivative
(LiborMarketModel* lmm, int s, int i, bool lt, bool mc, bool cv) : 
// T=lmm->getTenorStructure()[s] is the continuous time to expiry
Option(lmm->getTenorStructure()[s],lt,mc,cv),
LMM(lmm),
LPG(new LiborPathGenerator(lmm,s,i)), n(lmm->getDimension()), t(s)
{  }


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
	Real delta=LMM->getDeltas()[0];
	// size of time step: 0.05 or at most 250 time steps in all.
	Real dt = max(0.05, t*delta/250.0);
	int nSteps = (int) delta/dt;
	
	return new LmmLattice2F(fl,t,nSteps,false);
}


Real 
LiborDerivative::
latticeForwardPrice()
{
    LmmLattice* theLattice = getDefaultLattice();
	Real latticePrice = Option::latticeForwardPrice(theLattice);
	delete theLattice;
	return latticePrice;
}		


void 
LiborDerivative::
testPrice()
{ 
	 // all prices forward prices at time T_n
     Real aPrice,               // analytic price
          mcPrice,              // Monte carlo price
          cvmcPrice,            // Monte carlo price with control variate
	      lPrice,               // Price in the default lattice
	      epsilon=0.00000001;   // replacement for zero denominator
         
	  // message
	  cout << "\n\n\n\n\n" << *this << *LMM
           << "\nEffective dimension of the simulation = " << effectiveDimension()
	       << endl << endl
	       << "\nForward prices: "  << endl;
		
	  // PRICES: 
	
	  // analytic
      aPrice=analyticForwardPrice();
	  cout << "Analytic: " << aPrice << endl;
	
	  // lattice 
	  if(hasLattice_){
		  
	      lPrice=latticeForwardPrice();
	      cout << "Default lattice (driftless LMM only): " << lPrice << " , relative error: "
		       << relativeError(aPrice,lPrice,epsilon) << "%" << endl;
      }
		
	  // Monte Carlo
	  int nPath=10000;
	  
	  if(hasMonteCarlo_){
		  
          mcPrice=monteCarloForwardPrice(nPath);
	      cout << "Monte Carlo, " << nPath << " paths: "<< mcPrice << " , relative error: " 
		       << relativeError(aPrice,mcPrice,epsilon) << "%" << endl;
	  }
	  
	  if(hasControlVariate_){
		  
	      cvmcPrice=controlledMonteCarloForwardPrice(nPath);
          cout << "Monte Carlo with control variate: " << cvmcPrice << ", relative error: " 
		       << relativeError(aPrice,cvmcPrice,epsilon) << "%"<< endl;
	  }
		   
	  		   
      // QMC dynamics
	  LMM->switchToQMC();
	  
	  if(hasMonteCarlo_){
		  
          mcPrice=monteCarloForwardPrice(nPath);
          cout << "Quasi Monte Carlo: " << mcPrice << ", relative error: " 
		       << relativeError(aPrice,mcPrice,epsilon) << "%" << endl;
	  }
		   
	  LMM->restartSobolGenerator();
	  
	  if(hasControlVariate_){
		  
	      cvmcPrice=controlledMonteCarloForwardPrice(nPath);		   
		  cout << "Quasi Monte Carlo with control variate: " << cvmcPrice << ", relative error: " 
		   << relativeError(aPrice,cvmcPrice,epsilon) << "%"<< endl;
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
LiborDerivative(lmm,min(k+1,lmm->getDimension()-1),k,false,true,true),              
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
forwardPayoffAlongCurrentPath()
{                         
	Real X_iT_i,h,f;
    X_iT_i=LMM->XL(i,i);                     // Libor X_i(T_i)
    h=max(X_iT_i-delta_i*kappa,0.0);         // payoff at time T_{i+1}
    f=LMM->H_ii(i+1);                        // H_{i+1}(T_{i+1}) accrual T_{i+1}->T_n at time T_{i+1}
     // move this from time T_{i+1} to time T_n
    return f*h;	
}

	 

Real 
Caplet::
controlVariateMean(){ return (LMM->X_i0(i))*(LMM->H_i0(i+1)); }
 
	
Real
Caplet::	
controlVariateAlongCurrentPath() 
{
    Real X_iT_i=LMM->XL(i,i);             // Libor X_i(T_i)
    return X_iT_i*(LMM->H_it(i+1,i));         
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
// Libors L_j, j>=p needed until time t   
LiborDerivative(lmm,t_,p_,true,true,true),                                          
p(p_), q(q_), t(t_),                                  
kappa(strike)
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
forwardPayoffAlongCurrentPath()
{
	Real S_pqT,B_pqT,h,f;
	S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
    B_pqT=LMM->B_pq(p,q,t);
    h=B_pqT*max(S_pqT-kappa,0.0);                 // payoff at time T_t
    f=LMM->H_ii(t);                               // H_t(T_t) accrual T_t->T_n at time T_t
     // move this from time T_t to time T_n
    return f*h;	
}

	 
Real 
Swaption::
controlVariateMean() 
{ 
	Real fp=LMM->H_i0(p),        // H_p(0)
		 fq=LMM->H_i0(q);        // H_q(0)
    return fp-fq;
} 
 
	 
Real
Swaption::
controlVariateAlongCurrentPath() 
{ 
    Real fp,fq;
	fp=LMM->H_it(p,t),                            // H_p(T_t)
	fq=LMM->H_it(q,t);                            // H_q(T_t)
	return fp-fq;
} 


Real 
Swaption::
forwardPayoff(LmmNode* node)
{ 
	return node->forwardSwaptionPayoff(p,q,kappa); 
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
// Libors L_j, j>=D.p needed until time T_s
LiborDerivative(D->getLMM(),s,D->get_p(),true,true,true),     
B(D),
p(D->get_p()), q(D->get_q()), 
K(strike), t(s)
{   }
	
	
BondCall* 
BondCall::
sample(int n, int lmmType, int volType, int corrType)
{
	LiborMarketModel* 
	lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
	int p=n/3, q=2*n/3, t=p;
	RealArray1D c(q-p,p);
	for(int j=p;j<q;j++) c[j]=0.5+Random::U01();
	Bond* bond=new Bond(p,q,c,lmm);
	Real K=bond->cashPrice();
		
	return new BondCall(bond,K,t);
}
	 

BondCall* 
BondCall::
sampleCallOnZeroCouponBond(int n, int lmmType, int volType, int corrType)
{
	LiborMarketModel* 
	lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
	int i=n/2, t=i-1;
	Bond* B=new Bond(i,lmm);
	Real K=B->cashPrice();
		
	return new BondCall(B,K,t);
}

     
Real 
BondCall::
forwardPayoffAlongCurrentPath() 
{  
	Real F=B->forwardPrice(t),
		 Ht=LMM->H_ii(t);                // H_t(T_t)

	return max<Real>(F-K*Ht,0.0);
}
//<----------- INVESTIGATE ----------->
// if we use LMM instead of lmm=B->getLMM(), then the paths of the Libor model underlying B
// don't move even with our initilization of LiborDerivative::LMM as B->getLMM(), why is that???		

	 
Real 
BondCall::
controlVariateMean() { return B->forwardPrice(); }


Real
BondCall::
controlVariateAlongCurrentPath() 
{
    return B->forwardPrice(t); 
} 


Real 
BondCall::
forwardPayoff(LmmNode* node)
{	
	return node->forwardBondCallPayoff(this); 
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
   os << "Call on bond B with strike K =: " << K << *B;
}



MTGL_END_NAMESPACE(Martingale)