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
 


#include "Derivatives.h"
#include "FinMath.h"
#include "RandomObject.h"
#include "ControlledRandomVariable.h"
#include "PredictorCorrectorLMM.h"
#include "FastPredictorCorrectorLMM.h"
#include "DriftlessLMM.h"
#include "LowFactorDriftlessLMM.h"
#include "LiborMarketModel.h"          // Bond
#include <math.h>
#include <cmath>
#include <algorithm>
//#include <string>
#include <iostream>


MTGL_BEGIN_NAMESPACE(Martingale)


/*******************************************************************************
 *
 *                     DERIVATIVES
 *
 ******************************************************************************/


RandomVariable* 
Derivative::
forwardPayoff() { return new ForwardPayoff(this); } 
     

    
Real 
Derivative::
controlVariateMean() 
{
      std::cout << "Derivative.controlVariateMean():"
		   << endl << "no control variate implemented, aborting.";
      exit(0);
      return 0.0;    // keeps the compiler happy
} 
      
  

RealVector 
Derivative::
nextControlledForwardPayoff() 
{
     std::cout << "Derivative.nextControlledForwardPayoff():"
		  << endl << "no control variate implemented, aborting.";
     exit(0);
	 // keep the compiler happy
     RealVector v(2); 
     return v;    
} 
      

ControlledRandomVariable* 
Derivative::
controlledForwardPayoff() { return new ControlledForwardPayoff(this); } 

     
Real 
Derivative::
analyticForwardPrice() const
{
     std::cout << endl 
          << "Derivative.analyticForwardPrice():" << endl
          << "no analytic price implemented, aborting.";
     exit(0);
     return 0.0;    // keeps the compiler happy
}
      
      
Real 
Derivative::
monteCarloForwardPrice(int nPath){ return forwardPayoff()->expectation(nPath); }

	  
Real 
Derivative::
monteCarloForwardPrice(int nPath, string message){ return forwardPayoff()->expectation(nPath,message); }


Real 
Derivative::
controlledMonteCarloForwardPrice(int nPath){ return controlledForwardPayoff()->expectation(nPath); }


Real 
Derivative::ForwardPayoff::
nextValue() 
{ 
	return LD->nextForwardPayoff(); 
}
      
     
RealVector 
Derivative::ControlledForwardPayoff::
nextValue() 
{ 
	return LD->nextControlledForwardPayoff(); 
}   
		 

Real 
Derivative::ControlledForwardPayoff::
getControlVariateMean() 
{ 
	return LD->controlVariateMean(); 
}



// GLOBAL INSERTION

std::ostream& operator << 
(std::ostream& os, const Derivative& fl){ return fl.printSelf(os); }


/*******************************************************************************
 *
 *                          LIBOR DERIVATIVE
 *
 ******************************************************************************/


LiborDerivative::
LiborDerivative(LiborMarketModel* lmm, int t) : 
LMM(lmm), n(lmm->getDimension()), horizon(t) 
{  }
	

int 
LiborDerivative::
effectiveDimension() const { return LMM->effectiveDimension(0, horizon); }
	
	
std::ostream& 
LiborDerivative::
printSelf(std::ostream& os) const { return os << "Generic Libor Derivative"; } 
	
		

void 
LiborDerivative::
testPrice()
{ 
	 // all prices forward prices at time T_n
     Real aprice,               // analytic price
          mcprice,              // Monte carlo price
          cvmcprice,            // Monte carlo price with control variate
	      epsilon=0.00000001;   // replacement for zero denominator in relError
         
	 // correlation with control variate
	std::cout << "\n\n\n\n";
	printSelf(std::cout);                       // cout << *this linkage problems
	LMM->printSelf(std::cout);                  // cout << *LMM linkage problems
    std::cout << endl << "Effective dimension of the simulation = " 
	          << effectiveDimension()
	          << endl;
		
	  int nPath=20000;
      // prices
      aprice=analyticForwardPrice();
      mcprice=monteCarloForwardPrice(nPath);
	  cvmcprice=controlledMonteCarloForwardPrice(nPath);

      std::cout << "\nForward prices, " << nPath << " paths:" << endl
		   << "Analytic: " << aprice << endl
		   << "Monte Carlo: " << mcprice << " , relative error: " 
		   << relativeError(aprice,mcprice,epsilon) << "%" << endl
		   << "Monte Carlo with control variate: " << cvmcprice << ", relative error: " 
		   << relativeError(aprice,cvmcprice,epsilon) << "%"
		   << endl;
		   
      // QMC dynamics
	  LMM->switchToQMC();
		 
      // prices
      mcprice=monteCarloForwardPrice(nPath);
	  LMM->restartSobolGenerator();
	  cvmcprice=controlledMonteCarloForwardPrice(nPath);
 
      std::cout << "Quasi Monte Carlo, relative error: " 
		   << relativeError(aprice,mcprice,epsilon) << "%" << endl
		   << "Quasi Monte Carlo with control variate, relative error: " 
		   << relativeError(aprice,cvmcprice,epsilon) << "%"
		   << endl;

} // end testPrice



void 
LiborDerivative::
priceTest(int lmmType, int volType, int corrType) 
{
     Timer watch;
     setLMM(LiborMarketModel::sample(n,lmmType,volType,corrType)); 
     testPrice();
	 watch.stop();
	 watch.report(" ");

} // end priceTest
	


/*******************************************************************************
 *
 *                          CAPLETS
 *
 ******************************************************************************/


void 
Caplet::
setLMM(LiborMarketModel* lmm)
{ 
	LiborDerivative::setLMM(lmm); 
	// we must now adjust delta_i
	delta_i=(lmm->getDeltas())[i];
}
	

Caplet::
Caplet(int k, Real strike, LiborMarketModel* lmm) :
LiborDerivative(lmm,min(k+1,lmm->getDimension()-1)),               // Libors needed until time t=min(k+1,n-1)
i(k), 
kappa(strike),
delta_i((lmm->getDeltas())[i])
{   
	armControlVariate();
} 


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
	LMM->newPath(horizon,i);                      // needed Libors L_j, j>=i.
	Real X_iT_i,h,f;
    X_iT_i=LMM->XL(i,i);                          // Libor X_i(T_i)
    h=max(X_iT_i-delta_i*kappa,0.0);        // payoff at time T_{i+1}
    f=LMM->H_ii(i+1);                             // H_{i+1}(T_{i+1}) accrual T_{i+1}->T_n at time T_{i+1}
     // move this from time T_{i+1} to time T_n
    return f*h;	
}

	 

Real 
Caplet::
controlVariateMean(){ return (LMM->X_i0(i))*(LMM->H_i0(i+1)); }
 
	
RealVector 
Caplet::	
nextControlledForwardPayoff() 
{
	LMM->newPath(horizon,i);            // needed Libors L_j, j>=i
		 
    Real X_iT_i=LMM->XL(i,i),                        // Libor X_i(T_i)
		 h=max(X_iT_i-delta_i*kappa,0.0),      // payoff at time T_{i+1}
         f=LMM->H_ii(i+1);                           // H_{i+1}(T_{i+1}), accrual T_{i+1}->T_n at time T_{i+1}      
        
    RealVector v(2); 
    v[0]=h*f;                                        // forward transported payoff
    v[1]=X_iT_i*(LMM->H_it(i+1,i));                  // control variate
    return v;
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
   

std::ostream& 
Caplet::
printSelf(std::ostream& os) const
{ return os << "\nCaplet Cplt([T_"<<i<<",T_"<<i+1<<"],"<<kappa<<")"; }



/*******************************************************************************
 *
 *                          SWAPTIONS
 *
 ******************************************************************************/


Swaption::
Swaption(int _p, int _q, int _t, Real strike, LiborMarketModel* lmm) :
LiborDerivative(lmm,_t),                           // Libors needed until time T_t                  
p(_p), q(_q), t(_t),                                  
kappa(strike)
{   
	armControlVariate();
} 


Swaption* 
Swaption::
sample(int n, int lmmType, int volType, int corrType)
{
	LiborMarketModel* lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
	int p=n/3, q=n, t=p;
	Real kappa=lmm->swapRate(p,q);
	return new Swaption(p,q,t,kappa,lmm);
} 
	 
    

Real 
Swaption::
nextForwardPayoff()
{
    LMM->newPath(horizon,t);                      // needed Libors L_j, j>=t
	Real S_pqT,B_pqT,h,f;
	S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
    B_pqT=LMM->B_pq(p,q,t);
    h=B_pqT*max<Real>(S_pqT-kappa,0.0);           // payoff at time T_t
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
 
	 
RealVector 
Swaption::
nextControlledForwardPayoff() 
{
	LMM->newPath(horizon,t);            // needed Libors L_j, j>=t
		 
    Real S_pqT,B_pqT,h,ft,fp,fq;
	S_pqT=LMM->swapRate(p,q,t);                   // swaprate S_{p,q}(T_t)
	B_pqT=LMM->B_pq(p,q,t);
	h=B_pqT*max<Real>(S_pqT-kappa,0.0);           // payoff at time T_t
	ft=LMM->H_ii(t),                              // H_t(T_t) accrual to horizon from time t
	fp=LMM->H_it(p,t),                            // H_p(T_t)
	fq=LMM->H_it(q,t);                            // H_q(T_t)
        
    RealVector v(2); 
	v[0]=h*ft;                                     // forward transported payoff
	v[1]=fp-fq;                                    // control variate
	return v;
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
   

std::ostream& 
Swaption::
printSelf(std::ostream& os) const
{  return os << "\nSwaption Swpn(T_"<<t<<",[T_"<<p<<",T_"<<q<<"],"<<kappa<<")"; }

	

/*******************************************************************************
 *
 *                   CALL ON BONDS
 *
 ******************************************************************************/


void 
BondCall::	
setLMM(LiborMarketModel* lmm)
{ 
	LiborDerivative::setLMM(lmm); 
	// we must also reset the libor marke model for the underlying bond
	B->setLMM(lmm);      
}
	

BondCall::
BondCall(Bond* D, Real strike, int s) : 
LiborDerivative(D->getLMM(),s),      // Libors needed until time T_s
B(D),
p(D->get_p()), q(D->get_q()), 
K(strike), t(s)
{   
	armControlVariate();
}
	
	
BondCall* 
BondCall::
sample(int n, int lmmType, int volType, int corrType)
{
	LiborMarketModel* 
	lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
	int p=n/3, q=2*n/3, t=p;
	RealArray1D c(q-p,p);
	for(int j=p;j<q;j++) c[j]=0.5+Random::U01();
	Bond* B=new Bond(p,q,c,lmm);
	Real K=B->cashPrice();
		
	return new BondCall(B,K,t);
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
nextForwardPayoff() 
{  
	LiborMarketModel* lmm=B->getLMM();
	lmm->newPath(horizon,t);             // needed Libors L_j, j>=t
	Real F=B->forwardPrice(t),
		 Ht=lmm->H_ii(t);                // H_t(T_t)

	return max<Real>(F-K*Ht,0.0);
}
//<----------- INVESTIGATE ----------->
// if we use LMM instead of lmm=B->getLMM(), then the paths of the Libor model underlying B
// don't move even with our initilization of LiborDerivative::LMM as B->getLMM(), why is that???		

	 
Real 
BondCall::
controlVariateMean() { return B->forwardPrice(); }


RealVector 
BondCall::
nextControlledForwardPayoff() 
{
	LiborMarketModel* lmm=B->getLMM();
    lmm->newPath(horizon,t);              // needed Libors L_j, j>=t
	Real F=B->forwardPrice(t),
		 Ht=lmm->H_ii(t);                 // H_t(T_t)
        
    RealVector v(2); 
    v[0]=max<Real>(F-K*Ht,0.0);	       // forward transported payoff
    v[1]=B->forwardPrice(t);              // control variate
    return v;
} 
		 
   
Real 
BondCall::
analyticForwardPrice() const
{
    LiborMarketModel* lmm=B->getLMM();
	Real F=B->forwardPrice(),
	     H_t0=lmm->H_i0(t),                                          // H_t(0)
	     Q=F/H_t0,
         bondSigma=LMM->bondAggregateVolatility(B,t),                // aggregate volatility to T_t    
	     Nplus=FinMath::N(FinMath::d_plus(Q,K,bondSigma)),
	     Nminus=FinMath::N(FinMath::d_minus(Q,K,bondSigma));
	                                             
    return F*Nplus-K*H_t0*Nminus;  
} 
   
   
  
   
std::ostream& 
BondCall::
printSelf(std::ostream& os) const
{  
   os << "Call on bond B with strike K =: " << K << endl;
   return B->printSelf(os);                                         // linkage problems with os << *B
}



MTGL_END_NAMESPACE(Martingale)