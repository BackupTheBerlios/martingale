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

#include "LiborCalibrator.h"
#include "TypedefsMacros.h"
#include "FinMath.h"
#include "LiborFactorLoading.h"
#include "Optimizer.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <cmath>

//#include <string>
#include <math.h>


MTGL_BEGIN_NAMESPACE(Martingale)





/**********************************************************************************
 *
 *               Swaptions and Caplets
 *
 *********************************************************************************/



/** Write swaption data to stream: p, q, strike, forward price.
 *  We don't write the fields calibratedForwardPrice and error
 *  since we wan to use this function to write synthetic
 *  calibration data.
 */
std::ostream& operator << (std::ostream& os, const SwaptionData& swpn)
{
	os << swpn.p << "  " << swpn.q << "  " << swpn.strike << "  " 
	   << swpn.forwardPrice;
    return os;
} // end operator <<


/** Read swaption data from stream: p, q, strike, forward price. */
std::istream& operator >> (std::istream& is, SwaptionData& swpn)
{
	is >> swpn.p >> swpn.q >> swpn.strike >> swpn.forwardPrice;
    return is;
} // end operator <<

 
 

/** Write caplet data to stream: i, strike, forward price.
 */
std::ostream& operator << (std::ostream& os, const CapletData& cplt)
{
	os << cplt.i << "  " << cplt.strike << "  " << cplt.forwardPrice << "  ";
    return os;
} // end operator <<


/** Read caplet data from stream: i, strike, forward price.
 */
std::istream& operator >> (std::istream& is, CapletData& cplt)
{
	is >> cplt.i >> cplt.strike >> cplt.forwardPrice;
    return is;
} // end operator <<





/**********************************************************************************
 *
 *                LmmCalibrator
 *
 *********************************************************************************/
	
	
// ACCRUAL FACTORS, BONDS, SWAP RATES, .. AT TIME ZERO
	
/** H_0(0) */
Real 
LmmCalibrator::
H0()
{
    Real f=1.0;
	for(int k=0;k<n;k++)f*=(1+x[k]);
	return f;
}


/** H_i(0) */
Real 
LmmCalibrator::
H_i0(int i)
{
    Real f=1.0;
	if(i==n) return f;
		
	for(int k=i;k<n;k++)f*=(1+x[k]);
	return f;
}
     

/** B_i(0) */
Real 
LmmCalibrator::
B0(int i)
{ 
    Real f=1.0;
     // accumulate 1 from time t=0 to time t=T_i                    
     for(int j=0;j<i;j++)f*=(1.0+x[j]); 
     return 1.0/f;                                 // B_i(0)=1/f  
}   
 

/** S_{pq}(0) */
Real 
LmmCalibrator::
swapRate(int p, int q)
{
	Real f=1.0+x[q-1], S=delta[q-1];
    for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+x[k]); }
 
    return (f-1.0)/S;
} 


/** B_{pq}(0) */
Real 
LmmCalibrator::
B_pq(int p, int q)
{ 
    Real S=0.0, F=B0(q);
    for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+x[k]); }
    return S;
} 
	 

/** H_pq(0) */
Real 
LmmCalibrator::
H_pq(int p, int q){ return B_pq(p,q)*H0(); } 



// CAPLET AND SWAPTION PRICES


Real 
LmmCalibrator::
capletForwardPrice(int i, Real strike)
{
   Real delta_i=delta[i], 
        Li0=factorLoading->getInitialLibors()[i],                   // L_i(0)
        cSigma=capletAggregateVolatility(i),                        // aggregate volatility to expiry    
        Nplus=FinMath::N(FinMath::d_plus(Li0,strike,cSigma)),
        Nminus=FinMath::N(FinMath::d_minus(Li0,strike,cSigma)),
        f=H_i0(i+1);                                                 
 
   return delta_i*(Li0*Nplus-strike*Nminus)*f; 
}

   

Real 
LmmCalibrator::
swaptionForwardPrice(int p, int q, Real strike) 
{
	Real S_pq=swapRate(p,q),                                         // S_{p,q}(0)
         swpnSigma=swaptionAggregateVolatility(p,q),                 // aggregate volatility to T_p    
         Nplus=FinMath::N(FinMath::d_plus(S_pq,strike,swpnSigma)),
         Nminus=FinMath::N(FinMath::d_minus(S_pq,strike,swpnSigma)),
         f=H_pq(p,q);                                                // forward B_{p,q}(0)
  
    return f*(S_pq*Nplus-strike*Nminus);  
} 
   

void 
LmmCalibrator::
writeSyntheticData()
{
   // write caplets
   for(int i=2;i<n;i++){
		   
	   Real strike = factorLoading->getInitialLibors()[i];  // L_i(0)
	   Real forwardPrice = capletForwardPrice(i,strike);
		   
	   CapletData* caplet = new CapletData();
	   caplet->i=i;
	   caplet->strike=strike;
	   caplet->forwardPrice=forwardPrice;
		   
	   capletsOut << *caplet << endl;
   }
	   
   // write swaptions
   for(int p=2;p<n-1;p++)
   for(int q=p+2;q<n+1;q++){
		   
	   Real strike = swapRate(p,q);                               // S_pq(0)
	   Real forwardPrice = swaptionForwardPrice(p,q,strike);
		   
	   SwaptionData* swaption = new SwaptionData();
	   swaption->p=p;
	   swaption->q=q;
	   swaption->strike=strike;
	   swaption->forwardPrice=forwardPrice;
		   
	   swaptionsOut << *swaption << endl;
   } 
} // end writeSyntheticData
      
	
// CALIBRATION 

Real 
LmmCalibrator::
meanRelativeCalibrationError()
{
   int k=0;
   Real sum=0.0;
	   
   	// iterator is pointer to pointer to SwaptionData
	std::list<SwaptionData*>::const_iterator theSwaption;
	for(theSwaption=swaptions.begin(); theSwaption!=swaptions.end(); ++theSwaption)			
	{
		SwaptionData* currentSwaption=*theSwaption;
		sum+=abs(currentSwaption->error);
		k++;
	}

	// iterator is pointer to pointer to CapletData
	std::list<CapletData*>::const_iterator theCaplet;
	for(theCaplet=caplets.begin(); theCaplet!=caplets.end(); ++theCaplet)			
	{
		CapletData* currentCaplet=*theCaplet;
		sum+=abs(currentCaplet->error);
		k++;
	}
		
	return sum/k;
} // end meanRelativeCalibrationError
   
   
Real 
LmmCalibrator::
squaredCapletError()
{
	Real sumSquares=0.0;

	// iterator is pointer to pointer to CapletData
	std::list<CapletData*>::const_iterator theCaplet;
	for(theCaplet=caplets.begin(); theCaplet!=caplets.end(); ++theCaplet)			
	{
		CapletData* currentCaplet=*theCaplet;
			
		int i=currentCaplet->i;
		Real forwardPrice=currentCaplet->forwardPrice,
		     strike=currentCaplet->strike,
		     calibratedForwardPrice=capletForwardPrice(i,strike),
		     error=forwardPrice-calibratedForwardPrice;
			
		sumSquares+=error*error;
		currentCaplet->calibratedForwardPrice=calibratedForwardPrice;
		currentCaplet->error=error/(forwardPrice+0.0000000001);
	}
	return sumSquares;
} // end squaredCapletError
	
	
Real 
LmmCalibrator::
squaredSwaptionError()
{
	Real sumSquares=0.0;

	// iterator is pointer to pointer to SwaptionData
	std::list<SwaptionData*>::const_iterator theSwaption;
	for(theSwaption=swaptions.begin(); theSwaption!=swaptions.end(); ++theSwaption)			
	{
		SwaptionData* currentSwaption=*theSwaption;
			
		int p=currentSwaption->p,
		    q=currentSwaption->q;
		Real forwardPrice=currentSwaption->forwardPrice,
		     strike=currentSwaption->strike,
		     calibratedForwardPrice=swaptionForwardPrice(p,q,strike),
		     error=forwardPrice-calibratedForwardPrice;
			
		sumSquares+=error*error;
		currentSwaption->calibratedForwardPrice=calibratedForwardPrice;
		currentSwaption->error=error/(forwardPrice+0.0000000001);
	}
	return sumSquares;
} // end squaredCapletError
	
	

Real 
LmmCalibrator::
objectiveFunction(Real* x)
{
	factorLoading->setParameters(x);
	Real err=squaredCapletError();
	err+=squaredSwaptionError();
		
	return err;
}
	
     
Real 
LmmCalibrator::
calibrate(int nVals)
{
	readCaplets();
	readSwaptions();
	// initial guess for parameter values
	RealArray1D u(n+7);
			
    // initial guess for caplet vol scaling factors k[i]
	// iterator is pointer to pointer to CapletData
//<--------------------------BUG------------------------------>
// Can't have zero vols in DriftlessLMM since all the vols go into 
// into the computation of the caplet price, leads to nonpositive definite matrix.
// Must interpolate u[i] for the missing i.!!!
	std::list<CapletData*>::const_iterator theCaplet;
	for(theCaplet=caplets.begin(); theCaplet!=caplets.end(); ++theCaplet)			
	{
		CapletData* currentCaplet=*theCaplet;
		int i = currentCaplet->i;
		Real T_i=factorLoading->getT(i),
		     L_i0=x[i]/delta[i],
		     price=currentCaplet->forwardPrice,
		     strike=currentCaplet->strike;

		Real vol=FinMath::blackImpliedAggregateCallVolatility(L_i0,strike,price);
			
		VolSurface* volSurface=factorLoading->getVolSurface();
		Real int_sgsg=volSurface->integral_sgsg(0,T_i,T_i,T_i);
	
		u[i]=sqrt(vol/int_sgsg);
	}
	
		
	// now set a,b,c,d; alpha,beta,r_oo depending on the type of correlations
	// and volatilities
	switch(factorLoading->getVolSurfaceType()){
			
		case VolSurface::JR :
		u[n]=0.1; u[n+1]=0.7; u[n+2]=2.0; u[n+3]=0.3; break;
			
		case VolSurface::M :
		u[n]=1.0; u[n+1]=1.0; u[n+2]=1.0; u[n+3]=1.5; break;
			
		default : u[n]=1.0; u[n+1]=1.0; u[n+2]=1.0; u[n+3]=1.0;
	}
		
	switch(factorLoading->getCorrelationType()){
			
		case Correlations::JR :
		u[n+4]=0.2; u[n+4]=1.0; u[n+4]=1.0; break;
			
		case Correlations::CS :
		u[n+4]=1.0; u[n+5]=1.0; u[n+6]=1.0; break;
			
		default : u[n+4]=1.5; u[n+5]=0.05; u[n+6]=0.2; 
	}	
		
		
	RealArray1D du(n+7);
	for(int i=0;i<n+7;i++) du[i]=0.01;
		 
		 
	SobolLiborCalibrationOptimizer optimizer(n,this,u,nVals,du);
	const RealArray1D& w=optimizer.search();
		
	std::cout << "\n\n\na = " << w[n]
	          << "\nb = " << w[n+1]
	          << "\nc = " << w[n+2]
	          << "\nd = " << w[n+3]
	          << "\nalpha = " << w[n+4]
	          << "\nbeta = " << w[n+5]
	          << "\nr_oo = " << w[n+6]
	          << "\n\n\nMean relative calibration error: "
	          << meanRelativeCalibrationError();
		
	writeCaplets();
	writeSwaptions();
		
} // end calibrate
	




/**********************************************************************************
 *
 *              DriftlessLMM-Calibrator
 *
 *********************************************************************************/
	
	
Real 
DriftlessLmmCalibrator::
capletAggregateVolatility(int i)
{ 
	 Real T_i=factorLoading->getT(i);
	 UTRRealMatrix& R=factorLoading->logLiborCovariationMatrix(i,n,0,T_i).utrRoot();
	 RealVector x(n-i,i);
	 x[i]=1;
	 Real f=H_i0(i+1);
	 for(int j=i+1;j<n;j++){ Real Uj0=H_i0(j)-H_i0(j+1); x[j]=-Uj0/f; }
	 x*=R.transpose();
	 return x.norm();
} // end capletAggregateVolatility
	 
	 
Real 
DriftlessLmmCalibrator::
swaptionAggregateVolatility(int p, int q)
{
	 const RealArray1D& T=factorLoading->getTenorStructure();
     UTRRealMatrix& Q=factorLoading->logLiborCovariationMatrix(p,n,0,T[p]).utrRoot();
		 
	 RealVector x(n-p,p);
	 Real denom1=H_i0(p)-H_i0(q),
	      denom2=0;
	 for(int j=p;j<q;j++) denom2+=delta[j]*H_i0(j+1);
		
	 Real Up0=H_i0(p)-H_i0(p+1); x[p]=Up0/denom1;		 
	 for(int j=p+1;j<q;j++){
			 
		 Real Uj0=H_i0(j)-H_i0(j+1);
	     x[j]=Uj0/denom1-(T[j]-T[p])*Uj0/denom2;
	 }		 
	 for(int j=q;j<n;j++){
			 
	     Real Uj0=H_i0(j)-H_i0(j+1);
		 x[j]=-(T[q]-T[p])*Uj0/denom2;
	 }
	 x*=Q.transpose();
	 return x.norm();
} // end swaptionAggregateVolatility


// SYNTHETIC DATA GENERATION
void 
DriftlessLmmCalibrator::
writeSyntheticDataSample(int n, int volType=VolSurface::JR, int corrType=Correlations::CS)
{
	 // string conversion
	 Int_t nn(n);   
	 string  vols,                                // "CONST", "JR", "M"
	         corrs,                               // "JR", "CS"
	         dim=nn.toString(),
	         capletOutFile, swaptionOutFile;
		 
	 LiborFactorLoading* fl=LiborFactorLoading::sample(n,volType,corrType);
	 vols=fl->getVolSurface()->volSurfaceType();
	 corrs=fl->getCorrelations()->correlationType();
	 // called "InstrumentIn" since we'll be reading in from that
	 capletOutFile="SyntheticData/CapletsIn-DL-dim"+dim+"-"+vols+"-"+corrs+".txt";
     swaptionOutFile="SyntheticData/SwaptionsIn-DL-dim"+dim+"-"+vols+"-"+corrs+".txt";
		 
	LmmCalibrator* cal=new DriftlessLmmCalibrator
	 (fl,"CapletsIn","SwaptionsIn",capletOutFile.c_str(),swaptionOutFile.c_str());
		 
	 cal->writeSyntheticData();
		
} // end writeSyntheticDataSample
		 

void 
DriftlessLmmCalibrator::
writeSyntheticDataSample()
{
     for(int n=20;n<60;n+=10)
     for(int volType=0;volType<3;volType++)
	 for(int corrType=0;corrType<2;corrType++) writeSyntheticDataSample(n,volType,corrType);
     std::cout << "\n\nDone.";
}
	 
	
void 
DriftlessLmmCalibrator::
testCalibration(int nVals)
{
	int n=50;
	const char *CapletsIn="CapletsIn-DL-dim50-JR-CS.txt",
	           *SwaptionsIn="SwaptionsIn-DL-dim50-JR-CS.txt";
		
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(n,VolSurface::M,Correlations::JR);
		
	LmmCalibrator* cal=new DriftlessLmmCalibrator(fl,CapletsIn,SwaptionsIn);
	cal->calibrate(nVals);
}
	


/**********************************************************************************
 *
 *              PredictorCorrectorLMM-Calibrator
 *
 *********************************************************************************/

	
// CAPLET AND SWAPTION AGGREGATE VOLS
	
Real 
PredictorCorrectorLmmCalibrator::
capletAggregateVolatility(int i)
{ 
     Real T_i=factorLoading->getT(i);
	 Real volsqr=factorLoading->integral_sgi_sgj_rhoij(i,i,0,T_i);
	 return sqrt(volsqr);
} 
	 
	 
Real 
PredictorCorrectorLmmCalibrator::
swaptionAggregateVolatility(int p, int q)
{
	  const RealArray1D& T=factorLoading->getTenorStructure();
      UTRRealMatrix& 
      Q=factorLoading->logLiborCovariationMatrix(p,q,0,T[p]).utrRoot();
      RealVector x_pq(q-p,p);
      for(int j=p;j<q;j++) x_pq[j]=(B0(j)-B0(j+1))/B_pq(p,q);
	  x_pq*=Q.transpose();
	  return x_pq.norm();
} 

	 
void 
PredictorCorrectorLmmCalibrator::
writeSyntheticDataSample(int n, int volType=VolSurface::JR, int corrType=Correlations::CS)
{
	 // string conversion
	 Int_t nn(n);   
	 string  vols,                             // "CONST", "JR", "M"
	         corrs,                            // "JR", "CS"
	         dim=nn.toString(),
	         capletOutFile, swaptionOutFile;
		 
	 LiborFactorLoading* fl=LiborFactorLoading::sample(n,volType,corrType);
	 vols=fl->getVolSurface()->volSurfaceType();
	 corrs=fl->getCorrelations()->correlationType();
	 // called "InstrumentIn" since we'll be reading in from that
	 capletOutFile="SyntheticData/CapletsIn-DL-dim"+dim+"-"+vols+"-"+corrs+".txt";
     swaptionOutFile="SyntheticData/SwaptionsIn-DL-dim"+dim+"-"+vols+"-"+corrs+".txt";
		 
	 LmmCalibrator* cal=new PredictorCorrectorLmmCalibrator
	 (fl,"CapletsIn","SwaptionsIn",capletOutFile.c_str(),swaptionOutFile.c_str());
		 
	 cal->writeSyntheticData();
} 
		 

void 
PredictorCorrectorLmmCalibrator::
writeSyntheticDataSample()
{
     for(int n=20;n<60;n+=10)
     for(int volType=0;volType<3;volType++)
	 for(int corrType=0;corrType<2;corrType++) writeSyntheticDataSample(n,volType,corrType);
     std::cout << "\n\nDone.";
}


/**********************************************************************************
 *
 *               LmmOptimizer
 *
 *********************************************************************************/


bool 
SobolLiborCalibrationOptimizer::
isInDomain(Real* x)
{
	Real alpha=0.5, beta=0.5, r_oo=0.5;
	for(int i=0;i<n;i++) if(x[i]<0.02) return false;
			
	switch(cal->getFactorLoading()->getCorrelationType()){
			
		case Correlations::CS : 
				
		    alpha=x[n+4], beta=x[n+5], r_oo=x[n+6];     // book, 6.11.1    
		    if((alpha<beta)||(beta<0.0)||(r_oo<0.0)) return false;
			
		case Correlations::JR : 
				
		    beta=x[n+4]; 
		    if(beta<0.0) return false;
	}
	return true;
}
	

Real 
SobolLiborCalibrationOptimizer::
f(Real* x){ return cal->objectiveFunction(x); }
	

			

MTGL_END_NAMESPACE(Martingale)				
