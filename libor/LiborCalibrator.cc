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
#include "LiborMarketModel.h"
#include "Optimizer.h"
#include "Utils.h"
#include <iostream>
#include <list>
#include <string>
#include <cmath>
#include "QuasiMonteCarlo.h"

using std::ostream;
using std::istream;
using std::abs;



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
ostream& operator << (ostream& os, const SwaptionData& swpn)
{
	os << swpn.p << "  " << swpn.q << "  " << swpn.strike << "  " 
	   << swpn.forwardPrice;
    return os;
} // end operator <<


/** Read swaption data from stream: p, q, strike, forward price. */
istream& operator >> (istream& is, SwaptionData& swpn)
{
	is >> swpn.p >> swpn.q >> swpn.strike >> swpn.forwardPrice;
    return is;
} // end operator <<

 
 

/** Write caplet data to stream: i, strike, forward price.
 */
ostream& operator << (ostream& os, const CapletData& cplt)
{
	os << cplt.i << "  " << cplt.strike << "  " << cplt.forwardPrice << "  ";
    return os;
} // end operator <<


/** Read caplet data from stream: i, strike, forward price.
 */
istream& operator >> (istream& is, CapletData& cplt)
{
	is >> cplt.i >> cplt.strike >> cplt.forwardPrice;
    return is;
} // end operator <<





/**********************************************************************************
 *
 *                LmmCalibrator
 *
 *********************************************************************************/

LmmCalibrator::
LmmCalibrator
(LiborFactorLoading* fl,
 const char* capletsInFile="CapletsIn.txt",  const char* swaptionsInFile="SwaptionsIn.txt",
 const char* capletsOutFile="CapletsOut.txt", const char* swaptionsOutFile="SwaptionsOut.txt"
) :
factorLoading(fl),
n(fl->getDimension()),
X0(fl->getInitialXLibors()),
delta(fl->getDeltas()),
U0(n),
H0(n+1),
Bn0(0.0),
S_pq0(n+1,1),
capletsIn(capletsInFile),
swaptionsIn(swaptionsInFile),
logFile("CalibrationLog.txt"),
capletsOut(capletsOutFile),
swaptionsOut(swaptionsOutFile),
caplets(n-1,1),
swaptions(n-1,1)
{
	// H_i(0)
	H0[n]=1.0;
	for(int i=0;i<n;i++){
		
	   Real f=1.0;
	   for(int k=i;k<n;k++)f*=(1+X0[k]);
       H0[i]=f;
   }

   // U_i(0)
   for(int i=0;i<n;i++) U0[i]=H0[i]-H0[i+1];
	   
   // B_n(0)
   Real f=1.0;
   for(int k=0;k<n;k++) f*=(1+X0[k]);
   Bn0=1.0/f;
 
   // swap rates S_pq(0)
   for(int p=1;p<n;p++)
   for(int q=p+1;q<=n;q++)   
   {
	   
	   	Real f=1.0+X0[q-1], S=delta[q-1];
        for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+X0[k]); }
 
        S_pq0(p,q)=(f-1.0)/S;
	}
	  	
} // end constructor



// IO

void 
LmmCalibrator::
readCaplets()
{ 
	read<CapletData>(capletsIn,caplets); 
	logFile << "Read caplets:\n\n";
	write<CapletData>(logFile,caplets); 
	
}
	

void 
LmmCalibrator::
readSwaptions()
{ 
	read<SwaptionData>(swaptionsIn,swaptions);
	logFile << "\n\n\nRead swaptions:\n\n";
	write<SwaptionData>(logFile,swaptions); 
}


void 
LmmCalibrator::
recordCapletCalibration()
{
   for(int i=1;i<n;i++){
		   
	   // Caplet(i), price under current factorloaing parameters
	   Real strike = factorLoading->getInitialLibors()[i];  // L_i(0)
	   Real calPrice = capletForwardPrice(i,strike);	   
		   
	   CapletData* caplet = caplets[i];
	   caplet->calibratedForwardPrice=calPrice;
	   Real price=caplet->forwardPrice,
	        error=100.0*abs(price-calPrice)/price;
	   caplet->error=error;
   }
}
	
	
void 
LmmCalibrator::
writeCaplets()
{ 
	capletsOut << "i   strike    forwardPrice    calibratedForwardPrice    error(%)\n\n";
	recordCapletCalibration();
	write<CapletData>(capletsOut,caplets); 
}
	
	
void 
LmmCalibrator::
writeSwaptions()
{ 
	capletsOut << "p   q   strike    forwardPrice    calibratedForwardPrice    error(%)\n\n";
	write<SwaptionData>(swaptionsOut,swaptions); 
}



// CORRELATIONS

Real
LmmCalibrator::
rho(int i, int j){ return factorLoading->getRho()(i,j); }
	
// ACCRUAL FACTORS, BONDS, SWAP RATES, .. AT TIME ZERO

Real 
LmmCalibrator::
L_i0(int i){ return factorLoading->getInitialLibors()[i]; } 


/** H_i(0) */
Real 
LmmCalibrator::
H_i0(int i){  return H0[i];  }     // H_i=B_i/B_n  


/** H_{pq}(0) */
Real 
LmmCalibrator::
H_pq0(int p, int q)
{  
	Real sum=0;
	for(int k=p;k<q;k++) sum+=delta[k]*H0[k+1];
	return sum;
}


/** B_i(0) */
Real 
LmmCalibrator::
B_i0(int i){  return H0[i]*Bn0;  }     // H_i=B_i/B_n 


/** B_{pq}(0)=H_{pq}(0)*B_n(0) */
Real 
LmmCalibrator::
B_pq0(int p, int q){  return H_pq0(p,q)*Bn0; }


/** S_{pq}(0) */
Real 
LmmCalibrator::
swapRate(int p, int q){ return S_pq0(p,q); }



// CAPLET AND SWAPTION PRICES


Real 
LmmCalibrator::
capletForwardPrice(int i, Real strike, Real Sigma)
{
   Real delta_i=delta[i], 
        Li0=factorLoading->getInitialLibors()[i],                   // L_i(0)  
        Nplus=FinMath::N(FinMath::d_plus(Li0,strike,Sigma)),
        Nminus=FinMath::N(FinMath::d_minus(Li0,strike,Sigma)),
        f=H_i0(i+1);                                                 
 
   return delta_i*(Li0*Nplus-strike*Nminus)*f; 
}

   

Real 
LmmCalibrator::
swaptionForwardPrice(int i, Real strike, Real Sigma) 
{
	Real S_pq=swapRate(i,n),                                                                                
         Nplus=FinMath::N(FinMath::d_plus(S_pq,strike,Sigma)),
         Nminus=FinMath::N(FinMath::d_minus(S_pq,strike,Sigma)),
         f=H_pq0(i,n);                                                // forward B_{i,n}(0)
  
    return f*(S_pq*Nplus-strike*Nminus);  
} 

   

void 
LmmCalibrator::
writeSyntheticData()
{
   for(int i=1;i<n;i++){
		   
	   // Caplet(i)
	   Real strike = factorLoading->getInitialLibors()[i];  // L_i(0)
	   Real forwardPrice = capletForwardPrice(i,strike);
		   
	   CapletData* caplet = new CapletData();
	   caplet->i=i;
	   caplet->strike=strike;
	   caplet->forwardPrice=forwardPrice;
		   
	   capletsOut << *caplet << endl;

	   // Swaption(i,n)
	   strike = swapRate(i,n);                               
	   forwardPrice = swaptionForwardPrice(i,strike);
		   
	   SwaptionData* swaption = new SwaptionData();
	   swaption->p=i;
	   swaption->q=n;
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
   Real sum=0.0;	
   for(int i=1;i<n;i++){  
	
	    SwaptionData* currentSwaption=swaptions[i];
		sum+=abs(currentSwaption->error);

        //CapletData* currentCaplet=caplets[i];
		//sum+=abs(currentCaplet->error);
	}
		
	return sum/(2*(n-1));
} // end meanRelativeCalibrationError
   	




/**********************************************************************************
 *
 *              StandardLMM-Calibrator
 *
 *********************************************************************************/

StandardLmmCalibrator::
StandardLmmCalibrator
(LiborFactorLoading* fl,
 const char* capletsInFile="CapletsIn.txt", 
 const char* swaptionsInFile="SwaptionsIn.txt", 
 const char* capletsOutFile="CapletsOut.txt", 	 
 const char* swaptionsOutFile="SwaptionsOut.txt"
) :	
LmmCalibrator(fl,capletsInFile,swaptionsInFile,capletsOutFile,swaptionsOutFile),
capletImpliedSigma(fl->getDimension()-1,1),
cvMatrix(fl->getDimension()), 
x(fl->getDimension())
{    }



void 
StandardLmmCalibrator::
writeCapletImpliedSigmas()
{
    for(int i=1;i<n;i++){	
		
		CapletData* currentCaplet=caplets[i];
		Real L_i0=X0[i]/delta[i],
		     forwardPrice=currentCaplet->forwardPrice,
		     strike=currentCaplet->strike;
		
		// the blackImpliedAggregateCallVolatility solves the blackScholesFunction.
		// the capletForwardPrice multiplies this function with delta_i*H_i0(i+1)
		// we must divide by this to isolate the blackScholesFunction
		Real price=forwardPrice/(delta[i]*H_i0(i+1));
		capletImpliedSigma[i]=FinMath::blackImpliedAggregateCallVolatility(L_i0,strike,price); 
    }
}



void 
StandardLmmCalibrator::
writeCovariationMatrix(int p, int q)
{
	for(int i=p;i<q;i++)
	for(int j=i;j<q;j++){
		
		Real T_p=factorLoading->getT(p);
		cvMatrix(i,j)=factorLoading->integral_sgi_sgj_rhoij(i,j,0.0,T_p);
	}
}



Real
StandardLmmCalibrator::
objectiveFunction()
{
	Real strike, marketPrice, calibratedPrice, error=0.0, diff;
	// loop through Swaption(i,n), caplets matched exactly later
	for(int i=1;i<n;i++){ 
		
		SwaptionData* theSwaption=swaptions[i];
		marketPrice=theSwaption->forwardPrice; 
		strike=theSwaption->strike;
		calibratedPrice=swaptionForwardPrice(i,strike); 
		diff=calibratedPrice-marketPrice;
		error+=diff*diff;
		
		// write the information into the swaption
		(theSwaption->calibratedForwardPrice)=calibratedPrice;
		(theSwaption->error)=100*diff/marketPrice;
	}
	return error;
} // end objectiveFunction()
		


Real 
StandardLmmCalibrator::
objectiveFunction(const RealArray1D& x)
{
	factorLoading->setParameters(x);    // volsurface, correlations
	setScalingFactors();
	return objectiveFunction();
}




// Square root of xCx' with C=cvMatrix, x the field x, indices j,k=i,...,n-1.
Real
StandardLmmCalibrator::	
root_xCx(int i)
{
	UTRRealMatrix& C=cvMatrix;
//cerr << C;
	Real s=0.0, d=0.0;
	// use symmetry -- C only has the upper half of the covariance matrix.
	// diagonal:
	for(int j=i;j<n;j++) d+=x[j]*C(j,j)*x[j]; 
	// above diagonal
	for(int j=i;j<n;j++)
	for(int k=j+1;k<n;k++) s+=x[j]*C(j,k)*x[k];
    return sqrt(d+2*s);
}




// SYNTHETIC DATA GENERATION
void 
StandardLmmCalibrator::
writeSyntheticDataSample
(int n, int volType=VolSurface::JR, int corrType=Correlations::CS)
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
	 delete cal;
		
	 // called "InstrumentIn" since we'll be reading in from that
	 capletOutFile="SyntheticData/CapletsIn-PC-dim"+dim+"-"+vols+"-"+corrs+".txt";
     swaptionOutFile="SyntheticData/SwaptionsIn-PC-dim"+dim+"-"+vols+"-"+corrs+".txt";
	 
	 cal=new PredictorCorrectorLmmCalibrator
	 (fl,"CapletsIn","SwaptionsIn",capletOutFile.c_str(),swaptionOutFile.c_str());
	 
	 cal->writeSyntheticData();
	
	 
		
} // end writeSyntheticDataSample
		 


void 
StandardLmmCalibrator::
writeSyntheticDataSample()
{
     for(int n=20;n<60;n+=10)
     for(int volType=0;volType<3;volType++)
	 for(int corrType=0;corrType<2;corrType++) writeSyntheticDataSample(n,volType,corrType);
     std::cout << "\n\nDone.";
}
	 
	
void 
StandardLmmCalibrator::
testCalibration
(int nVals, int dim, 
 int dataLmmType, int dataVolType, int dataCorrType,
 int lmmType, int volType, int corrType)
{
	if(!((dim==20)||(dim==30)||(dim==40)||(dim==50))) dim=20; 
	// string conversion
	 Int_t Dim(dim);   
	 string  lmm, lmm_c,                                // "DL", "PC"
	         vols,  vols_c,                             // "CONST", "JR", "M"
	         corrs, corrs_c,                            // "JR", "CS"
	         dim_str=Dim.toString(),
	         capletInFile, swaptionInFile;
		    
	 // what the data came from
	 lmm=LiborMarketModel::lmmType(dataLmmType);
	 vols=VolSurface::volSurfaceType(dataVolType);
	 corrs=Correlations::correlationType(dataCorrType);
	
	 // what we calibrate
	 lmm_c=LiborMarketModel::lmmType(lmmType);
	 vols_c=VolSurface::volSurfaceType(volType);
	 corrs_c=Correlations::correlationType(corrType);

	 // the data files
	 capletInFile="SyntheticData/CapletsIn-"+lmm+"-dim"+dim_str+"-"+vols+"-"+corrs+".txt";
     swaptionInFile="SyntheticData/SwaptionsIn-"+lmm+"-dim"+dim_str+"-"+vols+"-"+corrs+".txt";
	 
	
	 LiborFactorLoading* fl=LiborFactorLoading::sample(dim,volType,corrType);

	 LmmCalibrator* cal;
	 switch(lmmType){
		 
		 case LiborMarketModel::PC  : 
			       cal=new PredictorCorrectorLmmCalibrator
		           (fl,capletInFile.c_str(),swaptionInFile.c_str(),"CapletsOut.txt","SwaptionsOut.txt");
		           std::cout << "\n\n\nCalibrating predictor-corrector LMM:"
		                     << "\nVolSurface: " << vols_c
		                     << "\nCorrelations: " << corrs_c
		                     << "\nDimension: " << dim
		                     << "\nto data "+capletInFile+", "+swaptionInFile;
		           break;
		 // default calibration to driftless LMM
		 default : cal=new DriftlessLmmCalibrator
		           (fl,capletInFile.c_str(),swaptionInFile.c_str(),"CapletsOut.txt","SwaptionsOut.txt");
		           std::cout << "\n\n\nCalibrating driftless LMM:"
		                     << "\nVolSurface: " << vols_c
		                     << "\nCorrelations: " << corrs_c
		                     << "\nDimension: " << dim
		                     << "\nto data "+capletInFile+", "+swaptionInFile;
	 }
	 cal->calibrate(nVals);
}



LiborFactorLoading* 
StandardLmmCalibrator::
calibrate(int nVals)
{
	readCaplets();
	readSwaptions();
	writeCapletImpliedSigmas();

	RealArray1D u(7);					
	// now set a,b,c,d; alpha,beta,r_oo depending on the type of correlations
	// and volatilities
	switch(factorLoading->getType()->volType){
			
		case VolSurface::JR :
		u[0]=0.1; u[1]=0.7; u[2]=2.0; u[3]=0.3; break;
			
		case VolSurface::M :
		u[0]=1.0; u[1]=1.0; u[2]=1.0; u[3]=1.5; break;
			
		default : u[0]=1.0; u[1]=1.0; u[2]=1.0; u[3]=1.0;
	}
		
	switch(factorLoading->getType()->corrType){
		
		case Correlations::JR :	
		u[4]=1.0; u[5]=0.2; u[6]=1.0; break;
			
		case Correlations::CS :
		u[4]=1.5; u[5]=0.05; u[6]=0.2; break;
			
		default : u[4]=1.0; u[5]=1.0; u[6]=1.0; 
	}	



	RealArray1D du(7);
	for(int i=0;i<7;i++) du[i]=0.02;
				 		 
	SobolLiborCalibrationOptimizer optimizer(this,u,nVals,du);
	const RealArray1D& w=optimizer.search();
	
	writeCaplets();
	writeSwaptions();
	
	Real meanError=meanRelativeCalibrationError();

	// report optimal parmeters
	std::cout << "\n\n\n\nOptimal parameters: " << endl << endl;
	std::cout << "\n\n\na = " << w[0]
	          << "\nb = " << w[1]
	          << "\nc = " << w[2]
	          << "\nd = " << w[3]
	          << "\nalpha = " << w[4]
	
	          << "\nbeta = " << w[5]
	          << "\nr_oo = " << w[6]
	          << "\n\n\nMean relative calibration error: "
	          << meanError;
	
    return factorLoading;
		
} // end calibrate






/**********************************************************************************
 *
 *              DriftlessLMM-Calibrator
 *
 *********************************************************************************/


// CAPLET AND SWAPTION FORWARD PRICES

Real 
DriftlessLmmCalibrator::
capletForwardPrice(int i, Real strike)
{
	// get the matrix R for both Caplet(i) and Swaption(i,n)
	writeCovariationMatrix(i,n);
	
	// Caplet(p), the vector x (book 6.8.2)
	x[i]=1.0;
	for(int j=i+1;j<n;j++) x[j]=-U0[j]/H0[i+1];
	
    // caplet volatility = sqrt(xCx')
	Real Sigma=root_xCx(i);
	return LmmCalibrator::capletForwardPrice(i,strike,Sigma);
}


Real 
DriftlessLmmCalibrator::
swaptionForwardPrice(int i, Real strike)
{
	// get the matrix R for both Caplet(i) and Swaption(i,n)
	writeCovariationMatrix(i,n);
	
	// Swaption(p,n), the vector x (book 6.8.3)
    const RealArray1D& T=factorLoading->getTenorStructure();
	Real denom1=H0[i]-1.0,                                      // 1.0=H_n(0)
	     denom2=0.0;
	for(int j=i;j<n;j++) denom2+=delta[j]*H0[j+1];
		
	x[i]=U0[i]/denom1;		 
	for(int j=i+1;j<n;j++) x[j]=U0[j]/denom1-(T[j]-T[i])*U0[j]/denom2;
	 		 
	// swaption volatility = sqrt(xCx')
    Real Sigma=root_xCx(i);
	return LmmCalibrator::swaptionForwardPrice(i,strike,Sigma);
}



void 
DriftlessLmmCalibrator::
setScalingFactors()
{	
	RealArray1D& c=factorLoading->getScalingFactors();
	// recursive computation starting from i=n-1, book, 6.11.3
    for(int i=n-1;i>0;i--){	
			
		Real T_i=factorLoading->getT(i);				
		// the B_{jk}, book, 6.11.1, upper half only
	    UTRRealMatrix B(n-1,i);
	    for(int j=i;j<n;j++)
	    for(int k=j;k<n;k++){ 
			
			Real T_j=factorLoading->getT(j),
			     T_k=factorLoading->getT(k);
			B(j,k)=factorLoading->getVolSurface()->integral_sgsg(0.0,T_i,T_j,T_k);
		}
		
		// the x_j=c_jG_j, j>i, at time t=0, book, 6.11.3
	    RealVector cG(n,i);                                    // so it exists for i=n-1, although superfluous
	    for(int j=i+1;j<n;j++) cG[j]=-c[j]*U0[j]/H0[i+1];
		
		// coefficients u,v,w of quadratic equation for k_i, book 6.11.3
		Real uu,vv,ww, diag;               
		uu=B(i,i);
		vv=0.0;
		for(int k=i+1;k<n;k++) vv+=cG[k]*rho(i,k)*B(i,k);
			
		ww=0.0; diag=0.0;
		// diagonal + 2*upper half
		for(int j=i+1;j<n;j++) diag+=cG[j]*cG[j]*B(j,j);
		for(int j=i+1;j<n;j++)
		for(int k=j+1;k<n;k++) ww+=cG[j]*cG[k]*rho(j,k)*B(j,k);
			
		ww*=2.0; ww+=diag;
		
		Real  SigmaSquare=capletImpliedSigma[i]*capletImpliedSigma[i],
		      q=vv*vv-uu*(ww-SigmaSquare);
		c[i]=(-vv+sqrt(q))/uu;
//cerr << "\n\nSigmaSquare = " << SigmaSquare << ", c["<<i<<"] = " << c[i];
		if(c[i]<0.0){
			
			cout << "\n\nDriftlessLmmCalibrator::setScalingFactors():"
			     << "\nCaplet prices force negative scaling factor c["<<i<<"]. Terminating.";
			exit(1);
		}
	} // end for i
} // setScalingFactors
		

// READ-WRITE TEST

void
DriftlessLmmCalibrator::
testIO()
{
	int n=50;
	const char *CapletsIn="CapletsIn-DL-dim50-JR-CS.txt",
	           *SwaptionsIn="SwaptionsIn-DL-dim50-JR-CS.txt";
		
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(n,VolSurface::M,Correlations::JR);
		
	LmmCalibrator* cal=new DriftlessLmmCalibrator(fl,CapletsIn,SwaptionsIn);
// cerr << "\n\n0";
	cal->readCaplets();
// cerr << "\n1";
	cal->readSwaptions();
// cerr << "\n2";
	cal->write<CapletData>(std::cout,cal->getCaplets()); 
	cal->writeCaplets();
// cerr << "\n3";
	cal->write<SwaptionData>(std::cout,cal->getSwaptions()); 
	cal->writeSwaptions();
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


// CAPLET AND SWAPTION FORWARD PRICES

Real 
PredictorCorrectorLmmCalibrator::
capletForwardPrice(int i, Real strike)
{
    Real T_i=factorLoading->getT(i), Sigma;

	// Caplet aggregate volatility
	Sigma=sqrt(factorLoading->integral_sgi_sgj_rhoij(i,i,0,T_i));
	return LmmCalibrator::capletForwardPrice(i,strike,Sigma);
}


Real 
PredictorCorrectorLmmCalibrator::
swaptionForwardPrice(int i, Real strike)
{
	// get the matrix R for both Caplet(i) and Swaption(i,n)
	writeCovariationMatrix(i,n);
	
	// Swaption(i,n), vector x, book 6.7
	for(int j=i+1;j<n;j++) x[j]=(B_i0(j)-B_i0(j+1))/B_pq0(i,n);
	 		 
	// swaption volatility = sqrt(xCx')
    Real Sigma=root_xCx(i)/swapRate(i,n);
	return LmmCalibrator::swaptionForwardPrice(i,strike,Sigma);
}



void 
PredictorCorrectorLmmCalibrator::
setScalingFactors()
{	
    RealArray1D& c=factorLoading->getScalingFactors();
	for(int i=1;i<n;i++){	
		
		Real T_i=factorLoading->getT(i);
        VolSurface* volSurface=factorLoading->getVolSurface();
		Real int_sgsg=volSurface->integral_sgsg(0.0,T_i,T_i,T_i);
	
		// impliedSigma=k_i^2*int_sgsg
		c[i]=capletImpliedSigma[i]/sqrt(int_sgsg);
	} // end for i
} // writeInitialScalingFactors


	
void 
PredictorCorrectorLmmCalibrator::
testCalibration(int nVals)
{
	int n=50;
	const char *CapletsIn="CapletsIn-DL-dim50-JR-CS.txt",
	           *SwaptionsIn="SwaptionsIn-DL-dim50-JR-CS.txt";
		
	LiborFactorLoading* 
	fl=LiborFactorLoading::sample(n,VolSurface::M,Correlations::JR);
		
	LmmCalibrator* cal=new PredictorCorrectorLmmCalibrator(fl,CapletsIn,SwaptionsIn);
	cal->calibrate(nVals);
}

 

/**********************************************************************************
 *
 *               LmmOptimizer
 *
 *********************************************************************************/


bool 
SobolLiborCalibrationOptimizer::
isInDomain(const RealArray1D& x) const
{
	Real a=x[0], b=x[1], c=x[2], d=x[3],
	     alpha=x[4], beta=x[5], r_oo=x[6],
	     sum;

	bool is_In=true;
	
	// volsurface a,b,c,d
	switch(cal->getFactorLoading()->getType()->volType){
			
		case VolSurface::M  :   is_In=is_In&&((a>0)&&(d>0.05));  break; // book, 6.11.1 
		case VolSurface::JR :   is_In=is_In&&(d*c+b>0.0); 
	}
	
	// correlations alpha, beta, rho
	switch(cal->getFactorLoading()->getType()->corrType){
			
		case Correlations::CS : sum=alpha/6+beta/3+log(r_oo);
		is_In=is_In&&((alpha>5*beta)&&(beta>0.0)&&(r_oo>0.05)&&(r_oo<0.99)&&(sum<0.0)); break;
		case Correlations::JR : is_In=is_In&&((beta>0.01)&&(beta<0.9)); 
	}
	return is_In;
}
	

Real 
SobolLiborCalibrationOptimizer::
f(const RealArray1D& x){ return cal->objectiveFunction(x); }



			

MTGL_END_NAMESPACE(Martingale)				
