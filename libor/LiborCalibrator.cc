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
         f=H_pq(i,n);                                                // forward B_{i,n}(0)
  
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

        CapletData* currentCaplet=caplets[i];
		sum+=abs(currentCaplet->error);
	}
		
	return sum/(2*(n-1));
} // end meanRelativeCalibrationError
   	

Real 
LmmCalibrator::
objectiveFunction(Real* x)
{
// cerr << "\n\n\nparameter vector for objective function:" << endl;
// for(int i=n;i<n+7;i++) cerr << "x[i] = " << x[i] << endl;
	factorLoading->setParameters(x);		
	return objectiveFunction();
}
	


Real 
LmmCalibrator::
calibrate(int nVals)
{
	readCaplets();
	//computeCapletMatrices();
	readSwaptions();
cerr << "\n\nRead data";
	//computeSwaptionMatrices();
	// initial guess for parameter values
	RealArray1D u(n+7);
			
    // initial guess for caplet vol scaling factors k[i]
	// iterator is pointer to pointer to CapletData
    for(int i=1;i<n;i++){				
		
		CapletData* currentCaplet=caplets[i];
		Real T_i=factorLoading->getT(i),
		     L_i0=x[i]/delta[i],
		     forwardPrice=currentCaplet->forwardPrice,
		     strike=currentCaplet->strike;
		
		// the blackImpliedAggregateCallVolatility solves the blackScholesFunction.
		// the capletForwardPrice multiplies this function with delta_i*H_i0(i+1)
		// we must devide by this to isolate the blackScholesFunction
		Real price=forwardPrice/(delta[i]*H_i0(i+1));

		Real vol=FinMath::blackImpliedAggregateCallVolatility(L_i0,strike,price);
			
		VolSurface* volSurface=factorLoading->getVolSurface();
		Real int_sgsg=volSurface->integral_sgsg(0.0,T_i,T_i,T_i);
	
		u[i]=vol/sqrt(int_sgsg);
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
		u[n+4]=0.2; u[n+5]=1.0; u[n+6]=1.0; break;
			
		case Correlations::CS :
		u[n+4]=1.5; u[n+5]=0.05; u[n+6]=0.2; break;
			
		default : u[n+4]=1.0; u[n+5]=1.0; u[n+6]=1.0; 
	}	
	
// print initialization parameters
for(int i=1;i<n;i++) cerr << "\nk["<<i<<"] = " << u[i];
cerr << "\na = " << u[n]
	 << "\nb = " << u[n+1]	
	 << "\nc = " << u[n+2]
     << "\nd = " << u[n+3]
	 << "\nalpha = " << u[n+4]	
	 << "\nbeta = " << u[n+5]
	 << "\nr_oo = " << u[n+6];

	RealArray1D du(n+7);
	for(int i=0;i<n+7;i++) du[i]=0.01;
		 		 
	SobolLiborCalibrationOptimizer optimizer(n+7,this,u,nVals,du);
	const RealArray1D& w=optimizer.search();
	
	writeCaplets();
	writeSwaptions();
	
	Real meanError=meanRelativeCalibrationError();

	// report optimal parmeters
	std::cout << "\n\n\n\nOptimal parameters: " << endl << endl;
	for(int i=1;i<n;i++) std::cout << "\nk["<<i<<"] = " << w[i];
	std::cout << "\n\n\na = " << w[n]
	          << "\nb = " << w[n+1]
	          << "\nc = " << w[n+2]
	          << "\nd = " << w[n+3]
	          << "\nalpha = " << w[n+4]
	          << "\nbeta = " << w[n+5]
	          << "\nr_oo = " << w[n+6]
	          << "\n\n\nMean relative calibration error: "
	          << meanError;
	
    return meanError;
		
} // end calibrate




/**********************************************************************************
 *
 *              StandardLMM-Calibrator
 *
 *********************************************************************************/



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


void 
StandardLmmCalibrator::
writeCovariationMatrixRoot(int p, int q)
{	
    UTRRealMatrix& A=cvMatrix;
	UTRRealMatrix& U=cvRoot;
	
	int m=q-1;                        // (m,m) lower left matrix corner
  
	// computation of Uij by reverse induction starting from i,j=n descending to i=j=p
	// get U_ij from relation row_i(U) dot row_j(U) = A_ij, j>=i, at that point all
	// U_rs with r>=i,s>=j (and not both equal) are already computed.
    for(int i=m;i>=p;i--)
    for(int j=m;j>=i;j--) {                      // j>=i
    
        // r_i(U).r_j(U)-U_ijU_jj
		Real sum=0.0;  
        for(int k=j+1;k<m;k++) sum += U(i,k)*U(j,k);   

   		Real Aij=A(i,j), R=Aij-sum;              // R-U_ijU_jj=A_ij-r_i(U).r_j(U)=0    
        if(j>i) U(i,j)=R/U(j,j);             // U[j][0]=U[j][j-j]=U_jj              
        else if(R>0) U(j,j)=sqrt(R);
        else{ cerr << "\n\nwriteCovariationMatrixRoot(): Matrix not positive definite:" 
			       << "\np = "<< p << ", q = " << q
                   << ", A["<<j<<"]["<<j<<"]-S="<< R 
			       << "\nTerminating.";
               exit(1); 
        } // end else
    } //end for i
} // end writeCovariationMatrixRoot



Real 
StandardLmmCalibrator::
capletForwardPrice(int i, Real strike)
{
	// wasteful to set both vols, but used only in data writing
	writeAggregateVolatilities(i,capletVol,swaptionVol);
	return LmmCalibrator::capletForwardPrice(i,strike,capletVol);
}


Real 
StandardLmmCalibrator::
swaptionForwardPrice(int i, Real strike)
{
	// wasteful to set both vols, but used only in data writing
	writeAggregateVolatilities(i,capletVol,swaptionVol);
	return LmmCalibrator::swaptionForwardPrice(i,strike,swaptionVol);
}



Real
StandardLmmCalibrator::
objectiveFunction()
{ 	
	Real strike, marketPrice, calibratedPrice, error=0.0, diff;
	// loop through Caplet(i) and Swaption(i,n)
	for(int i=1;i<n;i++){ 
		
		// set Sigma = capletVol, SwaptionVol 
		// from current parameters of the factor loading
		writeAggregateVolatilities(i,capletVol,swaptionVol); 
		
		CapletData* theCaplet=caplets[i];
		// read in
		marketPrice=theCaplet->forwardPrice; 
		strike=theCaplet->strike;
		// compute from current parameters of the factor loading
		calibratedPrice=LmmCalibrator::capletForwardPrice(i,strike,capletVol); 
		
		diff=calibratedPrice-marketPrice;
		error+=diff*diff;
		
		// write the information into the caplet
		(theCaplet->calibratedForwardPrice)=calibratedPrice;
		(theCaplet->error)=100*diff/marketPrice;
//		theCaplet->setCalibratedForwardPrice(calibratedPrice);
//		theCaplet->setError(100*diff/marketPrice);
		
		SwaptionData* theSwaption=swaptions[i];
		marketPrice=theSwaption->forwardPrice; 
		strike=theSwaption->strike;
		calibratedPrice=LmmCalibrator::swaptionForwardPrice(i,strike,swaptionVol); 

		diff=calibratedPrice-marketPrice;
		error+=diff*diff;
		
		// write the information into the swaption
		(theSwaption->calibratedForwardPrice)=calibratedPrice;
		(theSwaption->error)=100*diff/marketPrice;
//		theSwaption->setCalibratedForwardPrice(calibratedPrice);
//		theSwaption->setError(100*diff/marketPrice);
	}
	return error;
} // end objectiveFunction()
		
 


// Norm of the product Q'x of the vector x=x(p,n) with the transpose of the 
// matrix Q=cvRoot(p,n). Here x=x[i], Q=Q(i,j) with indices p<=i<=j<n.
Real
StandardLmmCalibrator::	
norm_Qtx(int p)
{
	UTRRealMatrix& Q=cvRoot;
	Real s=0.0;
	for(int i=p;i<n;i++){
		
		Real Qtx_i=0.0;                                   // [Q'x]_i
		for(int j=p;j<=i;j++) Qtx_i+=Q(j,i)*x[j];
		s+=Qtx_i*Qtx_i;
	}
    return sqrt(s);
}




// SYNTHETIC DATA GENERATION
void 
StandardLmmCalibrator::
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
 *              DriftlessLMM-Calibrator
 *
 *********************************************************************************/


// CAPLET AND SWAPTION AGGREGATE VOLS

void 
DriftlessLmmCalibrator::
writeAggregateVolatilities(int p, Real& cpltVol, Real& swpnVol)
{
    const RealArray1D& T=factorLoading->getTenorStructure();
	// get the matrix R for both Caplet(i) and Swaption(i,n)
	writeCovariationMatrix(p,n);
    writeCovariationMatrixRoot(p,n);       // matrix Q indices p<=i<=j<n
	
	// Caplet(p), the vector x (book 6.8.2)
	x[p]=1;
	Real f=H_i0(p+1);
	for(int j=p+1;j<n;j++){ Real Uj0=H_i0(j)-H_i0(j+1); x[j]=-Uj0/f; }
	
    // caplet volatility
	cpltVol=norm_Qtx(p);
		
	// Swaption(p,n), the vector x (book 6.8.3)
	 Real denom1=H_i0(p)-H_i0(n),
	      denom2=0;
	 for(int j=p;j<n;j++) denom2+=delta[j]*H_i0(j+1);
		
	 Real Up0=H_i0(p)-H_i0(p+1); x[p]=Up0/denom1;		 
	 for(int j=p+1;j<n;j++){
			 
		 Real Uj0=H_i0(j)-H_i0(j+1);
	     x[j]=Uj0/denom1-(T[j]-T[p])*Uj0/denom2;
	 }		 

	 // swaption volatility
     swpnVol=norm_Qtx(p);
	 
} // end setAggregateVolatilities
 

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
	
	

/**********************************************************************************
 *
 *              PredictorCorrectorLMM-Calibrator
 *
 *********************************************************************************/



// CAPLET AND SWAPTION AGGREGATE VOLS

void 
PredictorCorrectorLmmCalibrator::
writeAggregateVolatilities(int p, Real& cpltVol, Real& swpnVol)
{
    Real T_p=factorLoading->getT(p);
	// get the matrix R for both Caplet(i) and Swaption(i,n)
	writeCovariationMatrix(p,n);
    writeCovariationMatrixRoot(p,n);       // matrix Q indices p<=i<=j<n
	
	// Caplet(p)
	cpltVol=factorLoading->integral_sgi_sgj_rhoij(p,p,0,T_p);
    cpltVol=sqrt(cpltVol);
	
	// Swaption(p,n)
	for(int j=p+1;j<n;j++) x[j]=(B0(j)-B0(j+1))/B_pq(p,n);
    // swaption volatility
    swpnVol=norm_Qtx(p);

} // end setAggregateVolatilities
 

/**********************************************************************************
 *
 *               LmmOptimizer
 *
 *********************************************************************************/


bool 
SobolLiborCalibrationOptimizer::
isInDomain(Real* x)
{
	Real alpha=0.5, beta=0.5, r_oo=0.5, sum;
	for(int i=0;i<n;i++) if(x[i]<0.02) return false;
			
	switch(cal->getFactorLoading()->getCorrelationType()){
			
		case Correlations::CS : 
				
		    alpha=x[n+4]; beta=x[n+5]; r_oo=x[n+6];     // book, 6.11.1 
		    sum=alpha/6+beta/3+log(r_oo);
		    if((alpha<beta)||(beta<0.0)||(r_oo<0.0)||sum>0.0) return false;
			
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
