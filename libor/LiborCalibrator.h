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

#ifndef martingale_lmm_calibrator_h    
#define martingale_lmm_calibrator_h

#include "FinMath.h"
#include "LiborFactorLoading.h"
#include "Optimizer.h"
#include <iostream>
#include <fstream>
#include <list>
#include <string>

//#include <string>
#include <math.h>


MTGL_BEGIN_NAMESPACE(Martingale)



/**********************************************************************************
 *
 *               Swaptions and Caplets
 *
 *********************************************************************************/
 
/** Accrual interval \f$[T_p,T_q]\f$ (given by p,q), strike and analytic 
 *  forward swaption price.
 */
struct SwaptionData {
	
	int p,q;
	Real strike, forwardPrice, calibratedForwardPrice, error;
	
	/** All fields intialized with 0. */
	SwaptionData() :
	p(0), q(0), strike(0.0), forwardPrice(0.0), 
	calibratedForwardPrice(0.0), error(0.0)
    {   }
	
};



/** Write swaption data to stream: p, q, strike, forward price.
 */
std::ostream& operator << (std::ostream& os, SwaptionData* swpn)
{
	os << swpn->p << "  " << swpn->q << "  " << swpn->strike << "  " 
	   << swpn->forwardPrice;
    return os;
} // end operator <<


/** Read swaption data from stream: p, q, strike, forward price. */
std::istream& operator >> (std::istream& is, SwaptionData* swpn)
{
	is >> swpn->p >> swpn->q >> swpn->strike >> swpn->forwardPrice;
    return is;
} // end operator <<





/** Accrual interval \f$[T_p,T_{p+1}]\f$ (given by p), strike and analytic 
 *  forward caplet price.
 */
struct CapletData {
	
	int i;
	Real strike, forwardPrice, calibratedForwardPrice, error;
	
	/** All fields intialized with 0. */
	CapletData() :
	i(0), strike(0.0), forwardPrice(0.0), 
	calibratedForwardPrice(0.0), error(0.0)
    {   }
	
};
 
 

/** Write caplet data to stream: i, strike, forward price.
 */
std::ostream& operator << (std::ostream& os, CapletData* cplt)
{
	os << cplt->i << "  " << cplt->strike << "  " << cplt->forwardPrice << "  ";
    return os;
} // end operator <<


/** Read caplet data from stream: i, strike, forward price.
 */
std::istream& operator >> (std::istream& is, CapletData* cplt)
{
	is >> cplt->i >> cplt->strike >> cplt->forwardPrice;
    return is;
} // end operator <<




/**********************************************************************************
 *
 *               LiborCalibrator
 *
 *********************************************************************************/
	

/** Class which calibrates the factorloading of any given LMM to a set of
 *  caplet and swaption prices read from two separate files respectively.
 */

class LmmCalibrator {
	
	
protected:

	// reading data
	ifstream capletsIn;
	ifstream swaptionsIn;
	
	// writing synthetic data, calibration results
	ofstream capletsOut;	
	ofstream swaptionsOut;
	
	std::list<CapletData*> caplets;
	std::list<SwaptionData*> swaptions;
	
	LiborFactorLoading* factorLoading;
	int n;                                  // dimension of Libor process
	Real* x;                                // x[j]=X_j(0)
	Real* delta;                            // delta[j]=delta_j, accrual periods
	
	
public:
	
// CONSTRUCTOR
	
	LmmCalibrator
	(LiborFactorLoading* fl,
	 const char* capletsInFile="CapletsIn.txt",  
	 const char* swaptionsInFile="SwaptionsIn.txt",
	 const char* capletsOutFile="CapletsOut.txt",
	 const char* swaptionsOutFile="SwaptionsOut.txt"
	 ) :
	capletsIn(capletsInFile),
	swaptionsIn(swaptionsInFile),
	capletsOut(capletsOutFile),
	swaptionsOut(swaptionsOutFile),
	caplets(),
	swaptions(),
	factorLoading(fl),
	n(fl->getDimension()),
	x(fl->getInitialXLibors()),
	delta(fl->getDeltas())
    {  }
	
// ACCESSORS
	
	LiborFactorLoading* getFactorLoading(){ return factorLoading; }
	
	
// READING AND WRITING CAPLETS AND SWAPTIONS
// Instrument = CapletData, SwaptionData  

	/** Read Instruments from ifstream if and append to list dataList
	 */
	template<typename Instrument>
	void read(istream& is, std::list<Instrument*>& dataList)
    {
		Instrument* instrument = new Instrument();
		while(!is.eof()){  
			
			is >> instrument;
			dataList.push_back(instrument);
			instrument = new Instrument();
		}
		
		if(!is.eof()){ 
			
			cout << "\n\nLiborCalibrator.readCaplets(): "
			     << "failed to read all instruments. Terminating.";
			exit(EXIT_FAILURE);
		}
	    // the last caplet was allocated but nothing left to be read
		// delete it
		if(instrument) delete instrument;
			
    } // end read

	
    /** Write Instruments from list dataList to ofstream.
	 *  Writes all fields including calibratedForwardPrice and relative error.
	 */
	template<typename Instrument>
	void write(ostream& os, std::list<Instrument*> dataList)
    {
		// iterator is pointer to pointer to Instrument
		std::list<Instrument*>::const_iterator theInstrument;
		for(theInstrument=dataList.begin(); theInstrument!=dataList.end(); ++theInstrument)			
		{
			Instrument* currentInstrument=*theInstrument;
			os << currentInstrument;
			os << "  " << currentInstrument->calibratedForwardPrice
			   << "  " << currentInstrument->error << endl;
		}
		
		if(theInstrument!=dataList.end()){ 
			
			cout << "\n\nLiborCalibrator.write(): "
			     << "failed to write all instruments. Terminating.";
			exit(EXIT_FAILURE);
		}
    } // end write
	
	
	void readCaplets(){ read<CapletData>(capletsIn,caplets); }
	void readSwaptions(){ read<SwaptionData>(swaptionsIn,swaptions); }
	
	
	/** Writes caplets from list. 
	 */
	void writeCaplets()
	{ 
		capletsOut << "i   strike    forwardPrice    calibratedForwardPrice    error"
		           << endl << endl;
		write<CapletData>(capletsOut,caplets); 
	}
	
	
	/** Writes swaptions from list. 
	 */
	void writeSwaptions()
	{ 
		capletsOut << "p   q   strike    forwardPrice    calibratedForwardPrice    error"
		           << endl << endl;
		write<SwaptionData>(swaptionsOut,swaptions); 
	}

	
// ACCRUAL FACTORS, BONDS, SWAP RATES, .. AT TIME ZERO
	
	
	 /** H_0(0) */
     Real H0()
     {
         Real f=1.0;
	     for(int k=0;k<n;k++)f*=(1+x[k]);
	     return f;
     }
	 
	 /** H_i(0) */
     Real H_i0(int i)
     {
         Real f=1.0;
	     if(i==n) return f;
		
	     for(int k=i;k<n;k++)f*=(1+x[k]);
	     return f;
     }
     

     /** B_i(0) */
     Real B0(int i)
     { 
         Real f=1.0;
         // accumulate 1 from time t=0 to time t=T_i                    
         for(int j=0;j<i;j++)f*=(1.0+x[j]); 
         return 1.0/f;                                 // B_i(0)=1/f  
     }   
 


     /** S_{pq}(0) */
     Real swapRate(int p, int q)
     { 
        Real f=1.0+x[q-1], S=delta[q-1];
        for(int k=q-2;k>=p;k--){ S+=delta[k]*f; f*=(1.0+x[k]); }
 
        return (f-1.0)/S;
     } 


     
     /** B_{pq}(0) */
     Real B_pq(int p, int q)
     { 
         Real S=0.0, F=B0(q);
         for(int k=q-1;k>=p;k--){ S+=delta[k]*F; F*=(1.0+x[k]); }
         return S;
     } 
	 

     
     /** H_pq(0) */
     Real H_pq(int p, int q)
     { 
         return B_pq(p,q)*H0();
     } 


	
// CAPLET AND SWAPTION PRICES
	
   /** Volatility to expiration \f$T_i\f$.
	*/
   virtual Real capletAggregateVolatility(int i) = 0;
	
	
   /** Swaption on \f$[T_p,T_q]\f$, volatility to expiration \f$T_p\f$.
	*/
   virtual Real swaptionAggregateVolatility(int p, int q) = 0;


   /** Caplet on\f$[T_i,T_{i+1}]\f$.
	*/
   Real capletForwardPrice(int i, Real strike)
   {
	   Real delta_i=delta[i], 
	        Li0=factorLoading->getInitialTermStructure()[i],            // L_i(0)
            cSigma=capletAggregateVolatility(i),                        // aggregate volatility to expiry    
	        Nplus=FinMath::N(FinMath::d_plus(Li0,strike,cSigma)),
	        Nminus=FinMath::N(FinMath::d_minus(Li0,strike,cSigma)),
	        f=H_i0(i+1);                                                 
 
       return delta_i*(Li0*Nplus-strike*Nminus)*f; 
   }

   
   /** Swaption on\f$[T_i,T_{i+1}]\f$, expires at \f$T_p\f$.
	*/
   Real swaptionForwardPrice(int p, int q, Real strike) 
   {
       Real S_pq=swapRate(p,q),                                         // S_{p,q}(0)
            swpnSigma=swaptionAggregateVolatility(p,q),                 // aggregate volatility to T_p    
	        Nplus=FinMath::N(FinMath::d_plus(S_pq,strike,swpnSigma)),
	        Nminus=FinMath::N(FinMath::d_minus(S_pq,strike,swpnSigma)),
	        f=H_pq(p,q);                                                // forward B_{p,q}(0)
  
       return f*(S_pq*Nplus-strike*Nminus);  
   } 
   
   
// WRITING SYNTHETIC DATA
   
   /** Computes the analytic forward prices  of all at the money caplets and swaptions
    *  at least two periods in the future and writes them to the caplet and swaption
    *  outfiles in the format in which the calibrator reads them for calibration. 
    *  Note that these data cannot be read in the same pass of the 
    *  program since the outfiles are associated with the ofstreams.
    */  
   void writeSyntheticData()
   {
	   // write caplets
	   for(int i=2;i<n;i++){
		   
		   Real strike = factorLoading->getInitialTermStructure()[i];  // L_i(0)
		   Real forwardPrice = capletForwardPrice(i,strike);
		   
		   CapletData* caplet = new CapletData();
		   caplet->i=i;
		   caplet->strike=strike;
		   caplet->forwardPrice=forwardPrice;
		   
		   capletsOut << caplet << endl;
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
		   
		   swaptionsOut << swaption << endl;
	   } 
   } // end writeSyntheticData
   
   
// OBJECTIVE FUNCTION
	
	
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all caplets or swaptions.
	*  Parameters: the parameters of the current factorloading.
	*  Also writes the calibrated forward prices and relative errors into each caplet
    *  and swaption.
    */
   Real objectiveFunction(){ return squaredCapletError()+squaredSwaptionError(); }
   
   			
   
	
// CALIBRATION 
// Is done in the factorloadings which know their own parameters.
	

   /** Calibrates and returns pointer to calibrate factorloading.
    */
   LiborFactorLoading* calibrate(){ return factorLoading->calibrateSelf<LmmCalibrator>(this); }
   
   
   /** Mean relative error over all caplets and swaptions.
    *  Assumes the calibrated forward prices have already been written 
    *  into the caplets and swaptions.
    */
   Real meanRelativeCalibrationError()
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
   
   
private:
   
      
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all
    *  caplets. Calibrated price from the current parameters of the
    *  factor loading. Writes the calibrated forward price and
    *  relativec error into each caplet.
    */
   Real squaredCapletError()
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
	
	
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all
    *  swaptions. Calibrated price from the current parameters of the
    *  factor loading. Writes the calibrated forward price and
    *  relativec error into each caplet.
    */
   Real squaredSwaptionError()
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


		   
}; // end LiborCalibrator




/**********************************************************************************
 *
 *              DriftLessLMM-Calibrator
 *
 *********************************************************************************/


/** Calibrator for a {@link DriftlessLMM}.
 */
class DriftlessLmmCalibrator : public LmmCalibrator {
	
	
public:
	
	
	DriftlessLmmCalibrator
	(LiborFactorLoading* fl,
	 const char* capletsInFile="CapletsIn.txt", 
	 const char* swaptionsInFile="SwaptionsIn.txt", 
	 const char* capletsOutFile="CapletsOut.txt", 	 
	 const char* swaptionsOutFile="SwaptionsOut.txt"
	) :	
    LmmCalibrator(fl,capletsInFile,swaptionsInFile,capletsOutFile,swaptionsOutFile) 
    {    }
	
	
// CAPLET AND SWAPTION AGGREGATE VOLS
	

     Real capletAggregateVolatility(int i)
     { 
		 Real Ti=factorLoading->getTenorStructure()[i];
		 UTRMatrix<Real>& R=factorLoading->logLiborCovariationMatrix(i,n,0,Ti).utrRoot();
		 vector<Real> x(n-i,i);
		 x[i]=1;
		 Real f=H_i0(i+1);
		 for(int j=i+1;j<n;j++){ Real Uj0=H_i0(j)-H_i0(j+1); x[j]=-Uj0/f; }
		 x*=R.transpose();
		 return x.norm();
     } // end capletAggregateVolatility
	 
	 

     Real swaptionAggregateVolatility(int p, int q)
     { 
         Real* T=factorLoading->getTenorStructure();
		 UTRMatrix<Real>& Q=factorLoading->logLiborCovariationMatrix(p,n,0,T[p]).utrRoot();
		 
		 vector<Real> x(n-p,p);
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
	 
	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n to files with the following nomenclature:
	  *  <br>CapletsIn-dim(n)-DL-JR.txt<br>
	  *  is a set of caplet prices in dimension n based on the sample 
	  *  {@link JR_FactorLoading} and a {@link DriftlessLMM}.
	  *  The files are written to the directory SyntheticData in the src
	  *  directory.
	  *  
	  * @param n dimension of Libor process (number of accrual intervals).
	  * @param flType type of factorloading: LiborFactorLoading::CS,JR,CV.
	  */ 
	 static void writeSyntheticDataSample(int n, int flType)
     {
		 Int_t nn(n);   // for string conversion
		 string capletOutFile, swaptionOutFile, dim=nn.toString();
		 LiborFactorLoading* fl;
		 switch(flType){
			 
			 case LiborFactorLoading::JR : 
				  fl=JR_FactorLoading::sample(n);
				  capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-DL-JR.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-DL-JR.txt"; 
			      break;
			 case LiborFactorLoading::CV :
				  fl=ConstVolLiborFactorLoading::sample(n);
				  capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-DL-CV.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-DL-CV.txt"; 
			      break;
			 default : // Coffee-Shoenmakers
			      fl=CS_FactorLoading::sample(n);
			      capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-DL-CS.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-DL-CS.txt"; 
		 } // end switch
		 
		 LmmCalibrator* cal=new DriftlessLmmCalibrator
		 (fl,"CapletsIn","SwaptionsIn",capletOutFile.c_str(),swaptionOutFile.c_str());
		 
		 cal->writeSyntheticData();
		
	 } // end writeSyntheticDataSample
		 

	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n=20,30,40 for all three factorloading types.
	  */ 
	 static void writeSyntheticDataSample()
     {
	     for(int n=20;n<60;n+=10)
	     for(int flType=0;flType<3;flType++)writeSyntheticDataSample(n,flType);
	     cout << "\n\nDone.";
	 }
			 
			 
			 
	/** Allocates a sample factorloading of dimension n and
	 *  type flType (LiborFactorLoading::CS,JR,CV) reads 
	 *  caplet and swaption prices from the infiles and writes them
	 *  to the outfiles.
	 */
	static void test(int n, int flType)
    {
		LiborFactorLoading* fl;
		switch(flType)
	    {
			case LiborFactorLoading::JR : fl = JR_FactorLoading::sample(n); 
				                          break;
			case LiborFactorLoading::CV : fl = ConstVolLiborFactorLoading::sample(n); 
				                          break;
			default                     : fl = CS_FactorLoading::sample(n);
		}
		
		LmmCalibrator* cal=new DriftlessLmmCalibrator(fl);
		cal->readCaplets();
		cal->writeCaplets();
		cal->readSwaptions();
		cal->writeSwaptions();
	}
	

	
}; // end DriftlessLmmCalibrator




/**********************************************************************************
 *
 *              PredictorCorrectorLMM-Calibrator
 *
 *********************************************************************************/



/** Calibrator for a {@link PredictorCorrectorLMM}.
 */
class PredictorCorrectorLmmCalibrator : public LmmCalibrator {
	
	
public:
	
	
	PredictorCorrectorLmmCalibrator
	(LiborFactorLoading* fl,
	 const char* capletsInFile="CapletsIn.txt", 
	 const char* swaptionsInFile="SwaptionsIn.txt", 
	 const char* capletsOutFile="CapletsOut.txt", 	 
	 const char* swaptionsOutFile="SwaptionsOut.txt"
	) :	
    LmmCalibrator(fl,capletsInFile,swaptionsInFile,capletsOutFile,swaptionsOutFile) 
    {    }
	
	
// CAPLET AND SWAPTION AGGREGATE VOLS
	

     Real capletAggregateVolatility(int i)
     { 
         Real* T=factorLoading->getTenorStructure();
		 Real volsqr=factorLoading->integral_sgi_sgj_rhoij(i,i,0,T[i]);
		 return sqrt(volsqr);
     } 
	 
	 

     Real swaptionAggregateVolatility(int p, int q)
     { 
          Real* T=factorLoading->getTenorStructure();
		  UTRMatrix<Real>& 
		  Q=factorLoading->logLiborCovariationMatrix(p,q,0,T[p]).utrRoot();
		  vector<Real> x_pq(q-p,p);
		  for(int j=p;j<q;j++) x_pq[j]=(B0(j)-B0(j+1))/B_pq(p,q);
		  x_pq*=Q.transpose();
		  return x_pq.norm();
     } 

// SYNTHETIC DATA GENERATION
	 
	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n to files with the following nomenclature:
	  *  <br>CapletsIn-dim(n)-DL-JR.txt<br>
	  *  is a set of caplet prices in dimension n based on the sample 
	  *  {@link JR_FactorLoading} and a {@link DriftlessLMM}.
	  *  The files are written to the directory SyntheticData in the src
	  *  directory.
	  *  
	  * @param n dimension of Libor process (number of accrual intervals).
	  * @param flType type of factorloading: LiborFactorLoading::CS,JR,CV.
	  */ 
	 static void writeSyntheticDataSample(int n, int flType)
     {
		 Int_t nn(n);   // for string conversion
		 string capletOutFile, swaptionOutFile, dim=nn.toString();
		 LiborFactorLoading* fl;
		 switch(flType){
			 
			 case LiborFactorLoading::JR : 
				  fl=JR_FactorLoading::sample(n);
				  capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-PC-JR.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-PC-JR.txt"; 
			      break;
			 case LiborFactorLoading::CV :
				  fl=ConstVolLiborFactorLoading::sample(n);
				  capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-PC-CV.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-PC-CV.txt"; 
			      break;
			 default : // Coffee-Shoenmakers
			      fl=CS_FactorLoading::sample(n);
			      capletOutFile="SyntheticData/CapletsIn-dim"+dim+"-PC-CS.txt";
			      swaptionOutFile="SyntheticData/SwaptionsIn-dim"+dim+"-PC-CS.txt"; 
		 } // end switch
		 
		 LmmCalibrator* cal=new PredictorCorrectorLmmCalibrator
		 (fl,"CapletsIn","SwaptionsIn",capletOutFile.c_str(),swaptionOutFile.c_str());
		 
		 cal->writeSyntheticData();
		
	 } // end writeSyntheticDataSample
		 

	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n=20,30,40 for all three factorloading types.
	  */ 
	 static void writeSyntheticDataSample()
     {
	     for(int n=20;n<60;n+=10)
	     for(int flType=0;flType<3;flType++)writeSyntheticDataSample(n,flType);
	     cout << "\n\nDone.";
	 }
			 
			 
			 
	/** Allocates a sample factorloading of dimension n and
	 *  type flType (LiborFactorLoading::CS,JR,CV) reads 
	 *  caplet and swaption prices from the infiles and writes them
	 *  to the outfiles.
	 */
	static void test(int n, int flType)
    {
		LiborFactorLoading* fl;
		switch(flType)
	    {
			case LiborFactorLoading::JR : fl = JR_FactorLoading::sample(n); 
				                          break;
			case LiborFactorLoading::CV : fl = ConstVolLiborFactorLoading::sample(n); 
				                          break;
			default                     : fl = CS_FactorLoading::sample(n);
		}
		
		LmmCalibrator* cal=new PredictorCorrectorLmmCalibrator(fl);
		cal->readCaplets();
		cal->writeCaplets();
		cal->readSwaptions();
		cal->writeSwaptions();
	}
	
	
}; // end PredictorCorrectorLmmCalibrator



/**********************************************************************************
 *
 *               LmmOptimizer
 *
 *********************************************************************************/



/** Optimizer to calibrate a {@link CV_FactorLoading}.
 *  We use a BFGS optimizer.
 */
class CV_Optimizer : public BFGS {
	
protected:
	
	LmmCalibrator* cal;
	LiborFactorLoading* factorLoading;    // the factorloading calibrated by cal
	
public:
	
	/** @param n dimension (number of variables).
	 *  @param cal the LMM calibrator.
	 */
	LmmOptimizer(int n, LmmCalibrator* _cal) : Optimizer(n) 
	cal(_cal),
	factorLoading(cal->getFactorLoading())
	{  }
	
	
};


/** Optimizer for 
	
	





MTGL_END_NAMESPACE(Martingale)

#endif
