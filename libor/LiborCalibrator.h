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
 *               LmmOptimizer
 *
 *********************************************************************************/
 

// forward declaration
class LmmCalibrator;


/** SobolSearch optimizer to calibrate a {@link LiborFactorLoading}.
 */
class SobolLiborCalibrationOptimizer : public SobolSearch {
	
protected:
	
	LmmCalibrator* cal;
	
public:
	
	/** n+7 is the number of parameters for a factor loading.
	 *
	 *  @param n dimension (number of variables).
	 *  @param cal the LMM calibrator.
	 */
	SobolLiborCalibrationOptimizer
	(int n, LmmCalibrator* _cal, const RealArray1D& x0, int nVals, const RealArray1D& delta) : 
	SobolSearch(n,x0,nVals,delta,true),
	cal(_cal)
	{  }
	
	bool isInDomain(Real* x);
	
	Real f(Real* x);
	
}; // end SobolLiborCalibrationOptimizer
				




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
 *  We don't write the fields calibratedForwardPrice and error
 *  since we wan to use this function to write synthetic
 *  calibration data.
 */
std::ostream& operator << (std::ostream& os, const SwaptionData& swpn);


/** Read swaption data from stream: p, q, strike, forward price. */
std::istream& operator >> (std::istream& is, SwaptionData& swpn);




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
std::ostream& operator << (std::ostream& os, const CapletData& cplt);

/** Read caplet data from stream: i, strike, forward price.
 */
std::istream& operator >> (std::istream& is, CapletData& cplt);



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
	int n;                                   // dimension of Libor process
	// from the factorloading:
	const RealArray1D& x;                    // x[j]=X_j(0)
	const RealArray1D& delta;                // delta[j]=delta_j, accrual periods
	
	
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
    {  
		
		std::cout << "\n\n\nIntitial term structure: fl" << endl
		          << fl->getInitialXLibors()
		          << "\n\n\nIntitial term structure: calibrator" << endl
		          << x
		          << "\n\n\nDeltas: fl" << endl
		          << fl->getDeltas()
		          << "\n\n\nDeltas: calibrator" << endl
		          << delta;
		
	}
	
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
			
			is >> *instrument;
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
			os << *currentInstrument;
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
     Real H0();
	 
	 /** H_i(0) */
     Real H_i0(int i);
     
     /** B_i(0) */
     Real B0(int i);

     /** S_{pq}(0) */
     Real swapRate(int p, int q);
     
     /** B_{pq}(0) */
     Real B_pq(int p, int q);
	      
     /** H_pq(0) */
     Real H_pq(int p, int q);

	
// CAPLET AND SWAPTION PRICES
	
   /** Volatility to expiration \f$T_i\f$.
	*/
   virtual Real capletAggregateVolatility(int i) = 0;
	
	
   /** Swaption on \f$[T_p,T_q]\f$, volatility to expiration \f$T_p\f$.
	*/
   virtual Real swaptionAggregateVolatility(int p, int q) = 0;


   /** Caplet on\f$[T_i,T_{i+1}]\f$.
	*/
   Real capletForwardPrice(int i, Real strike);
   
   /** Swaption on\f$[T_i,T_{i+1}]\f$, expires at \f$T_p\f$.
	*/
   Real swaptionForwardPrice(int p, int q, Real strike);
   
   
// WRITING SYNTHETIC DATA
   
   /** Computes the analytic forward prices  of all at the money caplets and swaptions
    *  at least two periods in the future and writes them to the caplet and swaption
    *  outfiles in the format in which the calibrator reads them for calibration. 
    *  Note that these data cannot be read in the same pass of the 
    *  program since the outfiles are associated with the ofstreams.
    */  
   void writeSyntheticData();
   
   
// OBJECTIVE FUNCTION
	
	
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all caplets or swaptions.
	*  Parameters: the parameters of the current factorloading.
	*  Also writes the calibrated forward prices and relative errors into each caplet
    *  and swaption.
    */
   Real objectiveFunction(){ return squaredCapletError()+squaredSwaptionError(); }
   
   			
   
	
// CALIBRATION 

   
   /** Mean relative error over all caplets and swaptions.
    *  Assumes the calibrated forward prices have already been written 
    *  into the caplets and swaptions.
    */
   Real meanRelativeCalibrationError();
   
   
   	
   /** The objective function as a function of the parameter values x
	*  applied to the factor loading.
    */
   Real objectiveFunction(Real* x);
	
     
   /** Calibrates the factor loading with nVals evaluations of the objective function
	*  (sum squared error over all caplets and swaptions which have been read) and prints
	*  the result to the files "CapletsOut.txt", "SwaptionsOut.txt".
    *
	*  @return the means squared error.
	*/
   Real calibrate(int nVals);	

   
   
private:
   
      
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all
    *  caplets. Calibrated price from the current parameters of the
    *  factor loading. Writes the calibrated forward price and
    *  relativec error into each caplet.
    */
   Real squaredCapletError();
	
	
   /** Sum of (forwardPrice-calibratedForwardPrice)^2 over all
    *  swaptions. Calibrated price from the current parameters of the
    *  factor loading. Writes the calibrated forward price and
    *  relativec error into each caplet.
    */
   Real squaredSwaptionError();
	
	
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
	

     Real capletAggregateVolatility(int i);
	 
	 

     Real swaptionAggregateVolatility(int p, int q);


// SYNTHETIC DATA GENERATION
	 
	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n to files in the directory SyntheticData in the src
	  *  directory.
	  *  
	  * @param n dimension of Libor process (number of accrual intervals).
	  * @param volType type of {@link VolatilitySurface} VolSurface::JR,M,CONST.
	  * @param corrType type of {@link Correlations} Correlations::JR,CS.
	  */ 
	 static void writeSyntheticDataSample
	 (int n, int volType=VolSurface::JR, int corrType=Correlations::CS);
		 

	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n=20,30,40 for all factorloading types.
	  */ 
	 static void writeSyntheticDataSample();
	 
	
	/** Test in dimension 50. Calibrates DriftlessLMM to synthetic data 
	 *  produced by PredictorCorrectorLMM.
	 *
	 * @param nVals number of evaluations of the objective function.
	 */
	static void testCalibration(int nVals);
	

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
	

     Real capletAggregateVolatility(int i);
	 
	 

     Real swaptionAggregateVolatility(int p, int q);

	 
// SYNTHETIC DATA GENERATION
	 
	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n to files in the directory SyntheticData in the src
	  *  directory.
	  *  
	  * @param n dimension of Libor process (number of accrual intervals).
	  * @param volType type of {@link VolatilitySurface} VolSurface::JR,M,CONST.
	  * @param corrType type of {@link Correlations} Correlations::JR,CS.
	  */ 
	 static void writeSyntheticDataSample
	 (int n, int volType=VolSurface::JR, int corrType=Correlations::CS);
		 

	 /** Writes a sample of synthetic caplet and swaption prices in 
	  *  dimension n=20,30,40 for all factorloading types.
	  */ 
	 static void writeSyntheticDataSample();
	 
}; // end PredictorCorrectorLmmCalibrator



			



MTGL_END_NAMESPACE(Martingale)

#endif
