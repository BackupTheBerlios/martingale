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
	
	/** 7 parameters for volsurface and correlations.
	 *
	 *  @param cal the LMM calibrator.
	 */
	SobolLiborCalibrationOptimizer
	(LmmCalibrator* _cal, const RealArray1D& x0, int nVals, const RealArray1D& delta) : 
	SobolSearch(7,x0,nVals,delta,true),
	cal(_cal)
	{    }
	
	bool isInDomain(const RealArray1D& x) const;
	
	Real f(const RealArray1D& x);
	
	const RealArray1D& ssearch();
	
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
	

/** <p>Class which calibrates the parameters of the factorloading of any given LMM to a 
 *  set of caplet and swaption (forward) prices read from two separate files respectively.
 *  The factor loading is calibrated to all caplets Caplet(i) (on \f$[T_i,T_{i+1}]\f$) 
 *  and all coterminal swaptions Swaption(i,n) (swap on \f$[T_i,T_n]\f$) expiring at time
 *  \f$T_i\f$ respectively. The strikes can be chosen  arbitrarily.
 *
 * <p>The parameters on which the factor loading can depend must satisfy some restrictions
 * detailed here: {@link LiborFactorLoading}.
 */

class LmmCalibrator {
	
	
protected:
	
	LiborFactorLoading* factorLoading;
	int n;                                   // dimension of Libor process
	// from the factorloading:
	const RealArray1D& X0;                   // X0[j]=X_j(0)
	const RealArray1D& delta;                // delta[j]=delta_j, accrual periods
	 
	// caching some quantities
	RealVector U0;                           // U0[j]=U_j(0)
	RealVector H0;                           // H0[j]=H_j(0) 
	Real Bn0;                                // B_n(0)
	UTRRealMatrix S_pq0;                     // swaprates S_pq(0) 
	
	// reading data
	ifstream capletsIn;
	ifstream swaptionsIn;
	
	// writing synthetic data, calibration results
	ofstream logFile;
	ofstream capletsOut;	
	ofstream swaptionsOut;
	
	Array1D<CapletData*> caplets;
	Array1D<SwaptionData*> swaptions;
		
	
public:
	
// CONSTRUCTOR
	
	LmmCalibrator
	(LiborFactorLoading* fl,
	 const char* capletsInFile="CapletsIn.txt",  const char* swaptionsInFile="SwaptionsIn.txt",
	 const char* capletsOutFile="CapletsOut.txt", const char* swaptionsOutFile="SwaptionsOut.txt");
	
// ACCESSORS
	
	LiborFactorLoading* getFactorLoading(){ return factorLoading; }
	Array1D<CapletData*>& getCaplets(){ return caplets; }
	Array1D<SwaptionData*>& getSwaptions(){ return swaptions; }
	
// READING AND WRITING CAPLETS AND SWAPTIONS
// Instrument = CapletData, SwaptionData  

	/** Read Instruments from ifstream is into Array instruments.
	 */
	template<typename Instrument>
	void read(istream& is, Array1D<Instrument*>& instruments)
    {
		for(int i=1;i<n;i++){
		
			instruments[i] = new Instrument();
			is >> *(instruments[i]);
		} 	
    } // end read

	
    /** Write Instruments from array instruments to ofstream.
	 *  Writes all fields including calibratedForwardPrice and relative error.
	 */
	template<typename Instrument>
	void write(ostream& os, const Array1D<Instrument*>& instruments)
    {
		for(int i=1;i<n;i++){
			
			Instrument* currentInstrument=instruments[i];
			os << *currentInstrument;
			os << "  " << currentInstrument->calibratedForwardPrice
			   << "  " << currentInstrument->error << endl;
		}
    } // end write
	
	
	void readCaplets();
	void readSwaptions();
	
	/** Writes caplets from list. 
	 */
	void writeCaplets();
	
	/** Writes swaptions from list. 
	 */
	void writeSwaptions();
	
	
// CORRELATIONS
	
	/** The current correlations \f$\rho_{ij}\f$. */
	Real rho(int i, int j);
	
// ACCRUAL FACTORS, BONDS, SWAP RATES, .. AT TIME ZERO
	
	/** Libor L_i(0) */
    Real L_i0(int i);
	 
	/** H_i(0) */
    Real H_i0(int i);
     
    /** B_i(0) */
    Real B_i0(int i);
	
	/** H_pq(0)=B_pq(0)/B_n(0) */
	Real H_pq0(int p, int q);
	
	/** H_pq(0)=B_pq(0)/B_n(0) */
	Real B_pq0(int p, int q);
	
    /** S_{pq}(0) */
    Real swapRate(int p, int q);
     
	
// CAPLET AND SWAPTION PRICES

    /** Forward price of Caplet(i). Sigma is the aggregate volatility
	 *  \f$\Sigma_i^2(0,T_i)\f$ of the caplet to expiry. See book
	 *  6.8.2 and 6.6.
	 */
    Real capletForwardPrice(int i, Real strike, Real Sigma);
	
	
	/** Forward price of Swaption(i,n) (swap on \f$[T_i,T_n]\f$).
	 *  Sigma is the aggregate volatility \f$\Sigma_{i,n}^2(0,T_i)\f$ o
	 *  of the swap rate logarithm to expiry. See book 6.8.3 and 6.7.
	 */
    Real swaptionForwardPrice(int i, Real strike, Real Sigma);


   /** Caplet on\f$[T_i,T_{i+1}]\f$.
    *  Forward price computed from current parameters of the
	*  factor loading. Method is only needed to write synthetic data.
	*  Empty implementation returns 0. See also 
	*  {@link capletForwardPrice(int, Real, Real)}.
	*  How the volatility Sigma is computed depends on the LMM type
	*  and will have to be defined in the calibrator for each type of LMM.
	*/
   virtual Real capletForwardPrice(int i, Real strike){ return 0.0; }

   
   /** Swaption on\f$[T_i,T_n]\f$, expires at \f$T_p\f$.
    *  Forward price computed from current parameters of the
	*  factor loading. Method is only needed to write synthetic data.
	*  Empty implementation returns 0. See also 
    *  {@link swaptionForwardPrice(int, Real, Real)}.
	*  How the volatility Sigma is computed depends on the LMM type
    *  and will have to be defined in the calibrator for each type of LMM.
	*/
   virtual Real swaptionForwardPrice(int i, Real strike){ return 0.0; }
   
   
// WRITING SYNTHETIC DATA
   
   /** <p>Computes the analytic forward prices  of all at the money caplets and swaptions
    *  from the current parameters (state) of the factor loading and writes them to the caplet 
    *  and swaption outfiles in the format in which the calibrator reads them for calibration. 
    *  Note that these data cannot be read in the same pass of the 
    *  program since the outfiles are associated with the ofstreams.
    *
    * <p>This needs the appropriate implementations of the functions
    * {@link capletForwardPrice(int, Real)} and {@link swaptionForwardPrice(int, Real)}
    * which depend on the type of LMM and will have to be defined in the calibrator
    * for each type of LMM.
    */  
   void writeSyntheticData();
   
   
// OBJECTIVE FUNCTION
	
	
   /** Pricing error under the current parameters of the factor loading
    *  (all caplets and coterminal swaptions).
    */
   virtual Real objectiveFunction() = 0;
   
   
   	   	
    /** The objective function as a function of the parameter values x
	 *  applied to the factor loading.
     */
   virtual Real objectiveFunction(const RealArray1D& x) = 0;
   

	
// CALIBRATION 

   
   /** Mean relative pricing error over all caplets and swaptions.
    *  Assumes the calibrated forward prices and relative errors have already 
	*  been written into the caplet and swaption data objects. 
    */
   Real meanRelativeCalibrationError();
   

     
   /** Calibrates the factor loading with nVals evaluations of the objective function
	*  (sum squared error over all caplets and coterminal swaptions) and prints
	*  the result to the files "CapletsOut.txt", "SwaptionsOut.txt".
    *
	*  @return the calibrated factor loading.
	*/
   virtual LiborFactorLoading* calibrate(int nVals) = 0;
   
   
private:

	/** Compute the caplet forward prices under the current
	 *  parameters of the factorloading and write these together with
	 *  the calibration errors into the CapletData objects.
	 */
	void recordCapletCalibration();
   
	
}; // end LiborCalibrator




/**********************************************************************************
 *
 *              StandardLMM-Calibrator
 *
 *********************************************************************************/


/** <p>Calibrator for a Libor market model which follows a dynamics described in the book,
 *  that is, <i>driftless</i> or <i>predictor-corrector</i>. Implements everything except 
 *  computation of the aggregate volatilities to expiry of Libor and swap rate logarithms
 *  which differ in the two models. 
 *
 *  <p>Very likely other models based on the same factor loading
 *  types (deterministic volatilities and constant correlations) and possibly others will
 *  be able to use this class. It all depends wether there is a useable Black-Scholes 
 *  analytic approximation for the caplet and swaption forward prices based on a forecast
 *  of the aggregate volatility to expiry of Libor and swaprate logarithms.
 *  See book, 6.6, 6.7 (predictor corrector model) and 6.8.2, 6.8.3 (driftless model)
 *  as well as Appendix, C.4 for the general Black-Scholes type approximation. 
 *
 *  <p>Currently the calibrator assumes that <i>all</i> the forward
 *  prices of caplets Caplet(i) and Swaptions Swaption(i,n), i=1,2,...,n
 *  are in the files CapletsIn, SwaptionsIn in the form described in
 *  SyntheticData/Readme.html. Nothing is checked. The routine {@link #writeSyntheticData()} 
 *  produces output in this form.
 */
class StandardLmmCalibrator : public LmmCalibrator {
	
protected:
	
	RealVector capletImpliedSigma;           // implied caplet volatility to expiration Sigma_i(0,T_i)
	
	UTRRealMatrix cvMatrix;                  // workspace for covariation matrices
	
	Real capletVol, swaptionVol;             // current aggregate volatities (cache)
	RealVector x;                            // cache, vector x for caplet and swaption aggregate vols
	                                         // see book, 6.8.3, 6.8.2 (DriftlessLMM) and 6.7 (PredictorCorrectorLMM).
	
public:
		
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

	
	
// THE OBJECTIVE FUNCTION
	
	   	
     /** The objective function as a function of the parameter values x
	  *  applied to the factor loading.
      */
     virtual Real objectiveFunction(const RealArray1D& x);
	
	
	 Real objectiveFunction();
	
	
 	 /** From the current parameters for the volatility surface and correlations
      *  set the scaling factors of the factor loading \f$k[j]=c_j\f$ to the values
      *  which reproduce the caplet prices exactly. See book, 6.11.3.
      */
     virtual void setScalingFactors() = 0;
	
	
	

// CAPLET AND SWAPTION PRICES

   /** Caplet on\f$[T_i,T_{i+1}]\f$.
    *  Forward price computed from current parameters of the
	*  factorloading.
	*/
   virtual Real capletForwardPrice(int i, Real strike) = 0;
   
   /** Swaption on\f$[T_i,T_n]\f$, expires at \f$T_i\f$.
    *  Forward price computed from current parameters of the
	*  factorloading.
	*/
   virtual Real swaptionForwardPrice(int i, Real strike) = 0;


	
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
	 
	
	/** Tests calibration on the synthetic data in the directory SyntheticData.
	 *  This directory contains synthetic caplet and swaption prices in dimensions
	 *  20,30,40,50 (number of Libor accrual intervals). These data have been produced 
	 *  using the sample LMMs returned by the method
	 * {@link LiborMarketModel#sample(int,int,int,int)}. Following up how these are
	 * constructed one can determine what parameters are used in each case.
	 * We don't check that the parameters have legal values. Please make sure they do.
	 *
	 * @param nVals number of evaluations of the objective function (suggested:500).
	 * @param dim dimension of the LMM, number of accrual periods = 20,30,40,50.
	 * @param dataLmmType type of LMM that produced the data (LiborFactorLoading::DL,PC).
	 * @param dataVolSurfaceType type of VolSurface that produced the data (VolSurface::CONST,M,JR).
	 * @param dataCorrelationType type of Correlations that produced the data (Correlations::CS,JR).
	 * @param lmmType type of LMM to be calibrated: (LiborFactorLoading::DL,PC).
	 * @param volSurfaceType type of VolSurface to be calibrated (VolSurface::CONST,M,JR).
	 * @param correlationType type of Correlations to be calibrated (Correlations::CS,JR).
	 */
	static void testCalibration
	(int nVals, int dim, 
	 int dataLmmType, int dataVolSurfaceType, int dataCorrelationType,
	 int lmmType, int volSurfaceType, int correlationType
	 );
	
   
   /** Writes the log-Libor covariation matrix CV with indices
    * \f[CV_{ij}=\int_0^{T_p}\sigma_i(u)\sigma_j(u)\rho_{ij}du,\quad p\leq i\leq j<q,\f]
	* into the field cvMatrix with the same subscripts i,j=p,...,q-1; i<=j.
	*/
   void writeCovariationMatrix(int p, int q);
   
   
   /** Calibrates the factor loading with nVals evaluations of the objective function
	*  (sum squared error over all caplets and swaptions which have been read) and prints
	*  the result to the files "CapletsOut.txt", "SwaptionsOut.txt".
    *
	*  @return the calibrated factorLoading
	*/
   LiborFactorLoading* calibrate(int nVals);
   
  
protected:


	 /** The squareroot \f$\sqrt{xCx'=(Cx,x)\f$ where $C$ is the covariation 
	  *  matrix in the field cvMatrix with indices j,k=i,...,n-1 and x is the vector
	  *  in the field x indexed as x[j], j=i,...,n-1. This quantity is used to 
	  *  compute caplet and swaption aggregate volatilities and the fields
	  *  cvMatrix and x are then set approproately. See book 6.6, 6.7,
	  *  6.8.2, 6.8.3. Note that we don't take the detour \f$\Sigma=||R'x||\f$
	  *  using a pseudo square root of the covariation matrix C.
	  */
     Real root_xCx(int i);


private:

      
	 /** Computes the aggregate caplet volatilities \f$\ol\Sigma_i(0,T_i)\f$ to expiry
	  *  implied by the Black caplet formula from the caplet market prices and writes 
	  *  them into the vector capletImpliedSigma.
	  */
	 void writeCapletImpliedSigmas();
	 
	

}; // end StandardLmmCalibrator






/**********************************************************************************
 *
 *              DriftLessLMM-Calibrator
 *
 *********************************************************************************/


/** Calibrator for a {@link DriftlessLMM}.
 */
class DriftlessLmmCalibrator : public StandardLmmCalibrator {
	
public:
		
	 DriftlessLmmCalibrator
	 (LiborFactorLoading* fl,
	  const char* capletsInFile="CapletsIn.txt", 
	  const char* swaptionsInFile="SwaptionsIn.txt", 
	  const char* capletsOutFile="CapletsOut.txt", 	 
	  const char* swaptionsOutFile="SwaptionsOut.txt"
	 ) :	
     StandardLmmCalibrator(fl,capletsInFile,swaptionsInFile,capletsOutFile,swaptionsOutFile)
     {    }
	

// CAPLET AND SWAPTION PRICES

     /** Caplet on\f$[T_i,T_{i+1}]\f$.
      *  Forward price computed from current parameters of the
  	  *  factorloading.
	  */
     Real capletForwardPrice(int i, Real strike);
   
     /** Swaption on\f$[T_i,T_n]\f$, expires at \f$T_i\f$.
      *  Forward price computed from current parameters of the
	  *  factorloading.
	  */
     Real swaptionForwardPrice(int i, Real strike);

	 
 	 /** From the current parameters for the volatility surface and correlations
      *  set the scaling factors of the factor loading \f$k[j]=c_j\f$ to the values
	  *  which reproduce the caplet prices exactly. See book, 6.11.3.
	  */
     void setScalingFactors();
	
   
     /** Reads caplets and swaptions in dimension 50 from the files "CapletsIn-dim50-DL-JR-CS.txt"
      *  and "SwaptionsIn-dim50-DL-JR-CS.txt" in the current directory and writes them to the files
	  *  "CapletsOut.txt" and "SwaptionsOut.txt".
	  */
     static void testIO();
	 
	
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
class PredictorCorrectorLmmCalibrator : public StandardLmmCalibrator {
	
public:
	
	
	PredictorCorrectorLmmCalibrator
	(LiborFactorLoading* fl,
	 const char* capletsInFile="CapletsIn.txt", 
	 const char* swaptionsInFile="SwaptionsIn.txt", 
	 const char* capletsOutFile="CapletsOut.txt", 	 
	 const char* swaptionsOutFile="SwaptionsOut.txt"
	) :	
    StandardLmmCalibrator(fl,capletsInFile,swaptionsInFile,capletsOutFile,swaptionsOutFile)
    {    }

		
     /** Caplet on\f$[T_i,T_{i+1}]\f$.
      *  Forward price computed from current parameters of the
  	  *  factorloading.
	  */
     Real capletForwardPrice(int i, Real strike);
   
     /** Swaption on\f$[T_i,T_n]\f$, expires at \f$T_i\f$.
      *  Forward price computed from current parameters of the
	  *  factorloading.
	  */
     Real swaptionForwardPrice(int i, Real strike);

	
 	 /** From the current parameters for the volatility surface and correlations
      *  set the scaling factors of the factor loading \f$k[j]=c_j\f$ to the values
	  *  which reproduce the caplet prices exactly. See book, 6.11.3.
	  */
     void setScalingFactors();
	
	
	/** Test in dimension 50. Calibrates PredictorCorrectorLMM to synthetic data 
	 *  produced by PredictorCorrectorLMM.
	 *
	 * @param nVals number of evaluations of the objective function.
	 */
	static void testCalibration(int nVals);

	 
}; // end PredictorCorrectorLmmCalibrator



			



MTGL_END_NAMESPACE(Martingale)

#endif
