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

#ifndef martingale_testlmm_h    
#define martingale_testlmm_h

#include "TypedefsMacros.h"
#include "Derivatives.h"
#include "LiborFactorLoading.h"
#include "LiborMarketModel.h"            // for inclusion to main.cc

MTGL_BEGIN_NAMESPACE(Martingale)

/** A collection of free standing short test programs.
 */
 

/*******************************************************************************    
    
              TEST OF FACTORLOADINGS 
	
*******************************************************************************/
		 
/** Sets up a sample CS_FactorLoading in dimension n and then runs
 *  the self test {@link FactorLoading#selfTest}.
 */
void testFactorLoading(int n)
{
	LiborFactorLoading* fl=0;
	for(int volType=0;volType<3;volType++)
	for(int corrType=0;corrType<2;corrType++){
		
	   Timer watch;
	   watch.start();
	   fl=LiborFactorLoading::sample(n,volType,corrType);
	   fl->selfTest();
	   watch.stop();
	   watch.report(" ");
	}
	
} // end testFactorLoading


/** Allocates sample LiborFactorLoadings of type CS and JR in dimension n and sets
 *  Then runs the rank r factorization test {@link UTRMatrix#testFactorization(int r)}
 *  on all log-Libor covariation matrices 
 *  {@link LiborFactorLoading#logLiborCovariationMatrix(int t)}, t=0,1,...,n-3.
 */
static void testLiborFactorLoadingFactorization(int n, int r)
{

	LiborFactorLoading* fl=0;
	for(int volType=0;volType<3;volType++)
	for(int corrType=0;corrType<2;corrType++){
		
	   std::cout << "\nTesting Libor factor loading," 
		         << "\ncorrelation matrices low rank factorization:" 
		         << "\nVolatility surface type = VolSurface::" << volType
		         << "\nCorrelation type = Correlations::" << corrType;
	   Timer watch;
	   watch.start();
	   fl=LiborFactorLoading::sample(n,volType,corrType);
	   for(int t=0;t<n-2;t++){
		
		   const UTRRealMatrix& CV=fl->logLiborCovariationMatrix(t);
		   CV.testFactorization(r);
       }
	   watch.stop();
	   watch.report(" ");
	}
} // end factorizationTest

 
/*******************************************************************************
 *
 *                        LMM TESTS
 *
 ******************************************************************************/


/** Prints 20 paths of \f$X_{n-4}\$. All Libor market model types.
 */
void testLmmPaths(int n)
{
	for(int lmmType=0;lmmType<4;lmmType++)
	for(int volType=0;volType<3;volType++)
	for(int corrType=0;corrType<2;corrType++){
		
		LiborMarketModel* lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
		std::cout << "\n\n20 Libor paths, LMM type: " << lmm->modelType() 
		          << endl << endl;
		for(int path=0;path<2;path++){
			
			lmm->newPath();
		    for(int t=0;t<n-5;t++) std::cout << lmm->XL(n-4,t) << ", ";
			std::cout << lmm->XL(n-4,n-5) << endl;
		}
	}
}	

 
// SWAPTION PRICE

/** <p>Allocates  LMM in dimension n=30,50 and prices the at the money
 *  payer swaption \f$swpn(T_p,[T_p,T_q])\f$ exercisable at time \f$T_p\f$
 *  where \f$p=n/3, q=n\f$. 
 *
 *  <p>For each dimension the test is carried out for accrual intervals of length
 *  \f$\delta_j=0.1,0.5\f$. With increasing delta the predictor corrector 
 *  algorithm becomes less accurate.
 *
 *  <p>The Monte Carlo forward price of the swaption is computed from 20,000
 *  Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq p\f$. This is a good stress test for the LMM.
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testSwaptionPrice(int lmmType, int volType, int corrType)
{ 
    for(int n=30;n<70;n+=20){
		
		Swaption* swpn=Swaption::sample(n,lmmType,volType,corrType);
	    swpn->testPrice();
	}
       
} // end testSwaptionPrice




 
// CAPLET PRICE

/** <p>Allocates  LMM in dimension n=20,40,60,80 and prices 
 *  the at the money caplet \f$cplt([T_i,T_{i+1}])\f$,
 *  where \f$i=n/3\f$. This is intended to stress the predictor corrector algorithm
 *  as inaccuracies compound when the caplet payoff is transported forward to
 *  the horizon (accrual factors).
 *
 *  <p>For each dimension the test is carried out for accrual intervals of length
 *  \f$\delta=0.1,0.5\f$. With increasing delta the predictor corrector 
 *  algorithm becomes less accurate.
 *
 *  <p>The Monte Carlo forward price of this caplet at time \f$T_n\f$ is computed from 
 *  a sample of 20000 Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq i+1\f$.</p>
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCapletPrice(int lmmType, int volType, int corrType)
{ 
     for(int n=30;n<70;n+=20){
		
		LiborDerivative* cplt=Caplet::sample(n,lmmType,volType,corrType);
	    cplt->testPrice();
	}
} // end testCapletPrice


// BONDCALL PRICE

/** <p>Allocates  LMM in dimension n=30,50 and prices the at the money call
 *  on a bond along \f$[T_p,T_q]\f$ exercisable at time \f$T_p\f$ where 
 *  \f$p=n/3, q=2*n/3\f$. Coupons initialized randomly with \f$c_j\in[-0.5,1.5]\f$.
 *
 *  <p>For each dimension the test is carried out for accrual intervals of length
 *  \f$\delta_j=0.1,0.5\f$. With increasing delta the predictor corrector 
 *  algorithm becomes less accurate.
 *
 *  <p>The Monte Carlo forward price is computed from 20,000
 *  Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq p\f$. This is a good stress test for the LMM.
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCallOnBondPrice(int lmmType, int volType, int corrType)
{ 
    for(int n=30;n<70;n+=20){
		
		BondCall* bc=BondCall::sample(n,lmmType,volType,corrType);
	    bc->testPrice();
	}       
} // end testBondCallPrice



/** Same as {@link testCallOnBondCallPrice()} but the bond is now a
 *  zero coupon bond maturing at time \f$T_i\f$. The call expires at
 *  time \f$T_{i-1}\f$. This is a worst case for the assumptions of the
 *  analytic price formulas.
 *
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testCallOnZeroCouponBondPrice(int lmmType, int volType, int corrType)
{ 
    for(int n=30;n<60;n+=20)
	for(Real delta=0.25;delta<0.6;delta+=0.4){
		
		BondCall*  bc=BondCall::sampleCallOnZeroCouponBond(n,lmmType,volType,corrType);
	    bc->testPrice();
	}
       
} // end testCallOnZeroCouponBondPrice




/*******************************************************************************
 *
 *                  LIBORS WITH REDUCED NUMBER OF FACTORS
 *
 ******************************************************************************/




	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 