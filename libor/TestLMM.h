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




#include "Derivatives.h"


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
	LiborFactorLoading* fl;
	cout << "\nTesting CS_FactorLoading: " << endl;
	Timer watch;
	watch.start();
	fl=CS_FactorLoading::sample(n);
	watch.stop();
	watch.report("CS_FactorLoading set up.");

	fl->selfTest();
	
	cout << "\n\n\nTesting JR_FactorLoading: " << endl;
	watch.start();
	fl=JR_FactorLoading::sample(n);
	watch.stop();
	watch.report("JR_FactorLoading set up.");

	fl->selfTest();
	
} // end testFactorLoading


/** Allocates sample LiborFactorLoadings of type CS and JR in dimension n and sets
 *  Then runs the rank r factorization test {@link UTRMatrix#testFactorization(int r)}
 *  on all log-Libor covariation matrices 
 *  {@link LiborFactorLoading#logLiborCovariationMatrix(int t)}, t=0,1,...,n-3.
 */
static void testLiborFactorLoadingFactorization(int n, int r)
{

	cout << "\n\n  CS-FactorLoading:"  << endl;
	LiborFactorLoading* 
	fl = CS_FactorLoading::sample(n);
	for(int t=0;t<n-2;t++){
		
		UTRRealMatrix& CV=fl->logLiborCovariationMatrix(t);
		CV.testFactorization(r);
    }
	
	cout << "\n\n\n\n\n  JR-FactorLoading:" << endl;
	fl = JR_FactorLoading::sample(n);
	for(int t=0;t<n-2;t++){
		
		UTRRealMatrix& CV=fl->logLiborCovariationMatrix(t);
		CV.testFactorization(r);
    }
} // end factorizationTest

 
/*******************************************************************************
 *
 *                        LMM TESTS
 *
 ******************************************************************************/

/** Sets up a LognormalLMM in dimension n and prints
 *  the drift step vectors and covariationmatrixRoots.
 */
void testLognormalLMM(int n){ LognormalLMM::test(n); }


 
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
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testSwaptionPrice(bool PC=false, bool FPC=false, bool LS=true)
{ 
    for(int n=30;n<70;n+=20)
	for(Real delta=0.1;delta<0.6;delta+=0.4){
		
		Swaption* swpn=Swaption::sample(n,delta);
	    swpn->priceTest(delta,PC,FPC,LS);
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
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testCapletPrice(bool PC=false, bool FPC=false, bool LS=true)
{ 
    int n=80; Real delta=0.25;
	Caplet* cplt=Caplet::sample(n,delta);
	cplt->priceTest(delta,PC,FPC,LS);
       
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
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testCallOnBondPrice(bool PC=false, bool FPC=false, bool LS=true)
{ 
    for(int n=30;n<60;n+=20)
	for(Real delta=0.25;delta<0.6;delta+=0.4){
		
		BondCall*  bc=BondCall::sample(n,delta);
	    bc->priceTest(delta,PC,FPC,LS);
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
void testCallOnZeroCouponBondPrice(bool PC=false, bool FPC=false, bool LS=true)
{ 
    for(int n=30;n<60;n+=20)
	for(Real delta=0.25;delta<0.6;delta+=0.4){
		
		BondCall*  bc=BondCall::sampleCallOnZeroCouponBond(n,delta);
	    bc->priceTest(delta,PC,FPC,LS);
	}
       
} // end testCallOnZeroCouponBondPrice




/*******************************************************************************
 *
 *                  LIBORS WITH REDUCED NUMBER OF FACTORS
 *
 ******************************************************************************/




	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 