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
#include "LiborMarketModel.h"            // for inclusion to main.cc
#include "VolatilityAndCorrelation.h"    // for inclusion to main.cc

MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Test_LMM)

/** A collection of free standing short test programs.
 */
 

/*******************************************************************************    
    
              TEST OF FACTORLOADINGS 
	
*******************************************************************************/
		 
/** Sets up a sample CS_FactorLoading in dimension n and then runs
 *  the self test {@link FactorLoading#selfTest}.
 */
void testLiborFactorLoading(int n);


/** Allocates sample LiborFactorLoadings of type CS and JR in dimension n and sets
 *  Then runs the rank r factorization test {@link UTRMatrix#testFactorization(int r)}
 *  on all log-Libor covariation matrices 
 *  {@link LiborFactorLoading#logLiborCovariationMatrix(int t)}, t=0,1,...,n-3.
 */
void testLiborFactorLoadingFactorization(int n, int r);

 
/*******************************************************************************
 *
 *                        TESTS IN GENERAL LMM
 *
 ******************************************************************************/


/** Prints 20 paths of \f$X_{n-4}\$. All Libor market model types.
 */
void testLmmPaths(int n);
 
 
// SWAPTION PRICE

/** <p>Allocates  LMM in dimension n and prices the at the money
 *  payer swaption \f$swpn(T_p,[T_p,T_q])\f$ exercisable at time \f$T_p\f$
 *  where \f$p=n/3, q=n\f$. 
 *
 *  <p>The Monte Carlo forward price of the swaption is computed from 8192=2^{13}
 *  Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq p\f$. This is a good stress test for the LMM.
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testSwaptionPrice(int lmmType, int volType, int corrType);


// CAPLET PRICE

/** <p>Allocates  LMM in dimension n and prices 
 *  the at the money caplet \f$cplt([T_i,T_{i+1}])\f$,
 *  where \f$i=n/3\f$. This is intended to stress the predictor corrector algorithm
 *  as inaccuracies compound when the caplet payoff is transported forward to
 *  the horizon (accrual factors).
 *
 *  <p>The Monte Carlo forward price of this caplet at time \f$T_n\f$ is computed from 
 *  a sample of 8192=2^{13} Libor paths and compared to the analytic price both using 
 *  an MC and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq i+1\f$.</p>
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCapletPrice(int lmmType, int volType, int corrType);


// BONDCALL PRICE

/** <p>Allocates  LMM in dimension n (user supplied) and prices the at the money 
 *  call on a bond along \f$[T_p,T_q]\f$ exercisable at time \f$T_p\f$ where 
 *  \f$p=n/3, q=2*n/3\f$. Coupons initialized randomly with \f$c_j\in[-0.5,1.5]\f$.
 *
 *  <p>For each dimension the test is carried out for accrual intervals of length
 *  \f$\delta_j=0.1,0.5\f$. With increasing delta the predictor corrector 
 *  algorithm becomes less accurate.
 *
 *  <p>The Monte Carlo forward price is computed from 8192=2^{13}
 *  Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq p\f$. This is a good stress test for the LMM.
 *
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCallOnBondPrice(int lmmType, int volType, int corrType);



/** Same as {@link testCallOnBondCallPrice()} but the bond is now a
 *  zero coupon bond maturing at time \f$T_i\f$. The call expires at
 *  time \f$T_{i-1}\f$. This is a worst case for the assumptions of the
 *  analytic price formulas.
 *
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testCallOnZeroCouponBondPrice(int lmmType, int volType, int corrType);




/*******************************************************************************
 *
 *                  LATTICE TESTS
 *
 ******************************************************************************/

/** Asks for the number of factors and dimension n and then builds a 
 *  corresponding Lmmlattice with n-3 time steps (one time step per accrual interval)
 *  and runs the selftest on the lattice.
 */
void testLmmLattice();

 
/*******************************************************************************
 *
 *                        TESTS IN THE DRIFTLESS LMM
 *
 ******************************************************************************/


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testSwaption();


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testBermudanSwaption();


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCallOnBond();


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCallOnZeroCouponBond();


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCaplet();


/** Pricing of some Libor derivatives in the default driftless Libor Market 
 *  Model (constant volatility and JR correlations, book, 6.11.2).
 *  Analytic, default lattice, Monte Carlo and QMC prices with
 *  and without control variates. Asks for user input in a
 *  perpetual loop. Reports error relative to the analytic price. Note however 
 *  that the analytic price itself is an approximation.
 */
void testLiborDerivative();

	
	     
MTGL_END_NAMESPACE(Test_LMM)	
MTGL_END_NAMESPACE(Martingale)

#endif
 