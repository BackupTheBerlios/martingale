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
#include "BermudanOption.h"
#include "LiborFactorLoading.h"
#include "LiborMarketModel.h"            // for inclusion to main.cc
#include "LmmLattice.h"
#include <iostream>

using std::cout;
using std::endl;


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
void testLiborFactorLoading(int n)
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
		
	   cout << "\nTesting Libor factor loading," 
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
 *                        TESTS IN GENERAL LMM
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
		cout << "\n\n20 Libor paths, LMM type: " << *(lmm->getType()) 
		     << endl << endl;
		for(int path=0;path<2;path++){
			
			lmm->newPath();
		    for(int t=0;t<n-5;t++) std::cout << lmm->XL(n-4,t) << ", ";
			cout << lmm->XL(n-4,n-5) << endl;
		}
	}
}	

 
// SWAPTION PRICE

/** <p>Allocates  LMM in dimension n and prices the at the money
 *  payer swaption \f$swpn(T_p,[T_p,T_q])\f$ exercisable at time \f$T_p\f$
 *  where \f$p=n/3, q=n\f$. 
 *
 *  <p>The Monte Carlo forward price of the swaption is computed from 20,000
 *  Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq p\f$. This is a good stress test for the LMM.
 *
 * @param n dimension of the Libor process (number of compounding periods).
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testSwaptionPrice(int n, int lmmType, int volType, int corrType)
{ 
	Timer watch; watch.start();
	int p=n/3, q=n;
    Swaption* swpn=Swaption::sample(p,q,lmmType,volType,corrType);
	swpn->testPrice(20000);
	watch.stop(); watch.report(" ");	
} 




 
// CAPLET PRICE

/** <p>Allocates  LMM in dimension n and prices 
 *  the at the money caplet \f$cplt([T_i,T_{i+1}])\f$,
 *  where \f$i=n/3\f$. This is intended to stress the predictor corrector algorithm
 *  as inaccuracies compound when the caplet payoff is transported forward to
 *  the horizon (accrual factors).
 *
 *  <p>The Monte Carlo forward price of this caplet at time \f$T_n\f$ is computed from 
 *  a sample of 20000 Libor paths and compared to the analytic price both using an MC
 *  and a QMC dynamics. The forward transporting and discounting involves 
 *  all Libors \f$L_j, j\geq i+1\f$.</p>
 *
 * @param n dimension of the Libor process (number of compounding periods).
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCapletPrice(int n, int lmmType, int volType, int corrType)
{
	Timer watch; watch.start();
	LiborDerivative* cplt=Caplet::sample(n,lmmType,volType,corrType);
	cplt->testPrice(20000);
	watch.report(" ");
} 


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
 * @param n dimension of the Libor process (number of compounding periods).
 * @param lmmType type of Libor market model: {@link LiborMarketModel}::DL,PC,FPC.
 * @param volType type of volatility surface: {@link VolSurface}::JR,M,CONST.
 * @param corrType type of log-Libor correlations: {@link Correlations}::JR,CS.
 */
void testCallOnBondPrice(int n, int lmmType, int volType, int corrType)
{ 
	Timer watch; watch.start();
	BondCall* bc=BondCall::sample(n,lmmType,volType,corrType);
	bc->testPrice(20000);
	watch.report(" ");     
} 



/** Same as {@link testCallOnBondCallPrice()} but the bond is now a
 *  zero coupon bond maturing at time \f$T_i\f$. The call expires at
 *  time \f$T_{i-1}\f$. This is a worst case for the assumptions of the
 *  analytic price formulas.
 *
 * @param n dimension of the Libor process (number of compounding periods).
 * @param PC use a {@link PredictorCorrectorLMM} Libor market model.
 * @param FPC use a {@link FastPredictorCorrectorLMM} Libor market model.
 * @param LS use a {@link LightSpeedLMM} Libor market model.
 */
void testCallOnZeroCouponBondPrice(int n, int lmmType, int volType, int corrType)
{ 
	Timer watch; watch.start();		
	BondCall*  bc=BondCall::sampleCallOnZeroCouponBond(n,lmmType,volType,corrType);
	bc->testPrice(20000);
	watch.report(" ");
} 




/*******************************************************************************
 *
 *                  LATTICE TESTS
 *
 ******************************************************************************/

/** Asks for the number of factors and dimension n and then builds a 
 *  corresponding Lmmlattice with n-3 time steps (one time step per accrual interval)
 *  and runs the selftest on the lattice.
 */
void testLmmLattice()
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
	    cout << "\n\nBuilding and testing a lattice for the driftless libor market model:"
	         << "\nEnter dimension n of Libor process: n = ";
	    int n; cin >> n;
	    cout << "Enter number r of factors (2 or 3): r = ";
	    int r; cin >> r;
		if((r!=2)&&(r!=3))
			cout << "\n\n\nNumber of factors must be two or three.";
		else LmmLattice::test(r,n); 
		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}
} // end testLmmlattice



 
/*******************************************************************************
 *
 *                        TESTS IN THE DRIFTLESS LMM
 *
 ******************************************************************************/


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testSwaption()
{
    cout << "\nSwap interval: [T_p,T_q]."
		 << "\nEnter p = ";
	int p; cin >> p;
    cout << "Enter q = ";
    int q; cin >> q;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    Swaption* swpn = Swaption::sample(p,q); 
	swpn->testPrice(nPath);
}


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testBermudanSwaption()
{
    cout << "\nSwap interval: [T_p,T_q]."
		 << "\nEnter p = ";
	int p; cin >> p;
    cout << "Enter q = ";
    int q; cin >> q;
	cout << "Enter number of Libor paths for Monte Carlo simulation: ";
    int nPath; cin >> nPath;
	cout << "Enter number of training paths for the trigger:  ";
    int paths; cin >> paths;
	
	bool verbose=false;
    BermudanSwaption* bswpn = BermudanSwaption::sample(p,q,paths,verbose); 
	bswpn->testPrice(nPath);
}


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCallOnBond()
{
    cout << "\nCoupons will be intialized randomly."
	     << "\nInterval [T_p,T_q] for which coupons are received:"
		 << "\nEnter p = ";
	int p; cin >> p;
    cout << "Enter q = ";
    int q; cin >> q;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    BondCall* bondcall = BondCall::sample(p,q); 
	bondcall->testPrice(nPath);
}


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCallOnZeroCouponBond()
{
    cout << "\nZero coupon bond matures at T_p:"
		 << "\nEnter p = ";
	int p; cin >> p;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    BondCall* bondcall = BondCall::sampleCallOnZeroCouponBond(p);
	bondcall->testPrice(nPath);
}


/** Auxilliary function for {@link testLiborDerivative()}.*/
void testCaplet()
{
    cout << "\nCaplet on [T_i,T_{i+1}], i=n/3"
	     << "\nEnter dimension n of Libor process, n =";
	int n; cin >> n;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    Caplet* cplt = Caplet::sample(n); 
	cplt->testPrice(nPath);
}



/** Pricing of some Libor derivatives in the default driftless Libor Market 
 *  Model (constant volatility and JR correlations, book, 6.11.2).
 *  Analytic, default lattice, Monte Carlo and QMC prices with
 *  and without control variates. Asks for user input in a
 *  perpetual loop. Reports error relative to the analytic price. Note however 
 *  that the analytic price itself is an approximation.
 */
void testLiborDerivative()
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
	    cout << "\n\nPricing of at the money Libor derivatives"
		     << "\nin the default driftless Libor Market Model."
		     << "\nChoose derivative:"
		     << "\nCaplet (0)"
		     << "\nSwaption (1)"
		     << "\nBermudan swaption (2)"
		     << "\nCall on bond (3)"
		     << "\nCall on zero coupon bond (4)" 
		     << "\nDerivative = ";
		     
		int derivative; cin>>derivative;
		switch(derivative){
			
			case 0  : testCaplet(); break;
			case 1  : testSwaption(); break;
			case 2  : testBermudanSwaption(); break;
			case 3  : testCallOnBond(); break;
			case 4  : testCallOnZeroCouponBond(); break;
			default : 
			cout << "\nDerivative not recognized. Defaulting to swaption.";
			testSwaption();
		}
		
		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}
} // end testLiborDerivative


	
	     
MTGL_END_NAMESPACE(Test_LMM)	
MTGL_END_NAMESPACE(Martingale)

#endif
 