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
#include "LmmLattice.h"
#include "LatticeOption.h"

MTGL_BEGIN_NAMESPACE(Martingale)

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
		std::cout << "\n\n20 Libor paths, LMM type: " << lmm->getType() 
		          << endl << endl;
		for(int path=0;path<2;path++){
			
			lmm->newPath();
		    for(int t=0;t<n-5;t++) std::cout << lmm->XL(n-4,t) << ", ";
			std::cout << lmm->XL(n-4,n-5) << endl;
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
    Swaption* swpn=Swaption::sample(n,lmmType,volType,corrType);
	swpn->testPrice();
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
	cplt->testPrice();
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
	bc->testPrice();
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
	bc->testPrice();
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
	
	    std::cout << "\n\nBuilding and testing a lattice for the driftless libor market model:"
	              << "\nEnter dimension n of Libor process: n = ";
	    int n; std::cin >> n;
	    std::cout << "Enter number r of factors (2 or 3): r = ";
	    int r; std::cin >> r;
	    std::cout << "Light (h=5) or heavyweight (h=7) lattice? h = ";
     	int h; std::cin >> h;
	    switch(r*h){
		
		    case 10 : LiteLmmLattice2F::test(n); break;
	        case 15 : LiteLmmLattice3F::test(n); break;
		    case 14 : HeavyLmmLattice2F::test(n); break;
	        case 21 : HeavyLmmLattice3F::test(n); break;
			default : std::cout << "\n\n\nYou did not enter the parameters right.";
		}
		
		std::cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		std::cin >> do_again;
	}
} // end testLmmlattice


/** Swaption price in a lattices for the driftless Libor Market Model compared to 
 *  analytic and Monte Carlo prices.
 *  Asks for the type of lattice (light/heavy), number of factors, number of time steps per 
 *  Libor accrual interval and swap period [T_p,T_q]. Then sets up sample swaption 
 *  (constant volatility and JR correlations, book, 6.11.2) and computes the forward price.
 */
void testLatticeSwaption()
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
	    std::cout << "\n\nSwaption price in a lattice for the driftless libor market model:"
	              << "\nSwap interval: [T_p,T_q]."
		          << "\nEnter p = ";
	    int p; std::cin >> p;
		std::cout << "Enter q = ";
		int q; std::cin >> q;
	    std::cout << "Enter number r of factors (2 or 3): r = ";
	    int r; std::cin >> r;
	    std::cout << "Light (h=5) or heavyweight (h=7) lattice? h = ";
     	int h; std::cin >> h;
		std::cout << "Number of time steps in each Libor accrual interval = ";
		int nSteps; std::cin >> nSteps;
	    switch(r*h){
		
		    case 10 : LiteLatticeSwaption2F::test(p,q,nSteps); break;
	        case 15 : LiteLatticeSwaption3F::test(p,q,nSteps); break;
		    case 14 : HeavyLatticeSwaption2F::test(p,q,nSteps); break;
	        case 21 : HeavyLatticeSwaption3F::test(p,q,nSteps); break;
			default : std::cout << "\n\n\nYou did not enter the parameters right.";
		}
		
		std::cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		std::cin >> do_again;
	}
} // end testLatticeSwaption


	
	     
	
MTGL_END_NAMESPACE(Martingale)

#endif
 