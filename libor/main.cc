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


/* Created by Anjuta version 0.1.9 */
/*	This file will not be overwritten */


// #include "FinMath.h"
// #include "TestMatrix.h"
// #include "TestProbability.h"
// #include "TestLMM.h"
// #include "TestOptimizers.h"
// #include "LatticeOption.h"
// #include "TestFormula.h"
// #include "Examples.h"
// #include "DirichletProblem.h"
#include "LiborCalibrator.h"


using namespace Martingale;

/** Various examples and tests.
 */
int main()
{
	
/*******************************************************************************
 *
 *                        OPTIMIZATION
 *
*******************************************************************************/
	
	// testDownhillSimplex(5,500);
	// testBFGS(5,500);
	// testSobolSearch(5,500);
	
	
/*******************************************************************************
 *
 *                        RANDOM NUMBERS, VARIABLES, VECTORS
 *
*******************************************************************************/
	
    // testSTNCovarianceMatrix(600,50000);
    // testControlledSTN(1000000);
    // testPathFunctional(9,15,100000);
	
/*******************************************************************************
 *
 *                       MATRICES
 *
*******************************************************************************/
    
	// testMatrixInverse();
	// Examples::timeMatrixMultiply(1000,1);          
	// testMatrixMultiply();
	// testMatrixExponentials();
	// testComplexExponential();
	// TnT::factorizationTest(80);
	
/*******************************************************************************
 *
 *                        LIBOR PROCESS, LMM
 *
*******************************************************************************/
	
	// testFactorLoading(20);
	// testLiborFactorLoadingFactorization(60,3);    // rank 3 factorization

	// Examples::liborPathTiming(80,10000);         // 10000 paths in dimension 80
	// Examples::liborApproximationError(124, 120, 120, 20000);
	
	/* anomalies in high dimensions (n>=70) for all but DriftlessLMM, 
	bool LS=true,    // DriftlessLMM
	     PC=false,   // PredictorCorrectorLMM
		 FPC=false;  // FastPredictorCorrectorLMM
	//testCapletPrice(LS,PC,FPC); 
	//testSwaptionPrice(LS,PC,FPC);
	//testCallOnBondPrice(LS,PC,FPC);
	//testCallOnZeroCouponBondPrice(LS,PC,FPC);
	*/
	// testLiborMeans(30,25,20,50000);
	// testLognormalLMM(20);
	
	
	
/*******************************************************************************
 *
 *                        LIBOR CALIBRATION
 *
*******************************************************************************/

    // test the calibrator for the driftless LMM on a 
	// constant volatility LiborFactorLoading
    // DriftlessLmmCalibrator::test(20,LiborFactorLoading::CV);
	 DriftlessLmmCalibrator::writeSyntheticDataSample();
	
	// PredictorCorrectorLmmCalibrator::writeSyntheticDataSample();
	
	
/*******************************************************************************
 *
 *                        LMM LATTICES
 *
*******************************************************************************/
	
	
	 
	 // LatticeSwaption3F::test(15,17,30);
     // ConstVolLmmLattice2F::test(50);
	 
	 
	
/*******************************************************************************
 *
 *                 LATTICES FOR ASSET BASKETS
 *
*******************************************************************************/
	 
	 // BasketLattice3F::test(5,50);
	
	
/*******************************************************************************
 *
 *                       SOBOL
 *
*******************************************************************************/
	
	// testSobolSequence(5,16383);
	// testSobolInit(30);


/*******************************************************************************
 *
 *             STOCHASTIC PROCESSES, PATH FUNCTIONALS
 *
*******************************************************************************/	

	
	// int dim=10, T=10000, nPath=10000; 
	// Real dt=0.01;
	// Examples::brownianMotionInBall(dim,T,dt,nPath);
	
	// int t=20, T=100, nPath=500000;
	// VectorBrownianMotion::testPathFunctional(t,T,nPath);
	
	// int dim=30, T=10;
	// DirichletProblemExample::runExample(dim,T);
	


/*******************************************************************************
 *
 *             BOOK FORMULAS
 
 *
*******************************************************************************/	

	
	// testExponentialIntegralFormulas();
     
	return 0;
}
