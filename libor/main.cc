
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
// #include "TestFormulas.h"
// #include "Examples.h"
// #include "DirichletProblem.h"
// #include "VolatilityAndCorrelation.h"
 #include "LiborCalibrator.h"
 #include "LiborMarketModel.h"
 


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
	
	// testDownhillSimplex(5,5000);
	// testBFGS(5,500);
	// testSobolSearch(50,5000);
	
	
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
 *                       VOLATILITY AND CORRELATION
 *
*******************************************************************************/


    // testVolSurfaceIntegrals(1000000,0.5);
	
	
/*******************************************************************************
 *
 *                        LIBOR PROCESS, LMM
 *
*******************************************************************************/
	
	
	
	
	// testFactorLoading(20);
	// testLiborFactorLoadingFactorization(60,3);    // rank 3 factorization
	
	// testLmmPaths(20);

	// Examples::liborPathTiming(30,20000);         // 5000 paths in dimension 30
	
	
	// anomalies in high dimensions (n>=70) for all but DriftlessLMM,
	    // int lmmType=LiborMarketModel::DL;  // DL, LFDL, PC, FPC
	    // int volType=VolSurface::M;        // JR, M, CONST
        // int corrType=Correlations::CS;     // JR, CS
	// testCapletPrice(lmmType,volType,corrType); 
	// testSwaptionPrice(lmmType,volType,corrType);
	// testCallOnBondPrice(lmmType,volType,corrType);
	// testCallOnZeroCouponBondPrice(lmmType,volType,corrType);
	

	
	
	
/*******************************************************************************
 *
 *                        LIBOR CALIBRATION
 *
*******************************************************************************/

    // DriftlessLmmCalibrator::testIO();
	
    // Write the synthetic data in the directory SyntheticData
	// StandardLmmCalibrator::writeSyntheticDataSample();

	int nVals=500,                               // number of evaluations of the objective function
	    n=50,                                    // 20,30,40,50 dimension in which we calibrate
	    dataLmmType=LiborMarketModel::PC,        // type of LMM the data came from: DL,PC
		dataVolType=VolSurface::M,              // type of VolSurface used in the generation of data: CONST,JR,M
		dataCorrType=Correlations::JR,           // type of correlations used in the generation of data: JR,CS
		lmmType=LiborMarketModel::PC,            // type of LMM which will be calibrated: PC,DL
		volType=VolSurface::JR,                   // type of VolSurface which will be calibrated: CONST,JR,M
		corrType=Correlations::CS;               // type of correlations which will be calibrated: JR,CS
		
     StandardLmmCalibrator::testCalibration
	 (nVals,n,dataLmmType,dataVolType,dataCorrType,lmmType,volType,corrType);


	
	
/*******************************************************************************
 *
 *                        LMM LATTICES
 *
*******************************************************************************/
	
	
	 
	 // LatticeSwaption3F::test(15,19,30);
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

//<---------------missing-------------->
	// testSobolSequence(5,16383);
	// testSobolInit(30);


/*******************************************************************************
 *
 *             STOCHASTIC PROCESSES, PATH FUNCTIONALS
 *
*******************************************************************************/	


	 // int dim=10, T=10000, nPath=10000; 
	 // Real dt=0.01;
	 // Examples::brownianMotionInBall(dim,T,dt,nPath); <------------Subscript out of range exception--------->	
	
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
