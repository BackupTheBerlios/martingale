
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



// INCLUDES IN ORDER


#include "TestOptimizers.h"
#include "TestProbability.h"
#include "FinMath.h"
#include "Examples.h"
#include "TestMatrix.h"
#include "VolatilityAndCorrelation.h"
#include "TestLMM.h"
#include "LiborCalibrator.h"
#include "ExercisePlot.h"
#include "TestFormulas.h"
#include "DirichletProblem.h"
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
	
	// testDownhillSimplex(5,500);
	// testBFGS(5,500);
	// testSobolSearch(5,500);
	
	
/*******************************************************************************
 *
 *                        RANDOM NUMBERS, VARIABLES, VECTORS
 *
*******************************************************************************/
	
    // testSTNCovarianceMatrix(60,10000);
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
	// testSymmetricMatrixExponentials();
	// testComplexExponential();
	
	
	
	
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
	
	
	
	
	 // Test_LMM::testLiborFactorLoading(20);
	 // Test_LMM::testLiborFactorLoadingFactorization(60,3);    // rank 3 factorization
	
	 // Test_LMM::testLmmPaths(20);

	 // Examples::liborPathTiming(30,2000);         // 20000 paths in dimension 30
	
	
	// anomalies in high dimensions (n>=70) for all but DriftlessLMM,
/*
	     int lmmType=LiborMarketModel::PC;  // DL, LFDL, PC, FPC
	     int volType=VolSurface::JR;        // JR, M, CONST
         int corrType=Correlations::CS;     // JR, CS
*/
	// Test_LMM::testCapletPrice(lmmType,volType,corrType); 
	// Test_LMM::testSwaptionPrice(lmmType,volType,corrType);
	// Test_LMM::testCallOnBondPrice(lmmType,volType,corrType);
	// Test_LMM::testCallOnZeroCouponBondPrice(lmmType,volType,corrType);
	
	// Libor derivatives in the driftless LMM with constant volatilities, 
	// includes lattice pricing
    // Test_LMM::testLiborDerivative();
	
	

	
	
	
/*******************************************************************************
 *
 *                        LIBOR CALIBRATION
 *
*******************************************************************************/

    // DriftlessLmmCalibrator::testIO();
	
    // Write the synthetic data in the directory SyntheticData
	// StandardLmmCalibrator::writeSyntheticDataSample();
/*
	int nVals=500,                               // number of evaluations of the objective function
	    n=50,                                    // 20,30,40,50 dimension in which we calibrate
	    dataLmmType=LiborMarketModel::PC,        // type of LMM the data came from: DL,PC
		dataVolType=VolSurface::M,              // type of VolSurface used in the generation of data: CONST,JR,M
		dataCorrType=Correlations::JR,           // type of correlations used in the generation of data: JR,CS
		lmmType=LiborMarketModel::DL,            // type of LMM which will be calibrated: PC,DL
		volType=VolSurface::JR,                   // type of VolSurface which will be calibrated: CONST,JR,M
		corrType=Correlations::CS;               // type of correlations which will be calibrated: JR,CS
		
     StandardLmmCalibrator::testCalibration
	 (nVals,n,dataLmmType,dataVolType,dataCorrType,lmmType,volType,corrType);
*/

	
	
/*******************************************************************************
 *
 *                        LMM LATTICES
 *
*******************************************************************************/
	

     // Test_LMM::testLmmLattice();
     // plotBermudanExercise();
	 
	
	
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
	 // Examples::brownianMotionInBall(dim,T,dt,nPath); //<----Subscript out of range exception---->	
	

	 DirichletProblemExample::runExample();
	


/*******************************************************************************
 *
 *             BOOK FORMULAS
 
 *
*******************************************************************************/	

	
	// testExponentialIntegralFormulas();
     
	return 0;
}
