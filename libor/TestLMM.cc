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


#include "TestLMM.h"
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


void 
testLiborFactorLoading(int n)
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


void 
testLiborFactorLoadingFactorization(int n, int r)
{

	LiborFactorLoading* fl=0;
	for(int volType=0;volType<3;volType++)
	for(int corrType=0;corrType<2;corrType++){
		
	   printStars();	
	   cout << "\n\nTesting Libor factor loading," 
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

 

void 
testLmmPaths(int n)
{
	for(int lmmType=0;lmmType<4;lmmType++)
	for(int volType=0;volType<3;volType++)
	for(int corrType=0;corrType<2;corrType++){
		
		LiborMarketModel* lmm=LiborMarketModel::sample(n,lmmType,volType,corrType);
		printStars();
		cout << "\n\n5 Libor paths, LMM type: " << *(lmm->getType()) 
		     << endl << endl;
		for(int path=0;path<5;path++){
			
			lmm->newPath();
		    for(int t=0;t<n-5;t++) std::cout << lmm->XL(n-4,t) << ", ";
			cout << lmm->XL(n-4,n-5) << endl;
		}
	}
}	

 

void 
testSwaptionPrice(int lmmType, int volType, int corrType)
{ 
	int do_again=1;
	// main loop
	while(do_again==1){
	
		printStars();
	    cout << "\n\nPricing swaption along [T_p,T_n] "
	         << "in LMM of dimension n."
		     << "\n\nEnter p = ";
	    int p; cin >> p;
		cout << "Enter n = ";
		int n; cin >> n;
		
    	Timer watch; watch.start();
        Swaption* swpn=Swaption::sample(p,n,lmmType,volType,corrType);
	    swpn->testPrice(8192);
	    watch.stop(); watch.report(" ");

		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}
} 



void 
testCapletPrice(int lmmType, int volType, int corrType)
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
		printStars();
	    cout << "\n\nPricing caplet at i=n/3."
	         << "\nEnter dimension n of Libor process: n = ";
	    int n; cin >> n;
	
    	Timer watch; watch.start();
	    LiborDerivative* cplt=Caplet::sample(n,lmmType,volType,corrType);
	    cplt->testPrice(8192);
	    watch.report(" ");
		
		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}
} 



void 
testCallOnBondPrice(int lmmType, int volType, int corrType)
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
		printStars();
	    cout << "\n\nPricing at the money call on bond along [T_p,T_q] "
		     << "with random coupons c_j in [0.5,1.5]"
	         << "\n\nEnter p = ";
	    int p; cin >> p;
		cout << "Enter q = ";
	    int q; cin >> q;
		
	    Timer watch; watch.start();
	    BondCall* bc=BondCall::sample(p,q,lmmType,volType,corrType);
	    bc->testPrice(8192);
	    watch.report(" ");

		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}		
} 



void
testCallOnZeroCouponBondPrice(int lmmType, int volType, int corrType)
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
		printStars();
	    cout << "\n\nPricing at the money call expiring at T_{p-1}"
		     << "\non zero coupon bond maturing at T_p"
		     << "\nin LMM of dimension p+3."
	         << "\n\nEnter p = ";
	    int p; cin >> p;
		
	    Timer watch; watch.start();		
	    BondCall*  bc=BondCall::sampleCallOnZeroCouponBond(p,lmmType,volType,corrType);
	    bc->testPrice(8192);
	    watch.report(" ");

		cout << "\n\nDo another run (yes = 1, no = 0) do_again = ";
		cin >> do_again;
	}
} 



void
testLmmLattice()
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
		printStars();
	    cout << "\n\nBuilding and testing a lattice for the driftless libor market model:"
	         << "\n\nEnter dimension n of Libor process: n = ";
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



void
testSwaption()
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



void
testBermudanSwaption()
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



void
testCallOnBond()
{
    cout << "\nBond with random coupons c_j in [0.5,1.5] on [T_p,T_q]."
		 << "\n\nEnter p = ";
	int p; cin >> p;
    cout << "Enter q = ";
    int q; cin >> q;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    BondCall* bondcall = BondCall::sample(p,q); 
	bondcall->testPrice(nPath);
}



void
testCallOnZeroCouponBond()
{
    cout << "\nZero coupon bond matures at T_p:"
		 << "\n\nEnter p = ";
	int p; cin >> p;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    BondCall* bondcall = BondCall::sampleCallOnZeroCouponBond(p);
	bondcall->testPrice(nPath);
}



void
testCaplet()
{
    cout << "\nCaplet on [T_i,T_{i+1}], i=n/3 in LMM of dimension n."
	     << "\n\nEnter n = ";
	int n; cin >> n;
	cout << "Enter number of Libor paths: ";
    int nPath; cin >> nPath;
		
    Caplet* cplt = Caplet::sample(n); 
	cplt->testPrice(nPath);
}



void
testLiborDerivative()
{
	int do_again=1;
	// main loop
	while(do_again==1){
	
	    printStars();
		cout << "\n\nPricing of at the money Libor derivatives"
		     << "\nin the default driftless Libor Market Model."
		     << "\nChoose derivative:"
		     << "\nCaplet (0)"
		     << "\nSwaption (1)"
		     << "\nBermudan swaption (2)"
		     << "\nCall on bond (3)"
		     << "\nCall on zero coupon bond (4)" 
		     << "\n\nDerivative = ";
		     
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
