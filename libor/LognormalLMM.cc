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





#include "LognormalLMM.h"
using namespace Martingale;


/*******************************************************************************
 *
 *                                 LognormalLMM
 *
 ******************************************************************************/ 


	 
	 
	 // (X_p(T_t),...,X_{n-1}(T_t))\f$, 
     vector<Real>& LognormalLMM::XLvect(int t, int p) 
	 { 
		 XVec.setDimension(n-p);
		 XVec.setIndexBase(p);
		 UTRMatrix<Real>& eBtInv=eBInverse.getMatrix(t);
         // partial matrix-vector product eBtInv*eBY(t,*)
		 for(int j=p;j<n;j++){
			 
			 Real sum=0;
			 for(int k=j;k<n;k++) sum+=eBtInv(j,k)*eBY(t,k);
		     XVec[j]=exp(sum);
		 }
		 return XVec;
	 }
		 

	
//  CONSTRUCTOR
    
    LognormalLMM::
    LognormalLMM(LiborFactorLoading* fl) : LiborMarketModel(fl),
    Z(n), eBY(n), drifts(n),
	V(new Real[n]),
	eBInverse(n),
	covariationMatrixRoots(n-1),
	SG(new MonteCarloLiborDriver(n)),
	XVec(*(new vector<Real>(n)))
	{        
        // initialize eBY path array
        for(int j=0;j<n;j++) eBY(0,j)=log(x[j]);
			
		// Computation of the deterministic drift steps and 
		// volatility step covariance matrices
		
		// loop over accrual intervals [T_{t-1},T_t]
		for(int t=1;t<n;t++) { 
			
			Real T_t=Tc[t];
		    UTRMatrix<Real>& eBtInv=factorLoading->eBInverse(T_t,t);
			eBInverse.setMatrix(t,&eBtInv);
						
		    // It=integral_{T_t}^{T_{t+1}}eBu(s)ds
		    vector<Real> It(n-t,t); 

	        // Jt=integral_{T_t}^{T_{t+1}}eBnu(s)eBnu(s)'ds
		    UTRMatrix<Real> Jt(n-t,t); 

				
			// trapezoid rule with 4 points for the integral on [T_{t-1},T_t]
			Real s=Tc[t-1], dt=delta[t-1]/3;

			for(int m=0;m<4;m++){
				
				 // I(T_t)
			     vector<Real>& u=factorLoading->u(s,t);           // u(s)		
				 UTRMatrix<Real>& eBs=factorLoading->eB(s,t);     // eB(s)
				 u*=eBs;
				 u*=dt;
				 if((m==0)||(m==3)) u*=0.5;
				 It+=u;

				 UTRMatrix<Real>& nus=factorLoading->nu(s,t);     // nu(s)
				 eBs*=nus;                                        // eBnu(s)
				 UTRMatrix<Real>& eet=eBs.aat();                  // eBnu(s)eBnu(s)'
				 eet*=dt;
				 if((m==0)||(m==3)) eet*=0.5;
				 Jt+=eet;
				 
				 s+=dt;
			 }

			 
			// write the drift steps (drift vector is row of drifts matrix)
			for(int j=t;j<n;j++) drifts(t,j)=-It[j];     
			
            // Jt is the covariance matrix of the volatility step T_{t-1}->T_t
			// compute the upper triangular root for simulation
			covariationMatrixRoots.setMatrix(t-1,&(Jt.utrRoot()));
				
		} // end for t

		
	} // end constructor			
			
    
	
	LiborMarketModel* 
	LognormalLMM::sample(int n, Real delta)
    {
		LiborFactorLoading* 
		fl=LiborMarketModel::sampleFactorLoading(n,delta,LiborMarketModel::CS);
		return new LognormalLMM(fl);
	}

    

// WIENER INCREMENTS	
	
	
	void LognormalLMM::
	printWienerIncrements(int t, int T)
    {
         // recall that Z is an upper triangular array
		 for(int s=t;s<T;s++){
			 
             for(int k=s+1;k<n;k++){ cout << Z(s,k) << " "; 
				                     if(k==n-1) cout << endl; }
		}
	} // end printWienerIncrements
	

		
		
// TIME STEP

     // time step T_t -> T_{t+1} for Libors X_j, j>=p.
     void LognormalLMM::timeStep(int t, int p)
     {   
		 /* The matrix needed for the volatility step T_t->T_{t+1}.
          */
         UTRMatrix<Real>& R=covariationMatrixRoots.getMatrix(t); 
                    
         int q=max(t+1,p);   
         // only Libors L_j, j>=q make the step. To check the index shifts 
         // to zero based array indices below consider the case p=t+1 
         // (all Libors), then q=t+1.
         
         
         // volatility step vector V
         for(int j=q;j<n;j++)
         {
             V[j]=0;
             for(int k=j;k<n;k++) V[j]+=R(j,k)*Z(t,k);  
         }
 

         // drift steps are deterministic and precomputed
         for(int j=q;j<n;j++) eBY(t+1,j)=eBY(t,j)+drifts(t,j)+V[j];
                            
      
   }  // end timeStep
     

     	 
// SWAPTION AND CAPLET AGGREGATE VOLATILITIES (SIGMA)
	 

     Real LognormalLMM::
     capletAggregateVolatility(int i)
     { 
         Real volsqr=factorLoading->integral_sgi_sgj_rhoij(i,i,0,Tc[i]);
		 return sqrt(volsqr);
     } 
	 
	 

     Real LognormalLMM::
	 swaptionAggregateVolatility(int p, int q, Real T)
     { 
          UTRMatrix<Real>& 
		  Q=factorLoading->logLiborCovariationMatrix(p,q,0,T).utrRoot();
		  vector<Real> x_pq(q-p,p);
		  for(int j=p;j<q;j++) x_pq[j]=(B0(j)-B0(j+1))/B_pq(p,q);
		  x_pq*=Q.transpose();
		  return x_pq.norm();
     } 
    
     
     
     
// STRING MESSAGE
    
    string LognormalLMM::toString()
    {
		vector<Real> vols(n);
		for(int i=0;i<n;i++) vols[i]=vol(i); 

		ostringstream os;
		os << "\nLibor Market Model, random dynamics: " << SG->toString() << endl
		   << "fast predictor corrector simulation of true Libor dynamics." << endl
		   << factorLoading->toString() 
		   << "\n\nLibor volatilities:\n" << vols;
        
        return os.str();
     }
	 

	 
	 
// TEST PROGRAM
	 
void LognormalLMM::test(int n)
{ 
     Timer watch;
	
	 watch.start();
	 LiborFactorLoading* fl = CS_FactorLoading::sample(n);
	 LognormalLMM* lmm=new LognormalLMM(fl); 
	 
	 watch.report("Setup");
	
     for(int c=0;c<20;c++){
		 lmm->newPath();
	     cout << endl << endl << lmm->XLvect(n/2,n/2);
	 }
} // end test
	
	 

