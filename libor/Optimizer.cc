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



#include "Optimizer.h"
using namespace Martingale;

   
/*******************************************************************************
 *
 *                     		DOWNHILL SIMPLEX
 *
 ******************************************************************************/

 
// CONSTRUCTOR
    
    DownhillSimplex::
	DownhillSimplex(int n, RealArray1D& x, Real delta, int steps, bool _verbose) :
	Optimizer(n),
	vertices(n+1,n), sum(n), newVertex(n), y(n+1),
	fx(0.0), min(0), max(0), max2(0), nStep(steps),
	verbose(_verbose)
    {
        // vertices[n]=x
		for(int j=0;j<n;j++) vertices(n,j)=x[j]; 
		// vertices[j]=x+e_j, j<n
        for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){   
			
			vertices(i,j)=x[j];
            if(j==i)vertices(i,j)+=delta; 
        }
      
		// sum of all vertices
        for(int j=0;j<n;j++) sum[j]=(n+1)*x[j]+delta;       
            
    } // end constructor
	
	
	void DownhillSimplex::setInitialConditions()
    {
	    // function values at the vertices
        for(int i=0;i<n+1;i++) y[i]=f(vertices[i]);
        setHiLowIndices();
	}
              

	 
// REFLECTING, CONTRACTING, KEEPING TRACK OF HIGH-LOW POINTS	 
    

    void DownhillSimplex::
	setHiLowIndices()
    {
                    
            if(y[0]<y[1])
            { min=0; max2=0; max=1; }
            else 
            { min=1; max2=1; max=0; }
            
            for(int i=2;i<n+1;i++){ 
            
               fx=y[i];
               if(fx<y[min])min=i;                             // lowest
               else if(fx>y[max]){ max2=max; max=i; }          // highest
               else if(fx>y[max2])max2=i;                      // second highest
            }
    } // setHiLowIndices

    
    	
    Real DownhillSimplex::
	reflectVertex(int i, Real k)
    {       
        Real k1=(1+k)/n, k2=k1+k, w;
        for(int j=0;j<n;j++)newVertex[j]=k1*sum[j]-k2*vertices(i,j);
        w=f(newVertex);
        if(w<y[max])replaceWorst(w);
        return w;
    }   

	
    void DownhillSimplex::
	contract(int i)
    {
        for(int k=0;k<n+1;k++)
        for(int j=0;j<n;j++)
        vertices(k,j)=0.5*(vertices(i,j)+vertices(k,j));
        
        //update function values
        for(int j=0;j<n;j++)if(j!=i)y[j]=f(vertices[j]);
        setHiLowIndices();

        // update the barycenter
        for(int j=0;j<n;j++)sum[j]=0.5*((n+1)*vertices(i,j)+sum[j]);
        
    }  // end contract
	
    

    void DownhillSimplex::
	replaceWorst(Real w)
    {
        y[max]=w;
        // update center
        for(int j=0;j<n;j++) sum[j]=sum[j]-vertices(max,j)+newVertex[j]; 
		// copy newVertex into vertices[max]
		for(int j=0;j<n;j++) vertices(max,j)=newVertex[j];
        setHiLowIndices();
    }
	
	
		
// SEARCH FOR THE MINIMUM 

	
    RealArray1D& DownhillSimplex::search()
    {
        setInitialConditions();
		for(int step=0;step<nStep;step++){

            // min, max2, max are set, see constructor
            
            // reflect the worst vertex through the center of the
            // opposite face
            fx=reflectVertex(max,1.0);

            if(fx<y[min]) // extend in this direction by a factor of two
            { 
               fx=reflectVertex(max,-2.0);
            }
            else if(fx>y[max2])
            {  
               fx=reflectVertex(max,-0.5);
               if(fx>=y[max]) contract(min);
            }
        }
        
        if(verbose){
            
            // the state before the last step
            cout << "\n\nDownHillSimplex::search(): function minimum: " << y[min]
			     << "\nMinimizing vector:";
            for(int j=0;j<n;j++) cout << endl << vertices(min,j);
        }
		
		RealArray1D& xvec=*(new RealArray1D(n));
		for(int j=0;j<n;j++) xvec[j]=vertices(min,j);
        return xvec;
        
    } // end search

                   
 

/*******************************************************************************
 *
 *                        BFGS
 *
 ******************************************************************************/

    
	// definition of the static members
    Real BFGS::EPSX=0.00000001;
    Real BFGS::KF=0.0005;
    Real BFGS::EPSG=0.00000001;
    Real BFGS::EPS=0.000000001;



// CONSTRUCTOR


    BFGS::
	BFGS(int n, RealArray1D& _x, int _nVals, Real _stepmax, 
	     RealArray1D& _h, int _nRestarts, bool _verbose=false) :
    Optimizer(n),
	stepmax(_stepmax),
	nVals(_nVals), fVals(0), restarts(0), nRestarts(_nRestarts),
	x0(n),
	grad0(n),
	x(_x),
	xDelta(n),
	grad(n),
	z(n),
	d(n),
	gDelta(n),
	hdg(n),
	u(n),
	h(_h),
	fx0(0.0), fx(0.0),
	H(n,n),
	verbose(_verbose)
    {
        // approximate inverse Hessian intialized as identity matrix
        for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){ H[i][j]=0; if(j==i)H[i][j]=1.0; }
        
    }
	
	
           
// MIN SEARCH                

	
     RealArray1D& BFGS::search()
     {
		 setInitialConditions();
		 if(fVals==nVals/nRestarts) resetHessian();
         // main loop to termination
         while(true){
			            
            // line search from current point x in the current direction d
            lineSearch();
           
            // hessian update and next search direction
            bfgsUpdate();
            for(int i=0;i<n;i++){
                
                d[i]=0.0;
                for(int j=0;j<n;j++)d[i]-=H[i][j]*grad[j];
            }
            
            // CHECK TERMINATION 
			
            // returning from line search after using up all
            // function evaluations?
            if(fVals>nVals){
                
				if(verbose)
                cout << "\n\nBFGS.search(): number of function evaluations exhausted.";
                break;
            }

            
            // termination on zero gradient
            if(gradientIsZero()){
                
                if(verbose)
                cout << "\n\nBFGS.search(): zero gradient, fVals=" << fVals;
                break;
            }
            
            
            // termination on zero gradient
            if(relativeSize(xDelta,x)<EPSX){
				
                if(restarts<nRestarts){ 
					
					if(verbose)cout << "\n\nBFGS.search(): restarting";
					resetHessian(); restarts++; 
				}else{
				
                    if(verbose)cout << "\n\nBFGS.search(): no movement in x, fVals=" << fVals;
                    break;
				} // end if
            } // endif

         } // end main loop
		 
		// report results
		if(verbose){
            
            // the state before the last step
            cout << "\n\nBFGS::search(): function minimum: " << fx
			     << "\nMinimizing vector:";
            for(int j=0;j<n;j++) cout << endl << x[j];
        }
		
        return x;

   } // end search

	
	

    void BFGS::setInitialConditions()
    {            
        fx=f(x);   fVals++;
        grad=gradcdF(x);
    
        // initial line direction
        for(int j=0;j<n;j++)d[j]=-grad[j];

    } // end setInitialConditions

    
	Real BFGS::relativeSize(RealArray1D& x, RealArray1D& y)
    {
        Real rs=0, temp;
        for(int j=0;j<n;j++){ 
            temp=abs(x[j])/(1+abs(y[j]));
            if(temp>rs)rs=temp;
        } 
        return rs;
    }
    
   
    bool BFGS::gradientIsZero()
    {   
         Real q=0.0, tmp;
         for(int j=0;j<n;j++){
             
             tmp=abs(grad[j])*(1+abs(x[j]))/(1.0+fx);
             if(tmp>q)q=tmp;
         }
         return (q<EPSG);
    }
	
	
	void BFGS::resetHessian()
    {
		for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){ H[i][j]=0.0; if(i==j)H[i][j]=1.0; }
	}
    

// BFGS UPDATE
     
    void BFGS::bfgsUpdate()
    {            
         // (x-x0).(grad-grad0)
         Real dxdg=xDelta.dotProduct(gDelta);
            
         // skip update if dxdg is too small
         if(abs(dxdg)<sqrt(EPS*xDelta.norm()*gDelta.norm())) return;
                        
         // (grad-grad0)*H*(grad-grad0)
         Real dgHdg=gDelta.dotProduct(hdg);
            
         // u_j coordinates
         for(int j=0;j<n;j++)u[j]=xDelta[j]/dxdg-hdg[j]/dgHdg;
            
         // finally the update, H is symmetric
         for(int i=0;i<n;i++)
         for(int j=0;j<=i;j++){
                
             H[i][j]+=(xDelta[i]*xDelta[j]/dxdg)+
                      (hdg[i]*hdg[j]/dgHdg)+dgHdg*u[i]*u[j];
             H[j][i]=H[i][j];
         }
    } // end bfgsUpdate

	
	
// GRADIENT COMPUTATION  

	
	// forward difference
    RealArray1D& BFGS::gradF(RealArray1D& x, Real fx)
    {
        for(int j=0;j<n;j++){
            
			// u is used as workspace here
            for(int i=0;i<n;i++) u[i]=x[i];
            u[j]+=h[j];
            Real dxj=x[j]-u[j],
                   dyj=fx-f(u);
            
            fVals++; 
            grad[j]=dyj/dxj;
        }
		
        return grad;
    } // end gradF
	
            
    // central difference
    RealArray1D& BFGS::gradcdF(RealArray1D& x)
    {
        Real fx1,fx2;
        
        for(int j=0;j<n;j++){
            
            for(int i=0;i<n;i++)u[i]=x[i];       // u=x
            
            u[j]+=h[j];                          // u=x+h*e_j
            fx1=f(u);
            fVals++; 
            
            u[j]-=2*h[j];                        // u=x-h*e_j
            fx2=f(u);
            fVals++; 
            
            grad[j]=(fx1-fx2)/(2*h[j]);
        } // end j
		
        return grad;
    } // end gradcdF

	
	
// BACK TRACKING DURING LINE SEARCH                 

	Real BFGS::
	backTrack(Real t1, Real t2, Real f1, Real f2, Real k)
    {
         Real u,v,a,b,d,t;
         u=f1-fx0-t1*k;
         v=f2-fx0-t2*k;
         a=(u/(t1*t1)-v/(t2*t2))/(t1-t2);
         b=(-t2*u/(t1*t1)+t1*v/(t2*t2))/(t1-t2);
                    
         if(a==0.0) t=-k/(2.0*b);
         else
         {
              d=b*b-3*a*k;
              if(d<0.0)          t=0.5*t1;
              else if(b<=0)      t=(-b+sqrt(d))/(3.0*a);
              else               t=-k/(b+sqrt(d));
         } // end t
                    
         // ensure t<=t1/2
         if(t>0.5*t1) t=0.5*t1;
         return t;
    }
	
	
//  LINE SEARCH
    
	void BFGS::lineSearch()
    { 
        Real dNorm=d.norm(), f1=0.0, f2=0.0,
		     // relative size of directional vector to current point
		     rs=relativeSize(d,x);        
         
        // scale directional vector d to a maximum length of
        // stepmax, then do backward line search.
        // proper BFGS does not do that, only scales down to stepmax!!!!
        if(dNorm>stepmax)d.scale(stepmax/dNorm);
        
        Real k=grad.dotProduct(d);   // slope     
        if(k>=0) cout << "\n\nBFGS.lineSearch(): no descent, roundoff errors.";
        
        Real t_min=EPSX/rs,     // terminate if t<t_min
             t1=1.0,            // current line parameter t 
             t2=1.0,            // old line parameter t, lags by one step
             t=0.0;             // next line parameter t in the backtracking sequence
		
		int backtracklevel=1;   // number of current back track step
		
        // loop till termination
        bool done=false;
        while(!done){
            
            // step from x in direction d, 
            // write new point into z.
            for(int j=0;j<n;j++)z[j]=x[j]+t1*d[j];
            
            f1=f(z);              // f(x+t1*d)
            fVals++;              // count function evaluations
            
            // sufficient decrease in the value of f(x):
            if(f1<fx+KF*t1*k) break;
            
			if(backtracklevel==1)  // first backtrack: square approximation
            {   t=-k/(2.0*(f1-fx-k)); backtracklevel++; }
            else                   // subsequent backtracking
            {   t=backTrack(t1,t2,f1,f2,k); backtracklevel++; }
           
            
            // move the backtracking state forward
            t2=t1; f2=f1; 
            if(t<0.1*t1) t1=0.1*t1; else t1=t;     // prevent collapse of step size       
            if((fVals>nVals)||(t1<t_min)) done=true;

       } // end main loop
        
       // move the state forward to the new point x=z
       x0=x; fx0=fx;
       x=z;  fx=f1;
       grad0=grad;                  // save the old gradient
       grad=gradcdF(x);             // get the new gradient
	   
       for(int j=0;j<n;j++){
                
          xDelta[j]=x[j]-x0[j];
          gDelta[j]=grad[j]-grad0[j];
       }
       
       // hdg=H*dg
       for(int i=0;i<n;i++){
                
          hdg[i]=0.0;
          for(int j=0;j<n;j++)hdg[i]+=H[i][j]*gDelta[j];
       }
	   
    } // end lineSearch
	
	
   
/*******************************************************************************
 *
 *                   SOBOL SEARCH
 *
 ******************************************************************************/
	
	RealArray1D& SobolSearch::
	search()
    {		
		   if(verbose)
            cout << "\n\nSobolSearch::search(): dimension = " << n
		         << ", nPoints = " << nPoints;
		
		Real fx, fOpt=f(x);
		RealArray1D& xOpt=*(new RealArray1D(x));
		
		while(nPoints>0){
			
			for(int k=0;k<nPoints/2;k++){
				
			    // next quasi random point in the search window 
			    // intersected with the search domain.
		        Real* u=lds->nextPoint();
		        for(int i=0;i<n;i++) x[i]=xOpt[i]+d[i]*(1-2*u[i]);
			
		        while(!isInDomain(x)){
			
			         u=lds->nextPoint();
		             for(int i=0;i<n;i++) x[i]=xOpt[i]+d[i]*(1-2*u[i]);
		        }
			
				fx=f(x);
			    if(fx<fOpt){
					
					fOpt=fx; xOpt=x;
					if(verbose) cout << "\nmin = " << fx;
				}
			} // end for k
			nPoints/=2;
			for(int i=0;i<n;i++) d[i]*=q;
			lds->restart();
				
		} // end while
		
		// report results
		if(verbose){
            
            // the state before the last step
            cout << "\n\nFunction minimum: " << fx
			     << "\nMinimizing vector:";
            for(int j=0;j<n;j++) cout << endl << x[j];
        }
		
		return xOpt;
		
	} // end search
		
		    
	

	
	
// TEST FUNCTIONS
	
	   
    // the objective function in the example
	Real ObjectiveFunction::function_1(int n, Real* x)
    {
         Real u=0.0;
         for(int i=0;i<n;i++)u+=exp(-x[i]*x[i]);
         return 1/u;
     }
              


