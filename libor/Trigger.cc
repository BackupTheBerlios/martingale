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
spyqqqdia@yahoo.com

*/


#include "Trigger.h"
#include "BermudanOption.h"
#include "LiborMarketModel.h"
#include "Optimizer.h"


MTGL_BEGIN_NAMESPACE(Martingale)

using std::cout;
using std::endl;



/*******************************************************************************
 *
 *                     PjTrigger
 *
 ******************************************************************************/


PjTrigger::
PjTrigger
(LiborMarketModel* lmm, BermudanSwaption* swpn, bool vbose) :
LMM(lmm),
bswpn(swpn),
kappa(swpn->kappa),
nPath(swpn->nPath),
p(swpn->p),
q(swpn->q),
verbose(vbose),
path(swpn->nPath),
p1(q-p,p),
p2(q-p,p),
p3(q-p,p),
S(q-p,p)
{
	if(verbose)
	cout << "\n\nExercise trigger:";
	
	// index base t=p along time t
	for(int i=0;i<nPath;i++) path[i]=new RealArray2D(q-p,4,p,0);
	fillTrainingPaths();
	computeCoefficients();
}


PjTrigger::
~PjTrigger()
{
	for(int i=0;i<nPath;i++) delete path[i];
}


           
bool 
PjTrigger::
isTriggered(int t, int s)
{
	if(s<p)return false;
    if(s==q)return true;
          
    double f=LMM->L(s,s);
    if(s==q-1)return(f>kappa);
          
    double sw=LMM->swapRate(s+1,q,s),
           h=bswpn->currentForwardPayoff(s);
          
    return ((f-p3[s])*(sw-p2[s])>p1[s]*S[s])&&(sw>p2[s])&&(h>0);
}


int 
PjTrigger::
nextTriggered(int t)
{
	// Libors needed only up to time T_q
	LMM->newWienerIncrements(q);
	int s=t;
	while(!isTriggered(t,s)){ LMM->timeStep(s,p); s++; }
	
	return s;
}


bool 
PjTrigger::
exercise(int i, int t, const RealArray1D& x)
{
    RealArray2D& path_i=getTrainingPath(i);
    Real f=path_i[t][0],
         sw=path_i[t][1];
               
    return ((f-x[2])*(sw-x[1])>x[0]*S[t])&&(sw>x[1])&&(path_i[t][2]>0);
}
          
 

void 
PjTrigger::
fillTrainingPaths()
{
    if(verbose)
	cout << "\nComputing training paths.";
	
	for(int i=0;i<nPath;i++){
		
        RealArray2D& path_i=getTrainingPath(i);    
        LMM->newPath(q-1,p);
        for(int t=p;t<q-1;t++){  
                
			path_i[t][0]=LMM->L(t,t);
            path_i[t][1]=LMM->swapRate(t+1,q,t);
            path_i[t][2]=bswpn->currentForwardPayoff(t);
        } 
        // however for t=q-1 (no next swap rate) and we must set
        // the optimal exercise time s=t or s=t+1 (no exercise)
        int t=q-1;
        double f=LMM->L(t,t);
        path_i[t][0]=f;
        path_i[t][2]=bswpn->currentForwardPayoff(t);
        if(f>kappa) path_i[t][3]=t;
        else path_i[t][3]=t+1;
    } // end i
} // end fillPaths
 


Real 
PjTrigger::
objectiveFunction(int t, const RealArray1D& v)
{
    Real sum=0.0;
    for(int i=0;i<nPath;i++){
        
        RealArray2D& path_i=getTrainingPath(i);
        // exercise time along this path, note t<=q-2
        int s;
        if(exercise(i,t,v)) s=t;
        else s=(int)(path_i[t+1][3]); 
                        
        // s=n never exercised, otherwise payoff is in path_i[s][2]
        if(s<q) sum+=path_i[s][2];
    } // end i
//cerr << "\n\nf() = " << -sum;
    return -sum;
} 



              
void 
PjTrigger::
computeCoefficients()
{
    if(verbose)
	cout << "\nComputing coefficents parametrizing the exercise boundary.";
	    
    // common to all optimizers
    int nVals=4000;
    // initial search point, 
    // dynamically adjusted to be the optimum of the last search
    RealArray1D z(3), h(3);
	z[0]=2*kappa/3; z[1]=kappa; z[2]=0.75*kappa;
    h[0]=z[0]/16; h[1]=z[1]/16; h[2]=z[2]/16;
    Real stepmax=kappa/4;

    // backward recursive computation of the coefficient p_j(t)
    for(int t=q-2;t>=p;t--){
            
        if(verbose) cout << "\nExercise point t = " << t;
            
        LocalBFGS optimizer(this,t,z,nVals,stepmax,h,verbose);
		// very slow but thorough:
		// LocalSobolSearch optimizer(this,t,z,nVals,h,verbose);   
        const RealArray1D& pj=optimizer.search();
        p1[t]=pj[0];
        p2[t]=pj[1];
        p3[t]=pj[2];
            
        // for each path set the optimal exercise time s>=t starting
        // at time t
        for(int i=0;i<nPath;i++){
			
            RealArray2D& path_i=getTrainingPath(i);
            if(exercise(i,t,pj)) path_i[t][3]=t;
            else path_i[t][3]=path_i[t+1][3];
        }
            
        //reset starting point to last optimum
        z[0]=pj[0]; z[1]=pj[1]; z[2]=pj[2];
        h[0]=z[0]/4; h[1]=z[1]/4; h[2]=z[2]/4;
        stepmax=max(min(z[0],min(z[1],z[2])),0.01)/2;
      
    }// end t  
} // end computeCoefficients
 

// PARAMETER OPTIMIZATION AT EACH EXERCISE TIME t  

PjTrigger::LocalBFGS::
LocalBFGS
(PjTrigger* pj_trg, int s, const RealArray1D x, int nVals, 
 Real stepmax, const RealArray1D h, bool verbose
) :
// dimension 3 (p1,p2,p3) and 0 restarts of the optimizer
BFGS(x,nVals,stepmax,h,0,verbose), 
t(s), trg(pj_trg)
{  } 
       

Real 
PjTrigger::LocalBFGS::
f(const RealArray1D& v){ return trg->objectiveFunction(t,v); }



PjTrigger::LocalSobolSearch::
LocalSobolSearch
(PjTrigger* pj_trg, int s, const RealArray1D x, int nPoints, 
 const RealArray1D delta, bool verbose
) :
// dimension 3 (p1,p2,p3) 
SobolSearch(x,nPoints,delta,verbose), 
t(s), trg(pj_trg)
{  } 
       

Real 
PjTrigger::LocalSobolSearch::
f(const RealArray1D& v){ return trg->objectiveFunction(t,v); }



		
MTGL_END_NAMESPACE(Martingale)

