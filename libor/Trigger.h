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


#ifndef martingale_trigger_h    
#define martingale_trigger_h

#include "TypedefsMacros.h"
#include "Optimizer.h"        // base class BFGS


MTGL_BEGIN_NAMESPACE(Martingale)


// dependencies
class LiborMarketModel;
class BermudanSwaption;




/*! \file Trigger.h
 *  A trigger represents a family of stopping times \f$\tau_t\geq t\f$.
 *  This is needed for Bermudan option pricing. See book, 4.6.
 *  In this context the trigger is used to trigger option exercise.
 */
 
 

/**********************************************************************************
 *
 *                                     TRIGGER
 *
 *********************************************************************************/ 
 
/** Let T denote discrete time (number of time steps) to the horizon.
 *  The class trigger implements a family of integer valued stopping times
 *  \f[\tau_t:\Omega\to[t,T].\f]
 *  The equality \f$\tau_t=s\f$ means that s is the first time \f$s>t\f$ 
 *  at which a certain event occurs.
 */
class Trigger {
	
public:
	
/** Returns true if \f$s>t\f$ and the trigger event occurs at time s.
 *  In practice s is usually the first time after time t at which
 *  the trigger event occurs.
 */
virtual bool isTriggered(int t, int s) = 0;
	

}; // end Trigger



/*******************************************************************************
 *
 *                     PjTrigger
 *
 ******************************************************************************/


/** <p>Exercise trigger of a Bermudan swaption using the parametrization of the
 *  exercise boundary given by P. Jaeckel, "Monte Carlo Methods in Finance", 12.7. 
 */
class PjTrigger : public Trigger {
	
    
    LiborMarketModel* LMM;
    BermudanSwaption* bswpn;
   
    double kappa;       // swaption strike
   
    int nPath,          // number of training paths
        p,              // swaption exercise begins T_p
        q;              // swap ends T_q 
   
    bool verbose;       // messages during optimization
   
    Array1D<RealArray2D*> path; 
	// the components 0,1,2,3 of path_i=*path[i] are as follows:   
    // path_i[t][0].... Libor X(t,t) 
    // path_i[t][1].... swaprate S_{t+1,n}(t) 
    // path_i[t][2].... forward swap payoff from immediate exercise (>=0)
    // path_i[t][3].... time s>=t of optimal exercise starting at time t
    // all computed from Libor sample path i for p<=t<q
   
    RealArray1D p1,p2,p3;  // the coefficients p_j(t) in Jaeckel's trigger
	                       // at each exercise time t.
    RealArray1D S;         // S[t]=S_{t+1,n}(0), swaprate
	
public:
	
// ACCESSORS

/** The number of training paths.*/
int getPaths(){ return nPath; }
/** Swap period ends at \f$T_q\f$.*/
int get_q(){ return q; }
/** Array storing the information from training path i.*/
RealArray2D& getTrainingPath(int i){ return *(path[i]); }

	
// CONSTRUCTOR

/** @param lmm the Libor Market Model.
 *  @param swpn the Bermudan swaption.
 *  @param paths the number of training paths.
 *  @param vbose messages during initialization.
 */
PjTrigger
(LiborMarketModel* lmm, BermudanSwaption* swpn, bool vbose);

~PjTrigger();


           
// THE TRIGGER CONDITION

/** True if exercise is triggered at time <code>t</cod> false otherwise.
 *
 * @param t irrelevant.
 * @param s current time.
 */
bool isTriggered(int t, int s);


/** <p>True if the exercise condition is met under coefficients 
 *  <code>x</code> at time <code>t</code>  along training path 
 *  <code>i</code> , false otherwise. Exercises only if payoff>0. 
 */
bool exercise(int i, int t, const RealArray1D& x);
          
 
private:

/** Fills the training path arrays with all information necessary to
 *  decide exercise and compute the coefficients. 
 */
void fillTrainingPaths();
 
           
/** The pj-coefficients computed by the local optimizer.
 */      
void computeCoefficients();
 

// PARAMETER OPTIMIZATION AT EACH EXERCISE TIME t  

class LocalOptimizer : public BFGS {
           
	int t;                // time t at which the parameters are optimized
	PjTrigger* trg;       // the trigger which contains the optimizer
	
public:
       
       
   /** Optimizer for the parameters <code>p0,p1,p2</code>
    *  used to parametrize the exercise boundary.
	*
	* @param pj_trg the PjTrigger containing the optimizer.
    * @param s exercise time.
	* @param paths number of training paths.
	* @param x starting parameter vector.
    * @param nVals number of function evaluations.
    * @param stepmax maximum stepsize in line search.
    * @param h directional increments for gradient computation.
    * @param verbose messages during optimization.
    */
    LocalOptimizer
    (PjTrigger* pj_trg, int s, const RealArray1D x, int nVals, 
     Real stepmax, const RealArray1D h, bool verbose);
       
    // the function we want to optimize
	// forward payoff fom exercise for s>=t under parameters v
    Real f(const RealArray1D& v);

}; // end LocalOptimizer
 
};  // end PjTrigger

	

	
MTGL_END_NAMESPACE(Martingale)

#endif
 
