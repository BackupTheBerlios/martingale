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

#ifndef martingale_pathgenerator_h    
#define martingale_pathgenerator_h

#include "TypedefsMacros.h"
#include "LiborMarketModel.h"
#include "Trigger.h"


MTGL_BEGIN_NAMESPACE(Martingale)


/*! \file PathGenerator.h
 *  Controllers for the computation of paths (typically assets or interest 
 *  rates underlying options). The path generator sees to it that the minimum 
 *  number of variables is evolved only as far as needed.
 */



/** Path generators control the computation of asset price or interest rate 
 *  paths for the purpose of Monte Carlo option pricing. The path generator
 *  sees to it that the minimum number of variables is evolved only as far as 
 *  needed.
 */
class PathGenerator {
	
public:
	
	virtual void newPath() = 0;
		
}; // end PathGenerator


// LIBOR PATH GENERATORS


/** PathGenerator for Libor derivatives. Maintains a reference to a
 *  LiborMarketModel and forwards request for paths to the LMM.
 */
class LiborPathGenerator : public PathGenerator {
	
protected:
	
	LiborMarketModel* LMM;
	
public:
	
LiborPathGenerator(LiborMarketModel* lmm) : LMM(lmm) {  }
	
/** Discrete time t (continuous time T_t) until which the current
 *  Libor path was computed.
 */
virtual int getTime() = 0;

/** Sets the trigger for path generators which compute paths until a 
 *  stop is triggered. This is the default implementation: no action,
 *  appropriate for path generators computing paths to a fixed time.
 */
virtual void setTrigger(Trigger* trg){  }
	
};


/** PathGenerator evolving Libors \f$X_j, j\geq i\f$,
 *  until a fixed time t.
 */
class LiborPathsToFixedTime : public LiborPathGenerator {
	
	/** Libors are computed until time T_t.*/   int t;                     
	/** Libors evolved are L_j, j>=i.*/         int i;                     
	
public:
	
	/** @param lmm the underlying LiborMarketModel.
	 *  @param s time until which Libors are computed.
	 *  @param k Libors evolved are L_j, j>=k.
	 */
	LiborPathsToFixedTime(LiborMarketModel* lmm, int s, int k) :
	LiborPathGenerator(lmm), 
	t(s), i(k)
    {   }
	
	void newPath(){ LMM->newPath(t,i); }
	int getTime(){ return t; }
	
};


/** PathGenerator evolving Libors \f$X_j, j\geq i\f$,
 *  until a {@link Trigger} is triggered
 */
class LiborPathsToTriggerTime : public LiborPathGenerator {
	
    Trigger* trigger;  // trigger stops path computation
	int i;             // Libors evolved are L_j, j>=i. 
	int t;             // time at which path computation stops under all conditions.
	int s;             // time until which the last Libor path was computed.
	
public:
	
	/** Constructor initializes the trigger with null since typically
	 *  the entire derivative has to be initialized before the trigger
	 *  can be set. The trigger must be set calling 
	 *  <code>setTrigger()</code> as soon as everything is in place.
	 *
	 *  @param lmm the underlying LiborMarketModel.
	 *  @param r no Libor path computed past time \f$T_r\f$. 
	 *  @param k Libors evolved are L_j, j>=k.
	 */
	LiborPathsToTriggerTime(LiborMarketModel* lmm, int r, int k) :
	LiborPathGenerator(lmm), i(k), t(r)
    {   }
	
	void newPath();
	void setTrigger(Trigger* trg){ trigger=trg; }
	int getTime(){ return s; }
	
};


MTGL_END_NAMESPACE(Martingale)

#endif
 
