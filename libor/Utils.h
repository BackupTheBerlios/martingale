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

#ifndef martingale_utils_h
#define martingale_utils_h


#include <string>
#include <sstream>
#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)




/***************************************************************************************

                                  TIMING

***************************************************************************************/


class Timer {
	
	double before, 
	       after;   // system time in seconds
	
public:
	Timer(){ before=0; after=0; }
    void start(){ before=clock()/CLOCKS_PER_SEC; }
	void stop(){  after=clock()/CLOCKS_PER_SEC;  }
	
	/** Time to completion. */
	void report(string message);	
	
}; // end Timer


/***************************************************************************************

                STRING REPRESENTATION OF VARIOUS OBJECTS

***************************************************************************************/
 

/** encapsulation of a type with some useful utility functions
 *  (string conversion).
 */
template<class T>
class Type_t {

    T t;

public:
	
    Type_t(T t0): t(t0) { }
  
    string toString()
    {
         ostringstream os;
         os<<t;
         return os.str();
    }
	
}; // end Type_t



typedef Type_t<int> Int_t;
typedef Type_t<long> Long_t;
typedef Type_t<double> Double_t;
typedef Type_t<long double> LongDouble_t;
typedef Type_t<Real> Real_t;


/** General printing of all types which implement printSelf(std::ostream& os)
 */
template<class C>
std::ostream& operator << (std::ostream& os, const C& c)
{
	return c.printSelf(os);
}




/***************************************************************************************

                 Loop Status

***************************************************************************************/


/** Provides methods to report the progress of a loop and project
 *  time to completion.
 *
 * @author  Michael J. Meyer
 */
class LoopStatus {
        
    long  clocks_at_start,                 // clocks at beginning of loop
	      clocks_for_one_percent,          // clocks to complete 1% of the loop
	      seconds_left;
	
    string  title,                         // description of what's going on in the loop
	        message;                       // message printed
	
public:

	
/** Constructor
 *
 * @param loopTitle string descriptive of loop computation.
 */
LoopStatus(string loopTitle) :
//initialization list
clocks_at_start(clock()), clocks_for_one_percent(0),
seconds_left(1000000000), title(loopTitle),
message("timing the first 1% of the loop.")
{ std::cerr << "\n\n"+title+":\n"+message; }	


/** Clears last string from console
 */
void clear();
	

/** Reports current progress and projects time left from a loop
 * over N iterations when the current iteration is n. Progress report has to be 
 * called at iteration n=N/100 when it times the loop and uses this time to project 
 * time to completion every time it is called thereafter. Progress is reported to the console.
 *
 * @param n Current loop iteration.
 * @param N Total number of iterations in loop.
 */
void consoleReport(int n, int N);


}; //end LoopStatus




           
MTGL_END_NAMESPACE(Martingale)

#endif
