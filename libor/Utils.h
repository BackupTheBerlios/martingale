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

/*
 * Utils.h
 *
 * Created on February 3, 2003, 12:00 PM
 */


#include <string>
#include <sstream>
#include <iostream>
#include "TypedefsMacros.h"
#include "Matrices.h"


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
	
	void report(string message)
	{
	     int  totalSeconds=(int)(after-before),
	          minutes=totalSeconds/60,
	          seconds=totalSeconds%60;
       	 if(minutes==0)
	         std::cout<<"\n\n"<<message<<"\ntime: "<<seconds<<" seconds"<<endl;
	     else
	         std::cout<<"\n\n"<<message<<"\ntime: "<<minutes<<" minutes, "<<seconds<<" seconds"<<endl;
	} // end report
	
}; // end Timing


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



/** print vector
 */
template<class S>
std::ostream& operator << (std::ostream& os, const vector<S>& v)
{
	int d=v.getDimension();
    S* vdptr=v.getData();
	os << endl << "vector of dimension " << d << ":" << endl;
    for(int i=0;i<d-1;i++) os << vdptr[i] << ", ";
    os << vdptr[d-1];
    return os << endl << endl;
} // end operator <<



/** print lower triangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const LTRMatrix<S>& A)
{
	int dim=A.getDimension(), b=A.getIndexBase();
	os << endl << "Lower triangular matrix, dimension " << dim << ":" 
	   << endl << endl;
	for(int i=b;i<dim+b;i++){
		
		for(int j=b;j<i;j++) os << A(i,j) << ", ";
		os << A(i,i) << endl;
	}
    return os << endl << endl;
} // end operator <<



/** print upper triangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const UTRMatrix<S>& U)
{
	os << endl << "Transposed matrix:";
    return os << U.transpose() << endl << endl;
} // end operator <<


/** print rectangular matrix
 */
template<class S>
std::ostream& operator << (std::ostream& os, const Matrix<S>& A)
{
	int rows=A.getnRows(), 
	    cols=A.getnCols();
	S** D=A.getData();
	os << endl << "Rectangular " << rows << " by " << cols << " matrix:"
	   << endl << endl;
	for(int i=0;i<rows;i++){
		
	    for(int j=0;j<cols-1;j++) os << D[i][j] << ", ";
		os << D[i][cols-1] << endl;
	}
    return os << endl << endl;
} // end operator <<

	



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
void clear()
{		
	int m = message.length();
	for(int i=0;i<m;i++){ cerr<<"\b"; }
	for(int i=0;i<m;i++){ cerr<<" "; }
	for(int i=0;i<m;i++){ cerr<<"\b"; }
}
	

/** Reports current progress and projects time left from a loop
 * over N iterations when the current iteration is n. Progress report has to be 
 * called at iteration n=N/100 when it times the loop and uses this time to project 
 * time to completion every time it is called thereafter. Progress is reported to the console.
 *
 * @param n Current loop iteration.
 * @param N Total number of iterations in loop.
 */
void consoleReport(int n, int N)
{
	if(n<N/100) return;
	if(n==N/100) { 
   
            clocks_for_one_percent = clock() - clocks_at_start;
            clear();
    }
	
    // time to completion	
	double  percent_left = 100.0*(N-n)/N;
	long	secs_left    = (long)(percent_left*clocks_for_one_percent/CLOCKS_PER_SEC);
	
    if(secs_left<seconds_left-5){        //print a progress report 
		
		int
	        hrs     = secs_left/3600,
  	        mins    = (secs_left%3600)/60,
	        secs    = secs_left%60;
		// convert to string
        Int_t Hrs(hrs), Mins(mins), Secs(secs);
        
	    string
	        hours   = Hrs.toString(),
	    	minutes = Mins.toString(),
		    seconds = Secs.toString();
	
	    clear();
	
	    if(hrs>0) 
		    message = hours + " hr, " + minutes + " min, " + seconds + " sec.";
	    else if(mins>0)
		    message = minutes + " min, " + seconds + " sec.";
	    else 
		    message = seconds +" sec.";
	
	    std::cerr << message;
	    seconds_left=secs_left;
	} // endif
} // end consoleReport


}; //end LoopStatus




           
MTGL_END_NAMESPACE(Martingale)

#endif
