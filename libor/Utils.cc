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


#include "Utils.h"

MTGL_BEGIN_NAMESPACE(Martingale)



/***************************************************************************************

                                  TIMING

***************************************************************************************/

	
	void Timer::report(string message)
	{
	     int  totalSeconds=(int)(after-before),
	          minutes=totalSeconds/60,
	          seconds=totalSeconds%60;
       	 if(minutes==0)
	         std::cout<<"\n\n"<<message<<"\ntime: "<<seconds<<" seconds"<<endl;
	     else
	         std::cout<<"\n\n"<<message<<"\ntime: "<<minutes<<" minutes, "<<seconds<<" seconds"<<endl;
	} // end report





/***************************************************************************************

                 Loop Status

***************************************************************************************/


/** Clears last string from console
 */
void LoopStatus::clear()
{		
	int m = message.length();
	for(int i=0;i<m;i++){ cerr<<"\b"; }
	for(int i=0;i<m;i++){ cerr<<" "; }
	for(int i=0;i<m;i++){ cerr<<"\b"; }
}
	

void LoopStatus::consoleReport(int n, int N)
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



MTGL_END_NAMESPACE(Martingale)


