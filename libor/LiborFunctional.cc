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



#include "LiborFunctional.h"


MTGL_BEGIN_NAMESPACE(Martingale)

	
	
Real 
LiborFunctional::
H_pq(int p, int q, const RealArray1D& H, const RealArray1D& delta)
{
	Real sum=0.0;
	for(int j=p;j<q;j++) sum+=delta[j]*H[j+1];
			
	return sum;
}


Real 
LiborFunctional::
S_pq(int p, int q, const RealArray1D& H, const RealArray1D& delta)
{
	Real H_pq=0.0;
	for(int j=p;j<q;j++) H_pq+=delta[j]*H[j+1];
			
	return (H[p]-H[q])/H_pq;
}



MTGL_END_NAMESPACE(Martingale)
