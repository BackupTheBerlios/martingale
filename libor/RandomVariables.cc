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


// the entire implementation is in the header


#include "RandomVariables.h"
using namespace Martingale;


/*******************************************************************************    
    
                        Standard Normal Vector 
						
*******************************************************************************/


	vector<Real> StandardNormalVector::nextValue() 
	{
	     int d=getDimension();
		 vector<Real> X(d);
	     for(int i=0;i<d;i++) X[i]=Random::sTN();
         return X;
    }



 

/*******************************************************************************    
    
                      Standard Normal Variable             
	
*******************************************************************************/

	// sample - control variate pair
	vector<Real> StandardNormalVariable::nextValue() 
	{
		vector<Real> v(2);
		v[0]=Random::sTN(); v[1]=v[0];
		return v;
    }


/*******************************************************************************    
    
                      Empirical Random Variable             
	
*******************************************************************************/

 
 
	void EmpiricalRandomVariable::test(int n)
    {
		Real* data_set=new Real[n];
		for(int i=0;i<n;i++) data_set[i]=i;
			
		EmpiricalRandomVariable X(data_set,n);
		for(int j=0;j<200;j++) cout << X.nextValue() << ", ";
	}
   
