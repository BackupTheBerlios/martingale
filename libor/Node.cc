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

#include "LmmLattice.h"
using namespace Martingale;



/**********************************************************************************
 *
 *            GENERAL NODES 
 *
 *********************************************************************************/


// GENERAL NODES

	void Node::
	printTransitionProbabilities()
	{
		cout << "\n\n\nTransition probabilities\n: ";
		std::list<Edge>::const_iterator theEdge;
		for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge) 
		{	
			Real p=theEdge->probability;
			cout << p << ", ";
		}		
     } // end printTransitionProbabilities
	 


/**********************************************************************************
 *
 *            GENERAL LMM NODES 
 *
 *********************************************************************************/
	 
	 
		
	LmmNode::
	LmmNode(LiborFactorLoading* fl, int t) : Node(t),
	factorLoading(fl), 
	n(factorLoading->getDimension()),
	H(n-t+1,t) 
	{  }
	 

	Real LmmNode::
	Hpq(int p, int q)
    {
		Real* delta=factorLoading->getDeltas();
		Real sum=0.0;
		for(int j=p;j<q;j++) sum+=delta[j]*H[j+1];
			
		return sum;
	}
	
	

/**********************************************************************************
 *
 *            LMM NODES 
 *
 *********************************************************************************/
	
	
	

	void LmmNode2F::
	checkState(int t, vector<Real>& V, vector<Real>& Z, Matrix<Real>& RQt)
    {
        bool we_have_a_problem=false;
		
		if(fabs(V[n-1]-V1)+fabs(V[n-2]-V2)>0.000001){
			
			we_have_a_problem=true;
			cout << "\n\n\nInversion failure"
			     << "\nMatrix RQt: " << RQt
			     << "\n V[n-1] = " << V[n-1] << ", V1 = " << V1
			     << "\n V[n-2] = " << V[n-2] << ", V2 = " << V2;
		}
		
		for(int j=t;j<n;j++)
        if(V[j]>10.0){ 
			
			we_have_a_problem=true;
			cout << "\n\n\nLarge V_j, t = " << t 
	               << "\nj = " << j << ", V_j = " << V[j]
	               << "\n V1 = " << V1 << ", V2 = " << V2
			       << "\n\nMatrix RQt: " << RQt
	               << "Vector Zt: " << Z;
			break;
		}
	            
	    if(we_have_a_problem) exit(0);
	} // checkState()
	

	
	
	

