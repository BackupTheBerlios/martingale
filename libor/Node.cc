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

#include "Node.h"
#include "Array.h"
#include "LiborFactorLoading.h"
#include <iostream>
#include <cstdlib>                        // exit()
using namespace Martingale;



/**********************************************************************************
 *
 *            GENERAL NODES 
 *
 *********************************************************************************/
 
Node:: 
~Node()
{
	std::list<Edge>::const_iterator theEdge;         // pointer to Edge
	for(theEdge=edges.begin(); theEdge!=edges.end(); ++theEdge)  delete &(*theEdge); 
}


void 
Node::
printTransitionProbabilities() const
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
LmmNode(LiborFactorLoading* fl, int s, int steps) : Node(s),
factorLoading(fl), 
n(factorLoading->getDimension()),
nSteps(steps),
H(n-get_t()+1,get_t()) 
{  }


Real 
LmmNode::
X(int j) const
{
	int t=get_t();             // node in (T_{t-1},T_t]
	if((t==n)&&(j==n-1)) return (H[j]-1.0);
	if((t<n)&&(j>=t)&&(j<n)) return (H[j]-H[j+1])/H[j+1];
	// else
	cout << "\n\nLmmNode:X(j): j = " << j << " not in ["<<t<<","<<n-1<<"]"
	     << "\nTime step s = " << s
	     << "Terminating.";
	exit(0);		
}
	 

Real 
LmmNode::
Hpq(int p, int q) const
{
	int t=get_t();            // node in (T_{t-1},T_t]
	if((p<t)||(q<=p)){
		
		cout <<	"\n\nLmmNode::Hpq(p,q): t = " << t << ", p = " << p << ", q = " << q
		     << "\n Incompatible indices, terminating.";
		exit(0);
	}
		
	const RealArray1D& delta=factorLoading->getDeltas();
	Real sum=0.0;
	for(int j=p;j<q;j++) sum+=delta[j]*H[j+1];
			
	return sum;
}
	
	
void 
LmmNode::
printState() const
{
     std::cout << "\n\n Diagnostic: node at time step s = " << s
	           << "\nVector H: " << H
	           << "\nPrice pi = " << pi;
}


std::ostream& 
LmmNode::
printType(std::ostream& os) { return os << "general LMM node."; }


std::ostream& 
LmmNode2F::
printType(std::ostream& os) { return os << "2 factor LMM node."; }


std::ostream& 
LmmNode3F::
printType(std::ostream& os) { return os << "3 factor LMM node."; }




/**********************************************************************************
 *
 *            ASSET BASKET NODES 
 *
 *********************************************************************************/

	
void 
BasketNode::
printState() const
{
     std::cout << "\n\n Diagnostic: node at time step s = " << s
	           << "\nVector S: " << S
	           << "\nPrice pi = " << pi;
}


std::ostream& 
BasketNode::
printType(std::ostream& os) { return os << "general asset basket node."; }


std::ostream& 
BasketNode2F::
printType(std::ostream& os) { return os << "2 factor asset basket node."; }


std::ostream& 
BasketNode3F::
printType(std::ostream& os) { return os << "3 factor asset basket node."; }
	

	
	
	

