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


#ifndef martingale_node_h    
#define martingale_node_h

// implementation is in header
#include "TypedefsMacros.h"
#include "Array.h" 

MTGL_BEGIN_NAMESPACE(Martingale)



/*! \file Node.h
 *  <p>Nodes in a stochastic lattice for a standard Brownian motion Z.
 *  Emphasis is on smallest possible memory footprint. 
 */
 


/**********************************************************************************
 *
 *            NODES A STOCHASTIC LATTICE
 *
 *********************************************************************************/



/** <p>Nodes in a stochastic lattice for a standard Brownian motion Z.
 *  Emphasis is on smallest possible memory footprint. Nodes only store
 *  the array of edges and the state k of the Brownian motion 
 *  \f[Z_j=k[j]*a,\quad a=\sqrt{dt}\f]
 *  as a integer multiples of the tick size a. Nodes don't need to know 
 *  anything else, not at which time they live or even which lattice they 
 *  live in. The lattice knows that and all else such as the transition 
 *  probabilities. This avoids storing the same data in each of the million 
 *  nodes in the lattice.
 */
class StandardBrownianNode {
	
protected:
	
	/** the state Z_j=k[j]*a, a=sqrt(dt). */
	int* k_;

	/** Cache for option value at this node.*/
	Real  value;               
	
	/** List of edges originating at this node. */
	Array1D<StandardBrownianNode*> edges; 
	
public:
	
/** The state Z_j=k[j]*a, a=sqrt(dt). */
int* getIntegerTicks(){ return k_; }
	
/** Get the (option) value at this node. */
Real getValue() const { return value; }
/** SWet the (option) value at this node. */
void setValue(Real x){ value=x; }
	
/** List of edges */
Array1D<StandardBrownianNode*>& getEdges() { return edges; }
	
	
// CONSTRUCTOR
	
/** @param nEdge number of edges orginating from this node.
 */
StandardBrownianNode(int nEdge, const IntArray1D& k) : 
k_(new int[k.getDimension()]),
value(-1.0),
edges(nEdge)
{
	for(int i=0;i<k.getDimension();i++) k_[i]=k[i];
}
	
virtual ~StandardBrownianNode(){ delete[] k_; };


/** Diagnostic. Prints the time t, vector H and field pi.*/
std::ostream& printSelf(std::ostream& os)
{
	return
	os << "\n\nStandardBrownianNode, value V =  " << value;
}
	 
}; // end StandardBrownianNode




	
MTGL_END_NAMESPACE(Martingale)

#endif
 
