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

#ifndef martingale_pricing_h    
#define martingale_pricing_h


// function templates, fully definded in the header
// thus the need for the includes
#include "TypedefsMacros.h"
#include "Array.h"
#include <vector>               


MTGL_BEGIN_NAMESPACE(Martingale)

/*! \file Pricing.h
 * Standalone function templates for Monte Carlo and lattice pricing of options.
 * The purpose is to take the pricing code out of the class Option and also to
 * to get around rigid interfaces which have to be implemented.
 *
 * <p>The functions are parametrized both by the type of Lattice and the type of 
 * Option. In the case of lattice we have no choice since Lattice is itself a
 * class template which defines a family of types (parametrized by NodeType)
 * which have no common base type. The type of Option is made a template parameter
 * to gain more flexibility for the various types of Options: the concrete Option 
 * classes do not have to implement all the member functions which are used here.
 *
 * <p>For example the base class Option itself implements nothing and the subclasses
 * of Option which become the actual template parameters can implement these methods 
 * selectively to enable selected functionality: lattice pricing or Monte Carlo
 * pricing or both. The function templates are instantiated only at the point
 * of actual use (function call). The compiler then looks for the necessary pieces 
 * (here member functions) from the actual template parameters and if all are found 
 * the function comes into being otherwise a compiler error results.
 *
 * <p>This approach provides enormous flexibility. Classes no longer have to 
 * conform to well defined interfaces as in a traditional inheritance hierachy.
 * There is a price however: if you want to implement your own Option and
 * Lattice classes and make use of these function templates you will have to read the 
 * source code to understand what methods these classes have to provide and how 
 * they are used.
 */

/**********************************************************************************
 *
 *                       PRICING
 *
 *********************************************************************************/
 


/** <p>Standalone function templates for Monte Carlo and lattice pricing of options.
 */
namespace Pricing {
	
using std::vector;
	

/** <p>The forward price of theOption when priced in theLattice. Read the source code 
 *  to see what member functions the types LatticeType and OptionType must
 *  implement.
 *  <p>The routine is written under the assumption that the transition probabilities
 *  are computed by the lattice and are both time and state independent. This is rather
 *  restrictive but true for the lattices we consider. The code can easily be rewritten
 *  to cover the general case of both time and state dependent transition probabilities.
 *  Simply define <code>LatticeType::transitionProbability</code> to be a function of
 *  the edge i, time t and the node n. For reasons of speed this approach is not taken 
 *  here.
 *  <p>Please read the source code to determine the precise syntactic assumptions which 
 *  the code makes about the type of lattice <code>LatticeType</code>, the type of nodes
 *  within the lattice <code>LatticeType::NodeType</code> and the type 
 *  <code>OptionType</code> of option to be priced in the lattice
	.
 */
template<typename LatticeType, typename OptionType>
Real latticeForwardPrice(LatticeType* theLattice, OptionType* theOption)
{
	Real  expiration=theOption->getExpiration(),
	      dt = theLattice->getTimeStep();             // size of time step	
    int   timeSteps=theLattice->getTimeSteps(),       // number of time steps in lattice
	      expiry = (int)(expiration/dt);              // time steps to expiry
	
	// is the lattice built far enough:
	if(timeSteps*dt<expiration-0.00000001){
	   
		cout << "\n\nPricing::latticeForwardPrice:"
	        << "\nLattice not built all the way to expiration."
	        << "\nTerminating.";
		exit(1);
	}
	 
  	typedef typename LatticeType::NodeType NodeType;
	
	// run through the list of nodes at time T=expiration and set the known value pi_T=h_T
	vector<NodeType*>* nodes_at_expiry=theLattice->getNodeList(expiry);
	// *theNode is pointer to Node
    vector<NodeType*>::iterator theNode=nodes_at_expiry->begin();
	while(theNode!=nodes_at_expiry->end()) {
	
		NodeType* node=*theNode;
		node->setValue(theOption->forwardPayoff(node,theLattice,expiry));
		theNode++;
	}
	
	// backward computation through earlier nodes
	for(int t=expiry-1; t>=0;t--){
		
	    // can we exercise at time t (continuous time!)
		bool exercise_t = theOption->isExercisable(t*dt);
		
		vector<NodeType*>* nodes_t=theLattice->getNodeList(t);
	    theNode=nodes_t->begin();
    	while(theNode!=nodes_t->end()) {
	    
			NodeType* node = *theNode; 
			const Array1D<NodeType*> edges=node->getEdges();    // list of edges
			int nEdge = edges.getDimension();                   // number of edges
				
			Real V=0.0;           // option forward price at the node
			for(int i=0;i<nEdge;i++) {
	        
				 Real V_i=edges[i]->getValue();                   // value at node *edges[i]
				 Real p_i=theLattice->transitionProbability(i);   // transition probability along edge_i
				 V+=p_i*V_i;
	         }	
			 
			 if(!exercise_t) node->setValue(V);	
			 else  node->setValue( max( V, theOption->forwardPayoff(node,theLattice,t) ) );	
			 ++theNode;
	    } // end while(theNode)
			
	} // end for t
	
	// report the price at the root node
    NodeType* root=theLattice->getRoot();
	return root->getValue();
	
} // end 



/** The forward price of theOption computed from N sample payoffs
 *  compounded forward to the horizon.
 */
template<typename OptionType>
Real monteCarloForwardPrice(OptionType* theOption, int N)
{
	Real sum_payoff=0.0;
	for(int i=0;i<N;++i) sum_payoff+=theOption->nextForwardPayoff();
	return sum_payoff/N;	
} 


/** The beta coefficient of option payoff and control variate 
 *  computed from N forwardPayoff - controlVariate pairs (book, 2.8).
 *  See {@link #controlledMonteCarloForwardPrice} for the assumptions.
 */
template<typename OptionType>
Real betaCoefficient(OptionType* theOption, int N)
{
	// Recall that Cov(X,Y)=E(XY)-E(X)E(Y) and Var(X)=E(X^2)-E(X)^2.
	// X = Payoff, Y = ControlVariate
	Real X,Y,
	     sum_X=0.0, sum_Y=0.0,
	     sum_XX=0.0, sum_XY=0.0;
        
    for(int i=0;i<N;i++){
		
	   // next forwardPayoff - controlVariate pair
	   const RealArray1D& Z = theOption->nextControlledForwardPayoff();
	   X = Z[0];       // the payoff
	   Y = Z[1];       // the control variate
            
       sum_X+=X; sum_Y+=Y;
       sum_XX+=X*X; sum_XY+=X*Y;
    }
        
    // divide numerator and denominator by N^2
	return (N*sum_XY-sum_X*sum_Y)/(N*sum_XX-sum_X*sum_X);
        
}


/** The forward price of theOption computed from nPaths paths of thePathGenerator
 *  supplied by the option and using the control variate implemented for the option.
 *  See book, 2.8.
 */
template<typename OptionType>
Real controlledMonteCarloForwardPrice(OptionType* theOption, int nPaths)
{
	Real controlVariateMean=theOption->controlVariateMean();
	int M=1024;          // complete cycle for the Sobol generator
	Real beta=betaCoefficient(theOption,M);
	Real sum_payoff=0.0, 
	     payoff, 
	     controlVariate;
	
	for(int i=0;i<nPaths;++i){
		
	    // next forwardPayoff - controlVariate pair
	    const RealArray1D& Z = theOption->nextControlledForwardPayoff();
	    payoff = Z[0];       
	    controlVariate = Z[1];       
		sum_payoff+=payoff+beta*(controlVariateMean-controlVariate);
	}
	return sum_payoff/nPaths;	
} 



/** The correlation of the forward payoff of theOption with its control variate
 *  computed from N paths of thePathGenerator supplied by the option.
 */
template<class OptionType>
Real correlationWithControlVariate(OptionType* theOption, int N)
{
	Real x,         // forward payoff
	     y,         // corresponding control variate
	     sum_x=0.0, sum_xx=0.0, 
	     sum_y=0.0, sum_yy=0.0, sum_xy=0.0;
	for(int n=0;n<N;n++)
    {
	    // next forwardPayoff - controlVariate pair
	    const RealArray1D& Z = theOption->nextControlledForwardPayoff();
	    x = Z[0];       // the payoff
	    y = Z[1];       // the control variate    
		sum_x+=x; sum_y+=y;
        sum_xx+=x*x; sum_yy+=y*y; sum_xy+=x*y;
    }
    
	//when divided by N^4: 
    Real Var_xVar_y=((N*sum_xx-sum_x*sum_x)*(N*sum_yy-sum_y*sum_y));
    return (N*sum_xy-sum_x*sum_y)/sqrt(Var_xVar_y);
}


}; // end Pricing




MTGL_END_NAMESPACE(Martingale)

#endif
 
