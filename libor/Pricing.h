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
#include "Edge.h"
#include "PathGenerator.h"
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
 * There is a price however: a user who wishes to implement her own Option and
 * Lattice classes and make use of these function templates will have to read the 
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
	

/** The forward price of theOption when priced in thelattice.
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
		node->setPi(theOption->forwardPayoff(node));
		theNode++;
	}
	
	// backward computation through earlier nodes
	for(int t=expiry-1; t>=0;t--){
		
	    // can we exercise at time t?
		bool exercise_t = theOption->isExercisable(t);
		
		vector<NodeType*>* nodes_t=theLattice->getNodeList(t);
	    theNode=nodes_t->begin();
    	while(theNode!=nodes_t->end()) {
	    
			NodeType* node = *theNode; 
			vector<Edge*>* edges=node->getEdges();
				
			 Real E_t=0.0;              // E_t(pi_{t+1})
			 vector<Edge*>::const_iterator theEdge=edges->begin();
			 while(theEdge!=edges->end()) {
	        
		         Edge* edge = *theEdge;
				 Real p=edge->probability;          // transition prob. along the edge
				 Real pi_t1=edge->node->getPi();    // pi_{t+1} at target node of edge
				 E_t+=p*pi_t1;
				 ++theEdge;
	         }	
			 
			 if(!exercise_t) node->setPi(E_t);	
			 else            node->setPi( max( E_t, theOption->forwardPayoff(node) ) );	
			 ++theNode;
	    } // end while(theNode)
			
	} // end for t
	
	// report the price at the root node
    NodeType* root=theLattice->getRoot();
	return root->getPi();
	
} // end 



/** The forward price of theOption computed from nPaths paths of thePathGenerator.
 */
template<typename OptionType>
Real monteCarloForwardPrice(OptionType* theOption, int nPaths)
{
	PathGenerator* thePathGenerator=theOption->getPathGenerator();
	Real sum_payoff=0.0;
	for(int i=0;i<nPaths;++i){
		
		thePathGenerator->newPath();
		sum_payoff+=theOption->forwardPayoffAlongCurrentPath();
	}
	return sum_payoff/nPaths;
	
} 


/** The beta coefficient of option payoff and control variate 
 *  computed along N paths of thePathGenerator (book, 2.8).
 *  See {@link #controlledMonteCarloForwardPrice} for the assumptions.
 */
template<typename OptionType>
Real betaCoefficient(OptionType* theOption, int N)
{
	PathGenerator* thePathGenerator=theOption->getPathGenerator();
	// Recall that Cov(X,Y)=E(XY)-E(X)E(Y) and Var(X)=E(X^2)-E(X)^2.
	// X = Payoff, Y = ControlVariate
	Real X,Y,
	     sum_X=0.0, sum_Y=0.0,
	     sum_XX=0.0, sum_XY=0.0;
        
    for(int i=0;i<N;i++){
			
       thePathGenerator->newPath();
	   X = theOption->forwardPayoffAlongCurrentPath();
	   Y = theOption->controlVariateAlongCurrentPath();
            
       sum_X+=X; sum_Y+=Y;
       sum_XX+=X*X; sum_XY+=X*Y;
    }
        
    return (N*sum_XY-sum_X*sum_Y)/(N*sum_XX-sum_X*sum_X);
        
}


/** The forward price of theOption computed from nPaths paths of thePathGenerator
 *  supplied by the option and using the control variate implemented for the option.
 */
template<typename OptionType>
Real controlledMonteCarloForwardPrice(OptionType* theOption, int nPaths)
{
	PathGenerator* thePathGenerator=theOption->getPathGenerator();
	Real controlVariateMean=theOption->controlVariateMean();
	Real beta=betaCoefficient(theOption,nPaths/10);
	Real sum_payoff=0.0, 
	     payoff, 
	     controlVariate;
	
	int N=9*nPaths/10;
	for(int i=0;i<N;++i){
		
		thePathGenerator->newPath();
		payoff=theOption->forwardPayoffAlongCurrentPath();
		controlVariate=theOption->controlVariateAlongCurrentPath();
		sum_payoff+=payoff-beta*(controlVariate-controlVariateMean);
	}
	return sum_payoff/N;	
} 



/** The correlation of the forward payoff of theOption with its control variate
 *  computed from N paths of thePathGenerator supplied by the option.
 */
template<class OptionType>
Real correlationWithControlVariate(OptionType* theOption, int N)
{
    PathGenerator* thePathGenerator=theOption->getPathGenerator();
	Real x,         // forward payoff
	     y,         // corresponding control variate
	     sum_x=0.0, sum_xx=0.0, 
	     sum_y=0.0, sum_yy=0.0, sum_xy=0.0;
	for(int n=0;n<N;n++)
    {
        thePathGenerator->newPath();
		x = theOption->forwardPayoffAlongCurrentPath();
		y = theOption->controlVariateAlongCurrentPath();            
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
 
