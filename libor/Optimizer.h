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



#ifndef martingale_optimizer_h    
#define martingale_optimizer_h

#include "TypedefsMacros.h"
#include "Array.h"                            // direct members

MTGL_BEGIN_NAMESPACE(Martingale)


// we are using
class LowDiscrepancySequence;



/*! \file Optimizer.h
 * <p>BFGS, Downhill Simplex optimizers and global search along a Sobol sequence.
 * The BFGS optimizer is modified so as to
 * cope with objective functions \f$f(u)\f$ with step function like 
 * qualities. Our main applications are expectations
 * \f[f(u)=E[X(u)],\f]
 * where the random variable \f$X(u)\f$ depends on the parameter vector u.
 * Even if this function is smooth the smoothness is lost as soon as the expectation
 * is computed as a finite sample mean.
 *
 * <p>The optimizers follow the ideas in Numerical Recipes with some modifications. 
 * If this works for you be thankful. 
 * Multidimensional optimization is a nontrivial undertaking and a high quality implementation 
 * requires a significant amount of effort and expertise. I have no expertise in this field.
 *
 * <p>We have only a single application in mind: parameter optimization
 * for exercise triggers used to price Bermudan style options. Miraculously
 * the code which we have here can do that.
 */
   
/*******************************************************************************
 *
 *                     OPTIMIZER INTERFACE
 *
 ******************************************************************************/
 


    
/** Interface to all optimizers.
 */
class Optimizer {
	
protected:
    
    int n;                 // number of variables


public:
		
	/** The objective function. */
	virtual Real f(const RealArray1D& x) = 0;
    
    
    /** <p>Search for vector x minimizing the function f(x).
     *
     * @returns the minimizing vector x.
     */
    virtual const RealArray1D& search() = 0;
    
    
    /** @param n_ dimension of argument vector.
     */
    Optimizer(int n_) : n(n_) { }
	virtual ~Optimizer(){ }
	
    
}; // end Optimizer



   
/*******************************************************************************
 *
 *                     		DOWNHILL SIMPLEX
 *
 ******************************************************************************/


/** Downhill simplex optimizer as described in NR.
 *
 * @author  Michael J. Meyer
 */
class DownhillSimplex : public Optimizer {
    
  
   Array1D<RealArray1D*> vertices;         // vertices[i], i=0,...,n, are the 
                                           // vertices of the simplex
   RealArray1D sum;              // sum of all vertices
   RealArray1D newVertex;        // next vertex
   RealArray1D y;                // function values at the vertices
   
   Real fx;      // function at current point
   
   int min,      // index of minimum vertex
       max,      // index of vertex with highest function value
       max2,     // index of second highest vertex
       nStep;    // number of steps
   
   // The barycenter is updated dynamically whenever a vertex changes
   // since this is cheaper than doing a full computation in each step.
   
   bool verbose;
    
public:

   RealArray1D& vertex(int i){ return *(vertices[i]); }
   Real& vertex(int i, int j){ return (*(vertices[i]))[j]; }

// CONSTRUCTOR
    
    /** @param x starting point (initial simplex built around it).
     *  @param delta size of initial simplex.
     *  @param steps number of steps.
     *  @param vbose messages during computation.
     */
    DownhillSimplex(const RealArray1D& x, Real delta, int steps, bool vbose);
   
    ~DownhillSimplex();
	
	
    /** Initializes the vertices of the intial simplex with the function values.
	 *  We do this at thec start of the search since we can't call a virtual 
	 *  member frunction in the constructor.
	 */
	void setInitialConditions();

    /** Searches for and returns the minimizing vector. 
	 *  Implemented in the header because implementation in the .cc file
	 *  causes an inexplicable syntax error message.
	 */
    const RealArray1D& search();


	 
private:
 
	 
// REFLECTING, CONTRACTING, KEEPING TRACK OF HIGH-LOW POINTS	 
    
    /** Set the indices for best (min), worst (max) and second worst (max2)
     * vertex
     */
    void setHiLowIndices();
	
	
    /** <p>vertex[i] is reflected at the barycenter of the convex hull of the remaining 
	 *  vertices (the opposing face) by a factor k. If k is
     *  positive we land at the other side of the face, if negative we land
     *  at the same side of the face. Writes the resulting vertex into the field
     *  <code>newVertex</code>.</p>
     *
     * <p>If the new vertex is better than the worst existing (highest) vertex
     * this vertex is exchanged with the new vertex. Assumes that the variable 
     * center is kept current at the barycenter of the whole simplex.</p>
     *
     * @param k scaling of the distance of the moving vertex from the
     * center of the opposing face.
     * @returns function value at reflected vertex.
     * 
     */
    Real reflectVertex(int i, Real k);
	
	
    /** <p>Contracts the simplex by a factor of two in the direction of
     *  vertex i (which remains unaffected).</p>
     *
     * @param i number of vertex which remains fixed.
     * 
     */
    void contract(int i);
	
    
    /** <p>Exchanges vertices[max] (the worst vertex) with newVertex 
     *  and updates the barycenter and function values.
     *  The value w=f(newVertex) is passed as a parameter to
     *  save on function evaluations.</p>
     *
     * @param w must be f(newVertex).
     */
    void replaceWorst(Real w);
	
            
}; // end DownhillSimplex


	
/** Downhill Simplex with particular objective function given by
 *  a function pointer.
 */
class ConcreteDownhillSimplex : public DownhillSimplex {
	
	 Real (*of)(const RealArray1D& x);   // pointer to the objective function
	
public:
	
	 /** @param f pointer to objective function. 
	  */
	 ConcreteDownhillSimplex
	 (Real (*f)(const RealArray1D&), RealArray1D& x, Real delta, int steps, bool verbose) :
     DownhillSimplex(x,delta,steps,verbose),
	 of(f)
     {  
		 // now we can call f
		 setInitialConditions();
     }
	 
	 Real f(const RealArray1D& x){ return (*of)(x); }

}; 



/*******************************************************************************
 *
 *                        BFGS
 *
 ******************************************************************************/


/** <p>BFGS minimizer for multivariate functions as described in NR
 *  and with some modifications to cope with the lack of smoothness 
 *  of the objective function.
 *  Code clarity is preferred to efficiency. The code in NR is more compact
 *  and consequently more efficient. Here the view is taken that the bulk of
 *  the computational effort is spent on evaluations of the function to be
 *  minimized.</p>
 *
 * <p><b>Modification:</b> the BFGS algorithm uses a Newton step of a length
 * designed to hit the minimum quickly if the function is a parabola. This 
 * implementation always takes a time step of length <code>maxstep</code> in 
 * the current direction and then performs a line search backward toward the 
 * last point. Why do we do that? The functions we are trying to minimize are
 * of the form 
 * \f[f(u)=E[X(u)],\f]
 * where X(u) is a random variable depending on the parameter
 * vector u and E denotes the expectation as usual.</p>
 *
 * <p>The function f(u) may very well be smooth but unfortunately we cannot work 
 * with the true expectation but instead have to work with the Monte Carlo computed
 * sample mean
 * \f[f_N(u)=(X_1(u)+...+X_N(u))/N,\f]
 * where the \f$X_j(u)\f$ are independent observations of \f$X(u)\f$. 
 * This function will most likely 
 * not be smooth as a function of the parameter vector u. Instead the parameter u 
 * will have to move at least for a minimal threshhold value before the 
 * finite sample of \f$X(u)\f$ recognizes a difference. We should
 * therefore think of \f$f_N(u)\f$ as blocky, instead of a smooth  
 * parabola we will have a sequence of concentric flat terraces descending 
 * toward the minimum.</p>
 *
 * <p>Obviously we have to choose a step size which gets us off each terrace. 
 * Likewise we have to compute the function gradients with variable increments
 * of sufficient size lest we find ourselves with a zero gradient despite being 
 * far from the local minimum.</p>
 *
 * @author  Michael J. Meyer
 */
class BFGS : public Optimizer {
	
protected:
	
    // some constants 
    static Real
    /** Termination in line search if t falls below this threshold.*/
    EPSX,
    /** Proportional decrease in function value enforced at each step in line search.*/
    KF,
    /** Termination if norm of gradient falls below this threshold.*/
    EPSG,
    /** Small number to correct division by zero.*/
    EPS;

   
    Real stepmax;                 // maximum steplength in the line search
    int nVals,                    // maximum number of function evaluations
        fVals,                    // current number of function evaluations
	    restarts,                 // number of times the optimizer has been restarted
	    nRestarts;                // number of restarts budgeted

    // keep current
    RealArray1D 
	
	         x0,                  // the last point
             grad0,               // gradient at x0
             x,                   // the current point
             xDelta,              // x-x0
             grad,                // gradient at the current point x
             z,                   // next point to be computed (workspace)
             d,                   // current line direction
             gDelta,              // gradient difference grad-grad0
             hdg,                 // H(grad-grad0)
             u,                   // the additional vector in the bfgs update
             h;                   // directional increments to compute gradient
			 
    Real     fx0,                 // f(x0), function at the last point
             fx;                  // f(x), function at the current point

                 
    RealArray2D H;              // approximate inverse Hessian at the
                                  // current point.
    
    bool verbose;

    
public:

// ACCESSORS

    /** Dimension of the search space. */
    int getDimension() const { return n; }
    
    /** Current gradient.*/
    const RealArray1D& getGrad() const { return grad; }
    
    /** Current point.*/
    const RealArray1D& getX() const { return x; }
    
    /** Current direction.*/
    const RealArray1D& getD() const { return d; }
	

// CONSTRUCTOR
	


    /** 
     * @param u starting point of min search
     * @param vals maximum number of function evaluations allowed
     * @param maxstep maximum step length in line search
     * @param k variable increments dx_j in the gradient computation.
	 * @param resets number of times the optimizer is restarted.
     * @param vbose messages about results, default is <code>false</code.
     */
    BFGS(const RealArray1D& u, int vals, Real maxstep, 
	     const RealArray1D& k, int resets=2, bool vbose=false);
	
		    	
    /** Sets function value, gradient and initial direction by making calls to {@link Optimizer#f}. 
     *  We do this at the start of the search since we can't call a virtual member frunction in the 
	 *  constructor.
     */
    void setInitialConditions();
    
		
           
// MIN SEARCH                
      
     /** Unconstrained search for the minimum of the function 
      *  {@link Optimizer#f}. Vector is owned by this. Do not delete.
      *
      * @returns the minimizing vector <code>x</code>.
      */
     const RealArray1D& search();

	
	
    
private:
		

	
    /** Relative size of the vector x relative to the vector y
	 *  measured coordinate by coordinate rather than 
     *  through norms (global). Used for termination criteria.
     */
    Real relativeSize(RealArray1D& x, RealArray1D& y);
	
    
    /** Checks for termination on the size of 
     *  the gradient (approximately zero?)
     */     
    bool gradientIsZero();
    
    
    /** The bfgs update of the approximate inverse Hessian.
     */
    void bfgsUpdate();
	
	
	/** Resets the Hessian to the identity. This is a restart of the optimizer.
	 */
	void resetHessian();
	
	
// GRADIENT COMPUTATION  
    
    /** <p>Computes the gradient of {@link Optimizer#f} at the point 
     *  <code>x</code> by forward differencing against n other points 
     *  (x+h*e_j).</p>
     *
     *  <p>The gradient is crucial for this algorithm. If an exact form
     *  for the gradient is available by all means use that and override this 
     *  method.</p>
     *
     * <p>The default implementation writes the result into the field 
     * <code>gradient</code> and returns a reference to this field. 
     * No new memory is allocated. It is not necessary to follow this pattern 
     * if the method is overwritten.</p>
     *
     * @param x point at which the gradient is computed.
     * @param fx function value <code>f(x)</code>.
     */
    const RealArray1D& gradF(RealArray1D& x, Real fx);
	
   
    /** <p>Computes the gradient of {@link Optimizer#f} at the point 
     *  <code>x</code> by central finite differencing from 2n points 
     *  (x+-h*e_j).</p>
     *
     *  <p>The gradient is crucial for this algorithm. If an exact form
     *  for the gradient is available by all means use that and override this 
     *  method.</p>
     *
     * <p>The default implementation writes the result into the field 
     * <code>gradient</code> and returns a reference to this field. 
     * No new memory is allocated. It is not necessary to follow this pattern 
     * if the method is overwritten.</p>
     *
     * @param x point at which the gradient is computed.
     */
    const RealArray1D& gradcdF(RealArray1D& x);
	
	
// BACK TRACKING DURING LINE SEARCH                 
    
    /** Backtracking in the line search for {@link Optimizer#f}..
     *  Backtracking from a Newton step t which is too long after the first 
     *  simple backtrack. At this point we have the current line parameter t1
     *  and the line parameter t2 before that, the function values
     *  f0=f(x) at the start of the line and f1=f(x+t1*d), f2=f(x+t2*d) 
     *  as well as the slope k of the descent of f along the line
     *  from x in the direction of d.
     *
     * @returns the new and smaller value t of the line parameter.
     */
    Real backTrack(Real t1, Real t2, Real f1, Real f2, Real k);


//  LINE SEARCH
  
    
    /** <p>Line search for function minimum from current point in
     *  direction d. This algorithm assumes that d
     *  is a direction in which f decreases and then searches for
     *  a new point in this direction at which an acceptable decrease in the
     *  value of f occurs.</p>
     *
     * <p>A step of length maxstep tried first and a backtracking strategy 
     * employed if this step does not lead to a decrease in the value of f.
     * The resulting point is written into the field u.
     * The method does not update the fields, do this in the 
     * main loop immediately after the call to lineSearch.</p>
     */
    void lineSearch();

	
            
}; // end BFGS


	
/** BFGS with particular objective function given by
 *  a function pointer.
 */
class ConcreteBFGS : public BFGS {
	
	 Real (*of)(const RealArray1D& x);   // pointer to the objective function
	
public:
	
	 /** @param f pointer to objective function. 
	  */
	 ConcreteBFGS
	 (Real (*f)(const RealArray1D&), RealArray1D& x, int nVals, 
      Real stepmax, RealArray1D& h, int nRestarts=3, bool verbose=false) :
     BFGS(x,nVals,stepmax,h,nRestarts,verbose),
	 of(f)
     {   
		 // now we can call f
		 setInitialConditions();
	 }
	 
	 Real f(const RealArray1D& x){ return (*of)(x); }

}; 




/*******************************************************************************
 *
 *              SOBOL SEARCH 
 *
 ******************************************************************************/

/** Constrained global search using a Sobol sequence. Slow and dumb but thorough.
 *  Searches a rectangular window centered at the current best point with a Sobol
 *  sequence then restarts itself at the next best point with a slightly contracted 
 *  window and fewer search points.
 *
 * @author  Michael J. Meyer
 */
class SobolSearch: public Optimizer {
	
protected:
	
	int nPoints;                       // total number of search points
	Real q;                            // window contraction factor
	RealArray1D xOpt;                  // current best point
	RealArray1D x;                     // current point
	RealArray1D d;                     // window: all u with x_j-d_j < u_j < x_j+d_j
	
	bool verbose;
	
	LowDiscrepancySequence* lds;
	
	
public:
	
	/** @param x0 initial point.
	 *  @param nVals total number of function evaluations.
	 *  @param delta initial window all u with \f$x0_j-\delta_j<u_j<x0_j+\delta_j\f$.
	 *  @param vbose announce each new min during search.
	 */
	SobolSearch(const RealArray1D& x0, int nVals, const RealArray1D& delta, bool vbose=false);
	
	/** Wether or not the vector u is in the search domain.
	 *  This is the default implementation (true, unconstrained search).
	 */
	virtual bool isInDomain(const RealArray1D& u) const { return true; }
	
	const RealArray1D& search();
	
}; // end SobolSearch
	
	
	
	
/** SobolSearch with particular objective function given by
 *  a function pointer.
 */
class ConcreteSobolSearch : public SobolSearch {
	
	 Real (*of)(const RealArray1D&);   // pointer to the objective function
	
public:
	
	 /** @param f pointer to objective function. 
	  */
	 ConcreteSobolSearch
	 (Real (*f)(const RealArray1D&), const RealArray1D& x, 
	  int nVals, const RealArray1D& delta, bool verbose=false) :
     SobolSearch(x,nVals,delta,verbose),
	 of(f)
     {    }
	 
	 Real f(const RealArray1D& x){ return (*of)(x); }

};     
		
	
	
	
// TEST FUNCTIONS

/** Some examples of objective functions. */
namespace ObjectiveFunction {
	
	Real function_1(const RealArray1D&);
	
} // end namespace



MTGL_END_NAMESPACE(Martingale)

#endif
