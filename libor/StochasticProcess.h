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

#ifndef martingale_stochasticprocess_h
#define martingale_stochasticprocess_h

// template class fully defined in header
#include "Array.h"
#include "RandomObject.h"
#include "StochasticGenerator.h"

/*! \file StochasticProcess.h
 * Stochastic processes and associated objects such as stopping times, hitting times
 * and first exit times of euclidean regions, path functionals.
 */
 
// StoppingTime, Region, HittingTime, FirstExitTime, StochasticProcess, PathFunctional. 


MTGL_BEGIN_NAMESPACE(Martingale)



/********************************************************************************************
 *
 *                STOPPING TIMES, HITTING TIMES, FIRST EXIT TIMES
 * 
 ********************************************************************************************/
 
/** Stopping time for a stochastic Process. See book, 3.3.
 */
class StoppingTime {
	
protected:
	
	int T;  // number of time steps to horizon
	
public:
	
	/** @param T_oo number of time steps to horizon. */
	StoppingTime(int T_oo) : T(T_oo) {  }

    /** <p>Returns true if we must stop at time t,
     * false otherwise. In concrete implementations make sure
     * it returns true as soon as we hit the horizon.</p>
     */
    virtual bool stop(int t) = 0;
	
   
	/** The smallest time t such that stop(t) returns true.
     */
   int valueAlongCurrentPath()
   {
	   int t=0; while(!stop(t)&&(t<T)) t++;
	   return t;
   }
    
}; //end StoppingTime



/** A region (set) in a universe populated with elements of type RangeType.
 */
template<typename RangeType>
class Region {
	
public:

  /** The object x is either in the region or not.
   */
  virtual bool isMember(const RangeType& x) = 0;
  
  
}; //end Region


template<typename RangeType, typename ScalarType>
class StochasticProcess;     // defined below


/** Stopping time triggering as soon as the stochastic process X
 * hits a region D in its range. This range is populated with elements of
 * the vector type RangeType with components of the type ScalarType.
 * See {@link StochasticProcess}.
 *
 * @author  Michael J. Meyer
 */
template<typename RangeType, typename ScalarType>
class HittingTime : public StoppingTime {

   StochasticProcess<RangeType,ScalarType>* X;
   Region<RangeType>* D;
	
public:
   
   /** @param Y Process hitting region G.
    *  @param G Region hit by process Y.
    */
   HittingTime
   (StochasticProcess<RangeType,ScalarType>* Y, Region<RangeType>* G) : 
   StoppingTime(Y->getT()), X(Y), D(G) 
   {  }

   
   /** stop as soon as X(t) hits D or t=horizon.*/
   bool stop(int t){ return ((D->isMember(X->currentPath(t)))||(t==T)); }
          
}; //end HittingTime




/** Stopping time triggering as soon as the stochastic process X
 * exits a region D in its range. This range is populated with elements of
 * the vector type RangeType with components of the type ScalarType.
 * See {@link StochasticProcess}.
 *
 * @author  Michael J. Meyer
 */
template<typename RangeType, typename ScalarType>
class FirstExitTime : public StoppingTime {

   StochasticProcess<RangeType,ScalarType>* X;
   Region<RangeType>* D;
	
public:
   
   /** @param Y Process exiting region G.
    *  @param G Region exited by process Y.
    */
   FirstExitTime
   (StochasticProcess<RangeType,ScalarType>* Y, Region<RangeType>* G) : 
   StoppingTime(Y->getT()), X(Y), D(G) 
   {  }

   
   /** Stop as soon as X(t) leaves D or t=horizon.*/
   bool stop(int t){ return (!(D->isMember(X->currentPath(t)))||(t==T)); }
       
       
}; //end FirstExitTime



/********************************************************************************************
 *
 *                       REGIONS IN EUCLIDEAN SPACE
 * 
 ********************************************************************************************/ 

/** Region G in Euclidean space, RangeType=RealVector. 
 */
class EuclideanRegion : public Region< RealVector > {
	
protected:
	
	int dim;  // dimension
	
public:
	
	int getDimension() const { return dim; }
	
	/** @param d the dimension. */
	EuclideanRegion(int d) : dim(d) { }
	
	/** Projection of x onto the boundary of G. */
	virtual RealVector& boundaryProjection(const RealVector& x) = 0;
	
	/** Intersection of the straight line from u to v with the boundary  
	 *  of G. It is assumed that u is inside and v outside of G. 
	 *  Simply moves toward the boundary along this straight line using
	 *  continued bisection N times and then projects onto the boundary.
	 */
	RealVector& boundaryIntersection
	(const RealVector& u, const RealVector& v, int N=5) 
    {
		RealVector x(u), y(v), z(u);
		for(int i=0;i<N;i++){          
			
			z=x; z+=y; z*=0.5;
			if(isMember(z)) x=z; else y=z;
		}
		return boundaryProjection(y);
	} // end boundaryIntersection
			
}; // end EuclideanRegion


// BALL CENTERED AT THE ORIGIN

/** Ball in Euclidean space centered at the origin.
 */
class Ball : public EuclideanRegion {
	
	Real R;         // radius
	
public:
	
	Real getRadius() const { return R; }
	
	/** @param d dimension
	 *  @param r radius
	 */
	Ball(int d, Real r=1.0) : EuclideanRegion(d), R(r) {  }
	
	bool isMember(const RealVector& x)
    {
		Real sum=0;
		for(int i=0;i<dim;i++) sum+=x[i]*x[i];
		return(sum<R*R);
	}
	
	/** Projects x on the boundary of he ball. */
	RealVector& boundaryProjection(const RealVector& x)
    {
		RealVector& u=*(new RealVector(x));
		u*=R/x.norm();
		return u;
	}
	
}; // end Ball
	 
	
	
	


/********************************************************************************************
 *
 *                            STOCHASTIC PROCESS
 * 
 ********************************************************************************************/
 

/** <p> Stochastic process with arbitrary values of vector type <code>RangeType</code>
 *  with components of type <code>ScalarType</code> (similar to class {@link RandomObject},
 *  Dimension 1 with <code>RangeType=ScalarType</code> is allowed and both types could be
 *  <code>std::complex</code> for example). The only reference to these types occurs in the 
 *  method {@link #operator()(int,StoppingTime*)} which returns a pointer to 
 *  <code>RandomObject</code>.
 *  Unfortunately this method is fundamental in the theory of stochastic processes 
 *  (see book, 3.3) and so we cannot avoid the template parameters of 
 *  the class <code>RandomObject</code>.</p>
 *
 * <p><b>Path computation:</b> the basic procedure is to continue a path which 
 * has already been realized until time t from this time t to the horizon T
 * (branching at time t).</p>
 *
 * <p>Here we are simulating future scenarios conditional on the state at time t.
 * This is what is needed for the computation of conditional expectations 
 * conditioning on all information available at time t. The computation of 
 * entire paths is then merely a continuation from time t=0.</p>
 *
 * <p> Assuming a setup where paths are computed as a sequence of time steps
 * path computation is reduced to the single abstract method {@link #timeStep}.
 * If this is not the case the method {@link #timeStep} can be defined to be empty
 * and the path computation methods overridden.
 *
 * <p> Time is measured in discrete units (the time step dt). Integer time t
 * corresponds to continuous time t*dt. The case dt=1 is the special case of a 
 * sequential stochastic process.</p>
 *
 * @param RangeType range of process (possibly vectors).
 * @param ScalarType type of RangeType vector components.
 *
 * @author  Michael J. Meyer
 */
template<typename RangeType=Real, typename ScalarType=Real>
class StochasticProcess {
	
protected:
	
	int dim;      // dimension of range
    int T;        // number of time steps to horizon
	
public:
	
	/** Dimension of range.*/
    int getDimension() const { return dim; }
    
    /** Number of time steps to horizon.*/
    int getT() const { return T; }
	
	/** <p>The current path at time t.</p>
     *
     * @param t discrete time.
     */
    virtual RangeType& currentPath(int t) = 0;
 
        
// CONSTRUCTOR
	
	/** Constructor, scalar process, range dimension = 1.
     * @param T_oo number of time steps to horizon.
     */
    StochasticProcess(int T_oo) : dim(1), T(T_oo) {  }
    
    /** Constructor.
	 * @param d dimension of range.
     * @param T_oo number of time steps to horizon.
     */
    StochasticProcess(int d, int T_oo) : dim(d), T(T_oo) {  }
	
	virtual ~StochasticProcess(){ }

                                             
// TIME STEPS
	
    /** <p>Evolves a path from discrete time [0,t] to time t+1, that is, 
     *  computes \f$X_{t+1}\f$ from \f$X_u, u\leq t\f$.</p>
     *
     * @param t discrete time.
     */
    virtual void timeStep(int t) = 0;
	
   

    /** <p>Computes \f$X_s\f$ from \f$X_u, u\leq t\f$.</p>
     *
     *  <p> Sometimes this can be accomplished in a single step without stepping
     *  through intermediate times t+1,...,s-1. In this case path computations
     *  can often be sped up by sampling the path only at the times s which are 
     *  needed.</p>
     *
     *  <p>Default implementation: individual time steps, no efficiency gain.
     *  This is meant to be overridden where possible but a default is
     *  provided for convenience.</p>
     *
     * @param t current time.
     * @param s future time.
     */
    virtual void timeStep(int t, int s){  for(int u=t;u<s;u++)timeStep(u);  }
    

// PATHS

    /** <p>Computes a new path segment from time t to time s through all 
     *  intermediate times t+1,...,s-1. It is assumed that a path has been
     *  computed up to time t.
     *  Default implementation: calls to <code>timeStep</code>.</p>
     *
     * @param t current time.
     * @param s future time.
     */
    virtual void newPathSegment(int t, int s){  for(int u=t;u<s;u++)timeStep(u);  }
    
      
    
    /** <p>Computes a new path segment from time t to the random time tau&gt=t
     *  and returns the value of the time tau >= t when the path is stopped.
     *  Default implementation: calls to timeStep.</p>
     *
     * @param t current time
     * @param tau random future time
     */
    virtual int newPathSegment(int t, StoppingTime* tau)
    {
        int s=t;
        while(!(tau->stop(s))){ timeStep(s); s++; }
        return s;
     }
     
     
	/** <p>Computes a new path segment from time 0 to time t.
     * Default implementation: calls to timeStep.</p>
     *
     * @param t time when path is stopped.
     */
    virtual void newPathSegment(int t){  newPathSegment(0,t);  }
	
    
    /** <p>Computes a new path segment from time t=0 to the random time tau and
     * returns the value of the time tau when the path is stopped.
     * Default implementation: calls to timeStep.</p>
     *
     * @param tau random time when path is stopped.
     */
    virtual int newPathSegment(StoppingTime* tau){  return newPathSegment(0,tau);  }
	
	
	/** <p>Continues a path existing on [0,t] from time t to the horizon, 
     * that is, computes \f$X_s, s=t+1,\dots,T\f$ from \f$X_u, u\leq t\f$
     * (branching a path at time t). </p>
     *
     * Default implementation: calls to timeStep. 
     * This is not the most efficient method, override if fastest 
     * path computation is a concern.
     *
     * @param t time of branching.
     */
    virtual void newPathBranch(int t){ newPathSegment(t,T); }
	

    /** <p>Computes a new path from time t=0.</p>
     */
    virtual void newPath(){ newPathBranch(0); }
    

// SAMPLING
	
    /** <p>Computes the random object \f$X_\tau\f$, that is, X=this sampled at the
     * stopping time \f$\tau\f$ conditioned on the state of the stochastic
	 * process X at time t.</p>
     *
	 * @param t current discrete time.
     * @param tau stopping time (>=t) at which X is sampled.
     */
    RandomObject<RangeType,ScalarType>* operator()(int t, StoppingTime* tau)
    {
        return new SampledAt(this,t,tau);
	} 
    
    /** <p>Computes the random object \f$X_\tau\f$, that is, X=this sampled at the
     * stopping time \f$\tau\f$ as seen from time t=0 (no conditioning).</p>
     *
     * @param tau stopping time at which X is sampled.
     */
    RandomObject<RangeType,ScalarType>* operator()(StoppingTime* tau)
    {
        return new SampledAt(this,t,tau);
	} 
	
private:
	
	/** <p>The random object \f$X_\tau\f$, that is, X=this sampled at the
     * stopping time \f$\tau\f$ conditioned on the state of the stochastic
	 * process X at time t.</p>
     */
	class SampledAt : public RandomObject<RangeType,ScalarType> {
		
		StochasticProcess X;
		StoppingTime tau;
		int t;
		
	 public:
		
		SampledAt(StochasticProcess* Y, int s, StoppingTime* psi) : 
		RandomObject<RangeType,ScalarType>(Y->getDimension()), X(Y), tau(psi), t(s) 
		{  }
		RangeType nextValue(){  return X->currentPath(X->newPathSegment(t,tau)); }
		
	}; // end SampledAt
    
    
}; //end StochasticProcess



/********************************************************************************************
 *
 *                    VECTOR VALUED STOCHASTIC PROCESS
 * 
 ********************************************************************************************/

typedef StochasticProcess< RealVector, Real > VectorProcess;


/** <p>Vector valued process adapted to a Brownian filtration. Only the following aspect 
 *  is abstracted: the increments of such a process are likely to be computed by 
 *  transforming a vector Z with standard normal components (normalized increments of
 *  a Wiener process ie. a Brownian motion).</p> 
 *
 *  <p>A series of such vectors is needed to drive a path of the process. This class allocates
 *  space for these vectors and computes the standard normal deviates by applying the
 *  inverse cumulative normal distribution to uniform deviates in (0,1).
 *  Two types of uniform deviates are used:
 *  <ul>
 *      <li>Pseudo uniform deviates based on the Mersenne Twister 
            (Monte Carlo, "MC", this is the default).</li>
 *      <li>Quasi "uniform" deviates from the Sobol low discrepancy sequence
 *          (Quasi Monte Carlo, "QMC").</li>
 *  </ul>
 * The class provides a mechanism to switch between the two and overrides all the path 
 * generation methods in the class StochasticProcess.</p>
 *
 * @author  Michael J. Meyer 
 */
class BrownianVectorProcess : public VectorProcess {

protected:
	
	// rows addressed as Real*
	Array1D< RealVector* > path;      // *(path[t]) is the state of the process at time t.
	RealMatrix Z;                     // the row Z[t][] is the vector driving the time step t->t+1.	
	StochasticGenerator* SG;          // the generator computing the Z-matrix
	
public:
	
	/** The state of the process at time t.*/
	RealVector& currentPath(int t) { return *(path[t]); }

	/** Start the process off at the point x: X(0)=x. */
	void setStart(const RealVector& x){ currentPath(0)=x; }
 
// CONSTRUCTOR
	
	/** @param dim dimension.
	 *  @param T   number of time steps to the horizon.
	 */
	BrownianVectorProcess(int dim, int T) : VectorProcess(dim,T),
	path(T+1), Z(T,dim,0,0), 
	SG(new MonteCarloVectorDriver(dim))
    {   
		for(int t=0;t<=T;t++) path[t]=new RealVector(dim);
	}
	
	~BrownianVectorProcess()
	{  
		for(int t=0;t<=T;t++) delete path[t];
		delete SG; 
	}


// PATH GENERATION
		

    virtual void newPathSegment(int t, int s)
	{  
		SG->newWienerIncrements(t,s,Z);
		for(int u=t;u<s;u++)timeStep(u);  
	}
	
    /** <p>Computes a new path segment from time 0 to time t.
     * Default implementation: calls to timeStep.</p>
     *
     * @param t time when path is stopped.
     */
    virtual void newPathSegment(int t){  newPathSegment(0,t);  }
    
      
    virtual int newPathSegment(int t, StoppingTime* tau)
    {
        SG->newWienerIncrements(t,T,Z);
		int s=t;
        while(!(tau->stop(s))){ timeStep(s); s++; }
        return s;
     }
	 
	/** <p>Computes a new path segment from time t=0 to the random time \f$\tau\f$ and
     * returns the value of the time tau when the path is stopped.
     * Default implementation: calls to timeStep.</p>
     *
     * @param tau random time when path is stopped.
     */
    virtual int newPathSegment(StoppingTime* tau){  return newPathSegment(0,tau);  }
	
	/** Prints the current path up to time t as a sequence of vectors.
	 */
	void printCurrentPath(int t)
    {
		cout << "\n\nBrownian vector process, current path:\n";
		for(int s=0;s<=t;s++)
		cout << "\nt = " << s << ", X(t):" << currentPath(s);
	}
	
	 
// switching between Monte Carlo and Quasi Monte Carlo dynamics
	 
	/** Switches to quasi random dynamics based on Sobol sequence.
     */
    void switchToQMC() 
	{  
		if(SG) delete SG;
		SG = new SobolVectorDriver(dim,T);
	}
	
	
    /** Switches to pseudo random dynamics based on Mersenne Twister.
     */
    void switchToMC() 
	{ 
		if(SG) delete SG;
		SG = new MonteCarloVectorDriver(dim);
	}
	 
	 
 }; // BrownianVectorProcess
     
    
   


/********************************************************************************************
 *
 *                    REAL STOCHASTIC PROCESS
 * 
 ********************************************************************************************/

		
typedef StochasticProcess<Real,Real> ScalarProcess;


/** Similar to {@link BrownianVectorProcess} but in dimension one treated as a 
 *  separate case to avoid awkward syntax.
 *
 * @author  Michael J. Meyer 
 */
class BrownianScalarProcess : public ScalarProcess {
	
protected:
	
	RealVector path;          // path[t] is he state of the process at time t
	RealVector  Z;            // the row Z[t] is the standard normal increment driving 
	                          // the time step t->t+1.
	StochasticGenerator* SG;  // the generator computing the Z-vector
	
public:
	
	/** The path array. */
	RealVector& getPath() { return path; }
			
	/** The state of the process at time t.*/
	Real& currentPath(int t) { return path[t]; }
	
// Constructor
	
	/** @param T   number of time steps to the horizon.
	 */
	BrownianScalarProcess(int T) : ScalarProcess(T),
	path(T+1),
	Z(T), 
	SG(new MonteCarloScalarDriver())
    {   }
	
	~BrownianScalarProcess(){ delete SG; }

// the new path generation methods
	

    virtual void newPathSegment(int t, int s)
	{  
		SG->newWienerIncrements(t,s,Z);
		for(int u=t;u<s;u++)timeStep(u);  
	}
    
      
    virtual int newPathSegment(int t, StoppingTime* tau)
    {
        SG->newWienerIncrements(t,T,Z);
		int s=t;
        while(!(tau->stop(s))){ timeStep(s); s++; }
        return s;
     }
	 
// switching between Monte Carlo and Quasi Monte Carlo dynamics
	 
	/** Switches to quasi random dynamics based on Sobol sequence.
     */
    void switchToQMC() 
	{  
		if(SG) delete SG;
		SG = new SobolScalarDriver(T);
	}
	
	
    /** Switches to pseudo random dynamics based on Mersenne Twister.
     */
    void switchToMC() 
	{ 
		if(SG) delete SG;
		SG = new MonteCarloScalarDriver();
	}
	 
 }; // BrownianScalarProcess
     
    
   
		


/********************************************************************************************
 *
 *                              PATH FUNCTIONALS
 * 
 ********************************************************************************************/
 
/** <p>A path functional H is a random object (random variable with arbitrary not necessarily
 *  scalar range) which is a deterministic function \f$H=h(X)\f$ of the path of a stochastic 
 *  process X. In this context there is a notion of time and of the information \f$F_t\f$ available 
 *  at time t. This information is the observed path of the process up to time t.
 *  Thus a path functional can be conditioned on the state of the process at time t by 
 *  evaluating it along paths which branch from the observed path at time t (continue this
 *  path from time t to the horizon).
 *
 * <p>The class PathFunctional implements this conditioning. The method 
 *  {@link #conditionedAtTime(int t)} returns a random object \f$H_t\f$ which is the path 
 *  functional H conditioned on the state of the process at time t. 
 *  Conditional expectations \f$E_t[H]=E[H|F_t]\f$ then become ordinary 
 *  expectations \f$E[H_t]\f$. These conditional expectations are implemented, all other 
 *  methods delivered by the class RandomObject are also available using the method
 *  <code>conditionedAtTime()</code>.
 *
 * <p>It is assumed that paths of the process X are computed as a series of time steps
 * of equal size dt. With this time becomes an integer variable and discrete time t corresponds 
 * to continuous time t*dt. In other words discrete time t is the number of the time step.
 *
 * <p>The template parameters RangeType, ScalarType have the same meaning as in the case of 
 * a {@link RandomObject}. The functional is RangeType valued and RangeType is a vector type with
 * dim components of type ScalarType. The ScalarType only comes into play if covariances are 
 * computed, the dimension dim of the RangeType only comes into play if the covariance matrix 
 * is computed. The default values can be used otherwise. 
 *
 * <p>ProcessRangeType is the type of values of of the underlying stochastic process. 
 * The process and the path functional are assumed to have the same ScalarType.
 *
 * <p>A PathFunctional is fully specified as soon as the template parameters and the
 * method {@link #valueAlongCurrentPath} are defined.
 *
 * @param ProcessRangeType range of underlying stochastic process (possibly vectors).
 * @param RangeType range of functional (possibly vectors).
 * @param ScalarType common scalar type (type of vector components).
 *
 * @author  Michael J. Meyer
 */
template<typename ProcessRangeType=Real,         
         typename RangeType=Real, 
		 typename ScalarType=Real>
class PathFunctional {
	 
protected:
	 
	 // the process of which it is a functional
	 StochasticProcess<ProcessRangeType,ScalarType>* X;        
	 int dim;      // dimension of the RangeType of the functional.
	 
public: 
	
	 /** Dimension of functional RangeType. */
	 int getDimension(){ return dim; }
	 
	 /** The underlying stochastic process. */
	 StochasticProcess<ProcessRangeType,ScalarType>* getProcess(){ return X; }
	 
	 /** @param Y the underlying process. 
	  *  @param d dimension of RangeType of functional.
	  */
	 PathFunctional
	 (StochasticProcess<ProcessRangeType,ScalarType>* Y, int d=1) : 
	 X(Y), dim(d) {  }
	 
	 /** The value of the functional along the current path of the process.
	  */
	 virtual RangeType valueAlongCurrentPath() = 0;
	 
	 
	 /** <p>The path functional conditioned on the state of the process at time t
	  *  as a random object with values in RangeType.
	  *
	  *  <p>If RangeType==ScalarType==Real, this function returns a pointer
	  *  to RandomVariable and if RangeType==RealVector, ScalarType==Real, 
	  *  this function returns a pointer to RandomVector. Here vector is our own
	  *  light weigth vector class not std::vector.
	  */
	 RandomObject<RangeType,ScalarType>* 
	 conditionedAtTime(int t){ return new H_t(this,t); }
	 
	 /** Monte Carlo expectation of the functional conditioned on the state
	  *  of the process at time t.
	  * @param t time of conditioning.
	  * @param nPath number of paths for evaluation.
	  * @param progressReport periodically report time left to the console.
	  * @param message string identifying computation in case of a progress report.
	  */
	 RangeType conditionalExpectation
	 (int t, int nPath, bool progressReport=false, string message=" ")
     {
		 if(progressReport)
		 return conditionedAtTime(t)->expectation(nPath,message);
		 // otherwise
		 return conditionedAtTime(t)->expectation(nPath);
	 }

	 
	 
	 /** Monte Carlo conditional mean ([0]) and variance ([1]).
	  * @param t time of conditioning.
	  * @param nPath number of paths for evaluation.
	  * @param progressReport periodically reports time left to the console.
	  * @param message string identifying computation in case of a progress report.
	  */
	 RangeType* conditionalMeanAndVariance
	 (int t, int nPath, bool progressReport=false, string message=" ")
     {
		 if(progressReport)
		 return conditionedAtTime(t)->meanAndVariance(nPath,message);
		 // otherwise
		 return conditionedAtTime(t)->meanAndVariance(nPath);
	 }
	 
	 
	 
private: 
	 
	 /** The path functional conditioned on the state of the process at time t.
	  *  ScalarType and dimension of RangeType irrelevant, since no covariances
	  *  will be computed.
	  */
	 class H_t : public RandomObject<RangeType> {
		 
		 PathFunctional<ProcessRangeType,RangeType>* H;
		 int t;                      // time of conditioning
		 
	 public:
		 
		 /** @param Y underlying process.
		  *  @param s time of conditioning, default = 0.
		  */
		 H_t(PathFunctional<ProcessRangeType,RangeType>* G, int s=0) : 
		 RandomObject<RangeType>(), H(G), t(s)
	     {   }
		 
		 /** Next observation of the functional conditioned on the state of 
		  *  the process at time t.
		  */
		 RangeType nextValue()
	     {
			 StochasticProcess<ProcessRangeType,ScalarType>* X=H->getProcess();
			 X->newPathBranch(t);
			 return H->valueAlongCurrentPath();
		 }
	 }; // end H_t
	 
 }; // end PathFunctional
 
 
 /** A simple path functional H=h(X) of a vector process X of dimension at least two:
  *  h(X)=X_1(T)+X_2(T). Here ProcessRangeType==RealVector, and defaults for the
  *  functional: RangeType==ScalarType==Real.
  */
class SumFunctional : public PathFunctional< RealVector > {
			 
	VectorProcess* X; int T;		
public:
	SumFunctional(VectorProcess* Y, int T_oo) : 
	PathFunctional< RealVector >(Y), X(Y), T(T_oo) {  }
			 
	Real valueAlongCurrentPath()
    {
	    RealVector& x=X->currentPath(T);
		return x[0]+x[1];
	}
}; // end SumFunctional

MTGL_END_NAMESPACE(Martingale)

#endif

