/** WARANTY NOTICE AND COPYRIGHT
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



#include "FinMath.h"
#include <cmath>
#include <iostream>
using namespace Martingale;



/** Standalone methods to compute functions or solve eqations
 *  useful in basic financial mathematics.
 *
 * @author Michael J. Meyer
 */

    
 
 
/******************************************************************************* 

                           MULTINORMAL DENSITY
                                     
*******************************************************************************/ 
 
 /** The density of the standard multinormal distribution N(0,I).
  */
 Real 
 FinMath::
 multiNormalDensity(int dim, Real* z)
 {
     Real f=2.506628274631,    // sqrt(2pi)
          q=1.0;
     
     for(int k=0;k<dim;k++)q*=(f*std::exp(-z[k]*z[k]/2));
     return q;
 }

/******************************************************************************************************* 

                                       N(x), h_+, h_-
                                     
*******************************************************************************************************/ 

/**<p>The cumulative distribution function N(x)=Prob(X<=x) of a standard normal 
 *  variable X derived from the error function erf(x)=Prob(-sqrt(2)x<=X<=sqrt(2)x).
 * </p>
 *
 * @param x Any real number.
 */
  
 Real FinMath::N(Real x)
 { 
	 Real sqrt2=1.41421356;
	 return 0.5*(1.0+erf(x/sqrt2)); 
 }
  
 /** <p>The quantity d_+ in Margrabe's formula for the option to 
  * exchange assets S_1, S_2 (receive S_1 for kS_2). In this context:</p>
  *
  *<center> d_+=log(Q/k)/Sigma+Sigma/2 with Q=S_1(t)/S_2(t)</center>
  * 
  * <p>and Sigma is the aggregate volatility of log(Q) from current
  * time t to the horizon T, that is, the quadratic variation</p>
  *
  * <center>Sigma = &#139 log(Q) &#155_[t,T].</center>
  *
  * <p>If Q has constant annual volatility sigma this becomes</p>
  *
  * <center>Sigma=sigma*sqrt(tau)</center> 
  *
  * <p>where tau=T-t is time to expiry. In the case of the European call on 
  * S_1 with strike k the asset S_2 
  * is the zero coupon bond maturing at call expiration.</p>
  *
  * @param Q The quotient S_1(t)/S_2(t)
  * @param k Ratio of exchange between assets (receive S_1 for kS_2).
  * @param Sigma Quadratic variation &#139 log(Q) &#155_[t,T].
  */
 Real 
 FinMath::
 d_plus(Real Q, Real k, Real Sigma)
 { 
    return log(Q/k)/Sigma + Sigma/2; 
 }
 
 /** <p>The quantity d_- in Margrabe's formula for the option to 
  *  exchange assets. </p>
  *  <p>See {@link #d_plus} and replace "+Sigma/2" with "-Sigma/2".</p>
  **/
 Real 
 FinMath::
 d_minus(Real Q, Real k, Real Sigma)
 { 
    return log(Q/k)/Sigma - Sigma/2; 
 } 
   


/*******************************************************************************

                               BLACK-SCHOLES FUNCTION
                                     
*******************************************************************************/


 /**<p>Computes the function QN(h_+)-kN(h_-).</p>
 *
 * <p>Useful in application of the Black-Scholes formula or Margrabe's formula 
 * to price calls or options to exchange assets. For the option to 
 * exchange assets S_1,S_2 we have Q=S_1(t)/S_2(t) and</p>
 *
 * <center>Price(t,S_1,S_2)=QN(h_+)-kN(h_-)</center>
 *
 * <p> See {@link #d_plus} for information on the parameter Sigma.
 * In the special case of a call on S we have S_1=S, S_2 is the zero coupon 
 * bond maturing at call expiry, Q the forward price of S and the forward price 
 * of the call can be written as</p>
 *
 * <p><center>Call.ForwardPrice(t,S)=S_2*(QN(h_+)-kN(h_-))</center></p>
 *
 * @param Q The quotient S_1(t)/S_2(t)
 * @param k Ratio of exchange between assets (receive S_1 for kS_2).
 * @param Sigma Quadratic variation &#139 log(Q) &#155_[t,T].
 */    
 Real 
 FinMath::
 blackScholesFunction(Real Q, Real k, Real Sigma)
 {
   return Q*N(d_plus(Q,k,Sigma))-k*N(d_minus(Q,k,Sigma)); 
 } 
   
 /** <p>Derivative of the Black-Scholes function QN(h_+)-kN(h_-) with 
  *  respect to Sigma.</p>
  *
  * @param Q The quotient S_1(t)/S_2(t)
  * @param k Ratio of exchange between assets (receive S_1 for kS_2).
  * @param Sigma Quadratic variation &#139log(Q) &#155_[t,T].
  */    
 
 Real 
 FinMath::
 dBSF(Real Q, Real k, Real Sigma)
 {
   Real Sqrt2Pi=2.506628274631;
   return 2*k*exp(-d_minus(Q,k,Sigma)*d_minus(Q,k,Sigma)/2)/Sqrt2Pi; 
 } 
 
 
  /** <p><i>Discounted</i> Black-Scholes call price (note that current time may
   *  not be zero). Both the current time <code>t</code> and the time 
   *  <code>T</code> to expiry (from time zero) are needed since the formula 
   *  will be applied at time other than time zero and <code>S</code></p>
   *
   * @param S the <i>forward</i> asset price at expiry <code<T</code>.
   * @param K strike price
   * @param tau time to expiration
   * @param sigma annual asset volatility
   * @param B risk free bond B_T at expiry 
   */    
 Real 
 FinMath::
 bsDiscountedCallPrice(Real S, Real K, Real tau, Real sigma, Real B)
 {
   Real Sigma=sigma*sqrt(tau);                
   return blackScholesFunction(S,K,Sigma)/B; 
 } 
  
  /** <p><i>Discounted</i> Black-Scholes put price (note that current time may
   *  not be zero). Here we need the interest rate, the current time 
   *  <code>t</code> and the time <code>T</code> of expiration explicitly.
   *  The formula will be applied at times other than time zero.</p>
   *
   * @param S the <i>forward</i> asset price at expiry <code>T</code>.
   * @param K strike price
   * @param tau time to expiry
   * @param sigma annual asset volatility
   * @param B risk free bond B_T at expiry
   */    
 Real 
 FinMath::
 bsDiscountedPutPrice(Real S, Real K, Real tau, Real sigma, Real B)
 {
   Real Sigma=sigma*sqrt(tau);
   return (K+blackScholesFunction(S,K,Sigma)-S)/B; 
 }
   

/*******************************************************************************

                         IMPLIED SIGMAs
                                     
*******************************************************************************/



/** <p>Given \f$y,Q,k>0\f$ solves the equation 
 *  \f$blackScholesFunction(Q,k,\Sigma)=y\f$ for 
 *  \f$\Sigma\geq0\f$ using continued bisection.</p>
 *
 * <p> A solution \f$\Sigma=\Sigma(Q,k,y)\f$ exists if and only if \f$y<Q\f$.
 * To be used to compute implied volatilities from call forward prices y. 
 * Note that the solution \f$Sigma\f$ is not the implied annual volatility \f$sigma\f$, 
 * indeed these are related as \f$Sigma=sigma*\sqrt{\tau}\f$, where 
 * \f$\tau=T-t\f$ is time to expiry.</p>
 */
// Search restricted to [0.00001,1000], smaller or larger solutions cut off
// to interval endpoints (irrelevant anyway). The blackScholesFunction
// f(Q,k,Sigma) increases with Sigma. The algorithm preserves the relation<br> 
// f(Q,k,a) &lt=y&lt=f(Q,k,b).
Real 
FinMath::
blackImpliedAggregateCallVolatility(Real Q, Real k, Real y)
{
     Real a=0.00001,c,b=1000.0,error;
     int n=0;         
         
     if(blackScholesFunction(Q,k,a)>y)return a;
     if(y>=Q){
	  
	     std::cout << "\nFinMath::blackImpliedAggregateCallVolatility(): "
	               << "no solution.\n Returning 0.\n";
	     return 0.0;
     }
   
      while(n<100){
			
           c=0.5*(a+b);
           error=blackScholesFunction(Q,k,c)-y;
           if(fabs(error)<0.00000001) 
               return c;
           else if(error>0) b=c; else a=c;
           n++;
       } // end while
    
       return a;
  
} //end blackImpliedAggregateCallVolatility
    

 






