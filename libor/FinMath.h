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

#ifndef martingale_finmath_h    
#define martingale_finmath_h


#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(FinMath)	



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
Real multiNormalDensity(int dim, Real* z);


/******************************************************************************************************* 

                                       N(x), h_+, h_-
                                     
*******************************************************************************************************/ 

/**<p>The cumulative distribution function N(x)=Prob(X<=x) of a standard normal 
 *  variable X derived from the error function 
 * \f[erf(x)=Prob(-\sqrt{2}\,x\leq X\leq\sqrt{2}\,x).\f]
 * </p>
 *
 * @param x Any real number.
 */
Real N(Real x);

  
 /** <p>The quantity \f$d_+\f$ in Margrabe's formula for the option to 
  * exchange assets S_1, S_2 (receive S_1 for kS_2). In this context:</p>
  *
  * \f[d_+=log(Q/k)/\Sigma+\Sigma/2\f]
  *
  * where \f$Q=S_1(t)/S_2(t)\f$ and \f$\Sigma\f$ is the aggregate volatility of 
  * log(Q) from current time t to the horizon T. If Q has constant annual volatility 
  * \f$sigma\f$ this becomes
  *
  * \f[Sigma=\sigma\sqrt{T-t}.\f] 
  *
  * In the case of the European call on S_1 with strike k the asset S_2 
  * is the zero coupon bond maturing at call expiration.</p>
  *
  * @param Q The quotient S_1(t)/S_2(t)
  * @param k Ratio of exchange between assets (receive S_1 for kS_2).
  * @param Sigma aggregate volatility of Q on [t,T].
  */
Real d_plus(Real Q, Real k, Real Sigma);

 
 /** <p>The quantity \f$d_-\f$ in Margrabe's formula for the option to 
  * exchange assets S_1, S_2 (receive S_1 for kS_2). In this context:</p>
  *
  * \f[d_-=log(Q/k)/\Sigma-\Sigma/2\f]
  *
  * where \f$Q=S_1(t)/S_2(t)\f$ and \f$\Sigma\f$ is the aggregate volatility of 
  * log(Q) from current time t to the horizon T. If Q has constant annual volatility 
  * \f$sigma\f$ this becomes
  *
  * \f[Sigma=\sigma\sqrt{T-t}.\f] 
  *
  * In the case of the European call on S_1 with strike k the asset S_2 
  * is the zero coupon bond maturing at call expiration.</p>
  *
  * @param Q The quotient S_1(t)/S_2(t)
  * @param k Ratio of exchange between assets (receive S_1 for kS_2).
  * @param Sigma aggregate volatility of Q on [t,T].
  */
Real d_minus(Real Q, Real k, Real Sigma);



/*******************************************************************************

                               BLACK-SCHOLES FUNCTION
                                     
*******************************************************************************/


/**<p>Computes the function \f$Q*N(d_+)-k*N(d_-)\f$ with \f$d_\pm\f$ as in
 * {@link #d_plus}, {@link #d_minus}.</p>
 *
 * <p>Useful in application of the Black-Scholes formula or Margrabe's formula 
 * to price calls or options to exchange assets. For the option to 
 * exchange assets S_1,S_2 we have Q=S_1(t)/S_2(t) and</p>
 *
 * \f[Price(t,S_1,S_2)=Q*N(d_+)-k*N(d_-).\f]
 *
 * In the special case of a call on S we have S_1=S, S_2 is the zero coupon 
 * bond maturing at call expiry, Q the forward price of S and the forward price 
 * of the call can be written as
 *
 * \f[Call.ForwardPrice(t,S)=Q*N(d_+)-k*N(d_-)=.\f]
 *
 * @param Q The quotient S_1(t)/S_2(t)
 * @param k Ratio of exchange between assets (receive S_1 for kS_2).
 * @param Sigma aggregate volatility of Q on [t,T].
 */    
Real blackScholesFunction(Real Q, Real k, Real Sigma);

   
 /** <p>Derivative of the Black-Scholes function Q*N(d_+)-k*N(d_-) with 
  *  respect to Sigma.</p>
  *
  * @param Q The quotient S_1(t)/S_2(t)
  * @param k Ratio of exchange between assets (receive S_1 for kS_2).
  * @param Sigma aggregate volatility of Q on [t,T].
  */    
Real dBSF(Real Q, Real k, Real Sigma);

 
 
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
Real bsDiscountedCallPrice(Real S, Real K, Real tau, Real sigma, Real B);

  
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
Real bsDiscountedPutPrice(Real S, Real K, Real tau, Real sigma, Real B);
   

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
Real blackImpliedAggregateCallVolatility(Real Q, Real k, Real y);
    




MTGL_END_NAMESPACE(FinMath)
MTGL_END_NAMESPACE(Martingale)	

#endif