/* WARANTY NOTICE AND COPYRIGHTThis program is free software; you can redistribute it and/ormodify it under the terms of the GNU General Public Licenseas published by the Free Software Foundation; either version 2of the License, or (at your option) any later version.This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty ofMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See theGNU General Public License for more details.You should have received a copy of the GNU General Public Licensealong with this program; if not, write to the Free SoftwareFoundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.Copyright (C) Michael J. Meyermatmjm@mindspring.comspyqqqdia@yahoo.com*/


#include <math.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include "Random.h"
#include "Utils.h"


MTGL_BEGIN_NAMESPACE(Martingale)


// Mersenne twister mt19937 based uniform random number generator
boost::uniform_01<boost::mt19937,Real> uniform01(*(new boost::mt19937()));
   
// random draw from {1,-1} 
int Random::sign()
{ 
     static int S1=0,S2=0,S3=0;           // state saved accross function calls
	 int next=(S3*S2+19*S1+22)%331;       // next random integer
     S1=S2; S2=S3; S3=next;
     return ((next%2)==0)? 1:-1;  
} //end Sign
	
   
//  cumulative normal distribution function 
Real Random::N(Real x){ return 0.5*(1.0+erf(x/sqrt(2.0))); }
	

//  inverse normal cumulative distribution function 
Real Random::nInverse(Real x)
{
    //const Real SQRT_TWO_PI=2.5066282746310;
    
    // Coefficients for the rational approximation.
    const Real
            e_1 = -3.969683028665376e+01,
            e_2 =  2.209460984245205e+02,
            e_3 = -2.759285104469687e+02,
            e_4 =  1.383577518672690e+02,
            e_5 = -3.066479806614716e+01,
            e_6 =  2.506628277459239e+00;
    
    const Real
            f_1 = -5.447609879822406e+01,
            f_2 =  1.615858368580409e+02,
            f_3 = -1.556989798598866e+02,
            f_4 =  6.680131188771972e+01,
            f_5 = -1.328068155288572e+01;
    
    const Real
            g_1 = -7.784894002430293e-03,
            g_2 = -3.223964580411365e-01,
            g_3 = -2.400758277161838e+00,
            g_4 = -2.549732539343734e+00,
            g_5 =  4.374664141464968e+00,
            g_6 =  2.938163982698783e+00;
    
    const Real
            h_1 =  7.784695709041462e-03,
            h_2 =  3.224671290700398e-01,
            h_3 =  2.445134137142996e+00,
            h_4 =  3.754408661907416e+00;
    
    // split domain: ]-oo,x_l[, [x_l,x_u], ]x_u,+oo[
    const Real
            x_l   = 0.02425,
            x_u  = 1.0 - x_l;

    Real z, r;

    // in region ]-oo,x_l[
    if( x < x_l )
    {  
          z = sqrt(-2.0*log(x));
          z = (((((g_1*z+g_2)*z+g_3)*z+g_4)*z+g_5)*z+g_6) / 
              ((((h_1*z+h_2)*z+h_3)*z+h_4)*z+1.0);
    }
  
    // in region [x_l,x_u] 
    else if( x <= x_u )
    {  
         z = x - 0.5; r = z*z;
         z = (((((e_1*r+e_2)*r+e_3)*r+e_4)*r+e_5)*r+e_6)*z / 
             (((((f_1*r+f_2)*r+f_3)*r+f_4)*r+f_5)*r+1.0);
    }
 
    // in region ]x_u,+oo[
    else 
    {  
          z = sqrt(-2.0*log(1.0-x));
          z = -(((((g_1*z+g_2)*z+g_3)*z+g_4)*z+g_5)*z+g_6) /  
               ((((h_1*z+h_2)*z+h_3)*z+h_4)*z+1.0);
    }

    // Now |relative error| < 1.15e-9.  Do one iteration of Halley's third
    // order zero finder for full machine precision:
    /*
      r = (N(z) - x) * SQRT_TWO_PI * exp( 0.5 * z * z );	//	f(z)/df(z)
      z -= r/(1+0.5*z*r);
    */

    return z;

} // end nInverse
  


// simple uniform random number generator 
Real Random::U01()
{    
    static int R1=0,R2=0,R3=0;           // state saved accross function calls
    int next; Real f=1/331.0;
      
    for(int i=0;i<3;i++)
    {
        next=(R3*R2+19*R1+22)%331;        // next random integer
        R1=R2; R2=R3; R3=next;            // move state of the random integer sequence,
    }                                     // repeat three times to decorrelate adjacent U's
                                          // keep history three steps deep
    return f*(R3+f*(R2+f*R1));            // compute current U(n)
} //end U

   

// random integers in [0,n)
int Random::Uint(int n)
{
	return (int) n*U01();
}


// standard normal deviate 
Real Random::sTN(){ return nInverse(uniform01()); }



MTGL_END_NAMESPACE(Martingale)



   

 
