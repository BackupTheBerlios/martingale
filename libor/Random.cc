
#include "Random.h"
#include "Utils.h"
#include <cmath>


using std::erf;
using std::cerr;
using std::endl;

MTGL_BEGIN_NAMESPACE(Martingale)
MTGL_BEGIN_NAMESPACE(Random)


// static initialization
unsigned long 
MersenneTwister::mag01[2]={0x0, MATRIX_A};

MersenneTwister MT19937;


MersenneTwister::
MersenneTwister(unsigned long seed) : mt(new unsigned long[MT_N]), mti(MT_N)
{	
    // state vector intialization
    for (int i=0;i<MT_N;i++) {
		
        mt[i] = seed & 0xffff0000;
        seed = 69069 * seed + 1;
        mt[i] |= (seed & 0xffff0000) >> 16;
        seed = 69069 * seed + 1;
    }
}


Real 
MersenneTwister::
u01()
{
    unsigned long y;

    if (mti >= MT_N) { /* generate N words at one time */
        int kk;

        for (kk=0;kk<MT_D;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<MT_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk-MT_D] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[MT_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MT_N-1] = mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

	/* reals: (0,1)-interval */
    return ( ((double)y + 1.0) * 2.3283064359965952e-10 ); 
    /* return y; for integer generation */
} 


Real 
U01(){ return MT19937.u01(); } 


int 
Uint(int n){ return (int) (n*U01()); }


int 
sign()
{ 
     static int S1=0,S2=0,S3=0;           // state saved accross function calls
	 int next=(S3*S2+19*S1+22)%331;       // next random integer
     S1=S2; S2=S3; S3=next;
     return ((next%2)==0)? 1:-1;  
} //end Sign
	
   

Real 
N(Real x){ return 0.5*(1.0+erf(x/sqrt(2.0))); }
	


Real 
nInverse(Real x)
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
   

Real 
sTN(){ return nInverse(U01()); }




MTGL_END_NAMESPACE(Random)
MTGL_END_NAMESPACE(Martingale)



   

 
