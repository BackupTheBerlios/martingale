/***************************************************************************
                          BasisFunctions.h  -  description
                             -------------------
    begin                : Tue Jun 15 2004
    copyright            : (C) 2004 by 
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


 #include "Matrix.h"
 #include <string>
 using std::string;
 using namespace Martingale;

 
/** Basis functions on [-1,+1]. Usually orthonormal in \f$L^2([-1,+1], dt)\f$.
 *  Note that we are using <i>normalized</i> Lebesgue measure \f$d\mu=(1/2)dt\f$
 *  on [-1+1].
 **/
class BasisFunctions {

public:

   /** Computes the sequence \f$(\psi_0(x),\dots,\psi_m(x))\f$ of values of the
    *  basis functions.
    **/
   virtual RealArray1D values(Real x, int m)=0;

   /** Designation of basis.
    **/
   virtual string name()=0;

};


/** Orthonormal Legendre polynomials on [-1+1] with normalized Lebesgue measure. **/
class LegendreBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   string name(){ return "orthonormal Legendre basis"; }
};


/** Orthonormal Fourier basis on [-1+1] with normalized Lebesgue measure. **/
class FourierBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   string name(){ return "orthonormal Fourier basis"; }
};
