/***************************************************************************
                          BasisFunctions.cpp  -  description
                             -------------------
    begin                : Wed Jun 16 2004
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


#include "BasisFunctions.h"
#include <cmath>
using std::sin;
using std::cos;


RealArray1D
LegendreBasis::
values(Real x, int m)
{
	RealArray1D P(m+1);
	int r;
	P[0]=1.0;
	P[1]=sqrt(3)*x;
	for(int n=1;n<m;n++){

		r=2*n+1;
		P[n+1]=sqrt(r*(r+2))*x*P[n];
		P[n+1]-=n*sqrt(1.0+4.0/(r-2))*P[n-1];
		P[n+1]/=(n+1);
	}
	return P;
}


RealArray1D
FourierBasis::
values(Real x, int m)
{
   assert(m>=0);
   RealArray1D P(m+1);
   int k=0;
   // recall that we are using nomalized Lebesgue measure 0.5*dt on [-1,+1]
   // (affects the empirical coefficients)
   for(k=0;k<=m;k+=2) P[k]=2*cos(k*PI*x);
   for(k=1;k<=m;k+=2) P[k]=2*sin(k*PI*x);
   
	return P;
} 