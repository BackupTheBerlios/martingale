/***************************************************************************
                          SampleFunctions.h  -  description
                             -------------------
    begin                : Tue Jun 15 2004
    copyright            : (C) 2004 by Michael J. Meyer
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

#ifndef FunctionExamples_h
#define FunctionExamples_h

#include "TypedefsMacros.h"
#include "Array.h"
using namespace Martingale;



/** \f$f(t)=sin(2\pi t)\f$.
 **/
Real f_0(Real t);

/** \f$f(t)=20t^{1/3}exp(-6t)\f$.
 **/
Real f_1(Real t);

/** \f$f(t)=|t|\f$.
 **/
Real f_2(Real t);

/** \f$f(t)=|t|^{1/3}\f$.
 **/
Real f_3(Real t);

/** Sets regression parameters from user input.
 *
 * @param N index of last basis function.
 * @param n number of data points-1.
 * @param s array of data point x-coordinates.
 * @param y array of y_j=f(s_j).
 * @returns pointer to function to be expanded.
 *
 **/
RealFunction expansionDialog
(int& N, int& n, RealArray1D& s, RealArray1D& y);


#endif


