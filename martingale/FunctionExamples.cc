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


#include "FunctionExamples.h"
#include "Random.h"
#include <cmath>
#include <fstream>
#include <iomanip>    // setprecision

using std::ofstream;
using std::setprecision;
using std::sin;
using std::exp;
using std::log;
using namespace Martingale;


Real
f_0(Real t){ return sin(TWO_PI*t); }

Real
f_1(Real t){ return 5*t*exp(-9*t*t/2); }

Real
f_2(Real t){ return (t>0)? t:-t; }

Real
f_3(Real t){
   if(t==0) return 0.0;
   return (t>0)? exp(log(t)/3):exp(log(-t)/3);
}

RealFunction
expansionDialog
(int& N, int& n, RealArray1D& s, RealArray1D& y)
{
   RealFunction  f;
   cout << "Sample function f:" << endl
        << "f_0(t)=sin(2pi*t)...........[0]" << endl
        << "f_1(t)=5t*exp(-9t^2/2)......[1]" << endl
        << "f_2(t)=|t|..................[2]" << endl
        << "f_3(t)=|t|^{1/3}............[3]" << endl << endl
        << "Enter f=[0,1,2,3]=";
   int fnum; cin>>fnum;
   // *f is reference to pointer to sample function
   switch(fnum){
      case 1  : f=&f_1; break;
      case 2  : f=&f_2; break;
      case 3  : f=&f_3; break;
      default : f=&f_0;
   }
   cout << endl << "Enter number N+1 of basis functions, N=";
	cin>>N;
   cout << "Enter number n+1 of data points, n=";
	cin>>n;

	s.resize(n+1);
	y.resize(n+1);

	int random;
	cout << "Data points random (random=1) or evenly spaced (random=0), random=";
	cin>>random;

	// data abscissas s_j
	if(random==0)
	   for(int j=0;j<=n;j++) s[j]=-1.0+j*2.0/n;
   else
	   for(int j=0;j<=n;j++) s[j]=-1.0+2.0*Random::U01();

	int noisy;
	cout << "Function data noisy (noisy=1) or exact (noisy=0), noisy=";
	cin>>noisy;

	Real sigma=0.0;  // standard deviation of data noise
	if(noisy==1){

	   cout << "Standard deviation of noise, sigma=";
	   cin>>sigma;
	}
   // (*f) isreference to function pointer
	for(int j=0;j<=n;j++) y[j]=(*f)(s[j])+sigma*Random::sTN();

	// write log file
	ofstream lout("ExpansionLog.txt");
	lout << setprecision(3);
	lout << "Function expansion, session log." << endl << endl
	     << n+1 << " data points, ";
	if(random==1)
		lout << "randomly spaced." << endl;
	else
		lout << "evenly spaced." << endl;

	if(noisy==1)
		lout << "Function data noisy, "
	        << "standard deviation of noise" << sigma << endl;
	else
		lout << "Function data exact." << endl << endl;

   switch(fnum){
      case 1:  lout << "Function f(t)=5t*exp(-9t^2/2)"; break;
      case 2:  lout << "Function f(t)=|t|"; break;
      case 3:  lout << "Function f(t)=|t|^{1/3}"; break;
      default: lout << "Function f(t)=sin(2pi*t)";
   }

	lout << endl << endl 
        << "Session will write data files Expansion.dat, Function.dat\n"
	     << "for plotting with gnuplot, see doc/gnuplot_Readme.html";
	lout.close();

   return f;

} // ExpansionDialog






