/***************************************************************************
 *            gpr.cc
 *
 *  Sat May 29 16:54:15 2004
 *  Copyright  2004  cpp
 *  cpp@linux
 ****************************************************************************/

#include "gpr.h"
#include <fstream>
using std::ofstream;



GPR::
GPR(int dmax, RealArray1D& t, RealArray1D& w, BasisFunctions* bFcns):
regrType(GAUSSIAN),
N(dmax),
n(t.getDimension()-1),
s(t),
y(w),
basis(bFcns),
basis_name(basis->name()),
psi(N+1,n+1,0,0),
K(n+1), 
R(n+1),
mu(dmax+1),
a(dmax+1),
empCoeff(dmax+1)
{
   fullInitialization();
}



void
GPR::
fullInitialization()
{
   // intialize psi=(psi_k(s_j))_{0<=k<=N, 0<=j<=n}
	for(int j=0;j<=n;j++){

	   RealArray1D P=basisFunctionValues(s[j],N);
	   for(int k=0;k<=N;k++) psi(k,j)=P[k];
	}

	// allocate K=(K(s_i,s_j))
	for(int i=0;i<=n;i++)
	for(int j=i;j<=n;j++){

		Real sum=0.0;
		for(int k=0;k<=N;k++) sum+=psi(k,i)*psi(k,j);
		K(i,j)=sum;
	}
	LTRRealMatrix& rho=K.ltrRoot();
	R=rho;
	delete &rho;
	writeCoefficients();
}



Real
GPR::
expansion(Real t, int q)
{
	RealArray1D P=basisFunctionValues(t,q);
	Real f_q;
   if(regrType==GAUSSIAN)
	   for(int k=0;k<=q;k++) f_q+=a[k]*P[k];
   else
      // empirical regression
      for(int k=0;k<=q;k++) f_q+=empCoeff[k]*P[k];
	return f_q;
}


Real
GPR::
EA(int k)
{
  assert((0<=k)&&(k<=N));
	// compute the r[j]=R_{n+1,j}, 0<=j<=n,
	// equations 22,23 p11, gprs.ps
	RealArray1D r(n+1);
	r[0]=psi(k,0)/R(0,0);
	for(int j=1;j<=n;j++){

		Real sum=0.0;
		for(int i=0;i<j;i++)sum+=r[i]*R(j,i);

		r[j]=(psi(k,j)-sum)/R(j,j);
	}
  
  // the means mu_E[j]=E^P(E_j) of the evaluation functionals E_j
  // these are nonzero if P is not centered at zero
  RealArray1D mu_E(n+1);
  for(int j=0;j<=n;j++){

    Real sum=0.0;
    for(int k=0;k<=N;k++) sum+=mu[k]*psi(k,j);
    mu_E[j]=sum;
  }
     
	// compute the Z_j, 0<=j<=n
	RealArray1D Z(n+1);
	Z[0]=(y[0]-mu_E[0])/R(0,0);
	for(int j=1;j<=n;j++){
		
		Real sum=0.0;
		for(int k=0;k<j;k++) sum+=R(j,k)*Z[k];
			
		Z[j]=(y[j]-mu_E[j]-sum)/R(j,j);
	}

	// compute a_k=E(A_k), formula 20, p10, gprs.ps
  // note A_k has (unconditional) mean mu[k]
	Real a_k=mu[k];
	for(int j=0;j<=n;j++) a_k+=r[j]*Z[j];
	
	return a_k;
}



Real
GPR::
EC(int k)
{
   assert((0<=k)&&(k<=N));

	Real sum=0.0;
	for(int j=0;j<=n;j++) sum+=y[j]*psi(k,j);
	return sum/(n+1);
}


void
GPR::
writeCoefficients()
{
   for(int k=0;k<=N;k++){ a[k]=EA(k); empCoeff[k]=EC(k); }
}


void
GPR::
recomputeGaussianCoefficients()
{
   for(int k=0;k<=N;k++) a[k]=EA(k); 
}


string
GPR::
regressionType()
{
   if(regrType==GAUSSIAN) return "Gaussian";
   return "empirical";
}
			
		

void
GPR::
expansionData(RealFunction f, int q)
{		
	// write the data
	ofstream dout("FunctionData.txt");
	for(int j=0;j<=n;j++) dout << s[j] << "  " << y[j] << endl;
	dout.close();
	
	// write the expansions f_0,f_1,...,f_q evaluated at 801 points
	// to a file in gnuplot data format.
	ofstream fout("ExpansionData.txt");
	int m=800;
	for(int j=0;j<=m;j++){

	   Real t=-1.0+j*2.0/m;
	   fout << t << "  " << f(t) << "  ";
	   RealArray1D psi=basisFunctionValues(t,q);
	   Real sum=0.0;
	   for(int i=0;i<=q;i++){ sum+=a[i]*psi[i]; fout << sum << "  "; }
	   fout << endl;
   }
   fout.close();
}



void
GPR::
regressionExample()
{
   cout << "Gaussian regression, enter prameters."
        << endl << endl;

   // get center of Gaussian prior
   int ctr;
   Real r=0.0;
   cout << "\nMean of Gaussian prior P:" << endl
        << "origin....[0]" << endl
        << "empirical coefficients....[1]" << endl
        << "vector (r,r,...,r)....[2]"
        << "Default is [0]"
        << endl << endl
        << "Mean=[0,1,2]=";
   cin>>ctr;
   if (ctr==2){ cout << "Enter r=";  cin>>r; }

   int N;
   int n;
	RealArray1D s(1);
	RealArray1D y(1);            // resized in user dialogue
   // set parameters and select function
   RealFunction f=expansionDialog(N,n,s,y);
   BasisFunctions* bFcns=new FourierBasis();
   GPR gpr(N,s,y,bFcns);
   
   cout << endl
        << "Empirical [0] or Gaussian [1] regression? Regression=";
   int regr; cin>>regr;
   if(regr==0){ RegressionType rt=EMPIRICAL; gpr.setRegressionType(rt); }
   
   if(ctr==1) gpr.setPriorMeanToEmpiricalCoefficients();
   if(ctr==2) {
      RealArray1D nu(N+1); for(int k=0;k<=N;k++) nu[k]=r;
      gpr.setPriorMean(nu);
   }  
   gpr.expansionData(f,N);
   delete bFcns;
}


void
GPR::
orthoTest(int q, int m)
{
    // write the matrix B=(P_i(s_j))_{0<=i<=q,0<=j<=m}
	// s_j\in[-1,+1] evenly spaced.
	RealMatrix B(q+1,m+1,0,0);
	for(int j=0;j<=m;j++){

	   Real x=-1.0+j*2.0/m;
	   RealArray1D P=basisFunctionValues(x,q);
	   for(int i=0;i<=q;i++) B(i,j)=P[i];
	}
	// Scale B by 1/sqrt(n), this scales BB' by 1/n and then
	// orthonormality with Monte Carlo integration is equivalent to BB'=I.
	Real f=1.0/sqrt(m+1); B*=f;

	cout << "Orthonormality: matrix of L^2 inner products (psi_i,psi_j):"
        << endl << endl
	     << B.aat();
}


void
GPR::
basisFunctionData(int q, int m)
{
	ofstream fout("BasisFunctionData.txt");

	for(int j=0;j<=m;j++){
      
      Real t=-1.0+j*2.0/m;
	   fout << t << "  ";
      RealArray1D psi_t=basisFunctionValues(t,q);
	   for(int k=0;k<=q;k++) fout << psi_t[k] << "  ";
	   fout << endl;
   }
   fout.close();
}




void
GPR::
printBasisFunctions()
{
   cout << "Printing basis function values in gnuplot data format."
        << "Which basis?" << endl
        << "Legendre polynomials.........[0]" 
        << "Fourier basis................[1]" << endl
        << "Basis=";
   int b; cin>>b;

   BasisFunctions* bFncs;
   switch(b){

      case 0  : bFncs=new LegendreBasis(); break;
      default : bFncs=new FourierBasis();
   }
   int N=2;
   RealArray1D t(2);
   RealArray1D w(2);
	GPR gpr(N,t,w,bFncs);

   gpr.basisFunctionData(7,800);
}



