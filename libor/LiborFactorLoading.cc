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




#include "LiborFactorLoading.h"
#include "Utils.h"
#include <cmath>
using namespace Martingale;
using namespace std;




/*******************************************************************************
 *
 *                     Class: LiborFactorLoading
 * 
 ******************************************************************************/


 
    LiborFactorLoading::
    LiborFactorLoading(int dim, Real* L0, Real* deltas) :
	n(dim), 
	delta(new Real[n]), 
	Tc(new Real[n+1]),
	l(new Real[n]),
	x(new Real[n]),
	alpha(new Real[n]),
	beta(new Real[n])
    {   
	    // tenor structure, initial XLibors
		Tc[0]=0;
		for(int i=0;i<n;i++){ 
			
			delta[i]=deltas[i];
			Tc[i+1]=Tc[i]+delta[i];
			l[i]=L0[i];
			x[i]=delta[i]*l[i];
		}
		
		// drift linearization constants alpha_k, beta_k
        // for the X1-dynamics (includes X1_0, no index downshift j->j-1)
        for(int j=0;j<n;j++){
            
            Real r,s,a,b;
			if(j%2==0) { r=0.5; s=0.8; }    // alternate over/under estimation
			else       { r=0.1; s=0.1; }
		
			a=log(x[j])-r*log(3);
            b=log(x[j])+s*log(3);
            
            alpha[j]=(b*f(a)-a*f(b))/(b-a);
            beta[j] =(f(b)-f(a))/(b-a);
        }
	} // end constructor
	
	

			

// LOG-COVARIATION MATRICES AND DRIFT LINEARIZATION MATRICES

   

   UTRMatrix<Real>& LiborFactorLoading::
   getRho() const
   {
	   UTRMatrix<Real>& _rho=*(new UTRMatrix<Real>(n-1,1));
	   for(int i=1;i<n;i++)
	   for(int j=i;j<n;j++) _rho(i,j)=rho(i,j);
		   
	   return _rho;
   }
	
	
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int p,int q, Real t, Real T) const
   {
       int size=q-p;       // matrix size
	   UTRMatrix<Real>& logCVM=*(new UTRMatrix<Real>(size,p));

       for(int i=p;i<q;i++)
       for(int j=i;j<q;j++)
       logCVM(i,j)=integral_sgi_sgj_rhoij(i,j,t,T);

       return logCVM;
   }// end logCovariationMatrix
   

   
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrix(int t) const
   {
       Real Tt=Tc[t],         // T_t
	        Tt1=Tc[t+1];      // T_{t+1}
	   return logLiborCovariationMatrix(t+1,n,Tt,Tt1);
   }// end logCovariationMatrix
   

   
   UTRMatrix<Real>& LiborFactorLoading::
   logLiborCovariationMatrixRoot(int t) const
   {
	   return logLiborCovariationMatrix(t).utrRoot();
   }// end logCovariationMatrix


   Matrix<Real>& LiborFactorLoading::
   reducedRankLogLiborCovariationMatrixRoot(int t, int r) const
   {
	   return logLiborCovariationMatrix(t).rankReducedRoot(r);
   }// end logCovariationMatrix


   
   
   
// FUNCTIONS NEEDED FOR THE LOGNORMAL APPROXIMATION
   
   UTRMatrix<Real>& 
   LiborFactorLoading::A(Real t, int p) const
   {
       int size=n-p;       // matrix size
	   UTRMatrix<Real>& a=*(new UTRMatrix<Real>(size,p));

       for(int i=p;i<n;i++){
		   
		   a(i,i)=0;
           for(int j=i+1;j<n;j++)
           a(i,j)=beta[j]*sigma(i,t)*sigma(j,t)*rho(i,j);
	   }
       return a;
   }// end A
   
   
 
   UTRMatrix<Real>& 
   LiborFactorLoading::B(Real t, int p) const
   {
       int size=n-p;       // matrix size
	   UTRMatrix<Real>& a=*(new UTRMatrix<Real>(size,p));

       for(int i=p;i<n;i++){
		   
		   a(i,i)=0;
           for(int j=i+1;j<n;j++) 
		   a(i,j)=beta[j]*integral_sgi_sgj_rhoij(i,j,0,t);
	   }
       return a;
   }// end B
   


   UTRMatrix<Real>& 
   LiborFactorLoading::eB(Real t, int p) const
   {
	   return B(t,p).exp();
   }
   
 
   UTRMatrix<Real>& 
   LiborFactorLoading::eBInverse(Real t, int p) const
   {
	   UTRMatrix<Real>& Bt=B(t,p);
	   Bt*=(-1.0);
	   return Bt.exp();
   }
	   

   vector<Real>& 
   LiborFactorLoading::u(Real t, int i)
   {
	   vector<Real>& ujt=*(new vector<Real>(n-i,i));
	   for(int j=i;j<n;j++){
		   
		   Real sgjt=sigma(j,t),
		        s=0.5*sgjt*sgjt;
		   for(int k=j+1;k<n;k++) s+=alpha[k]*sgjt*sigma(k,t)*rho(j,k);
		   ujt[j]=s;
	   }
	   return ujt;
   }
   


   UTRMatrix<Real>& 
   LiborFactorLoading::C(Real t, int p)
   {
	   UTRMatrix<Real>& c=*(new UTRMatrix<Real>(n-p,p));
	   for(int j=p;j<n;j++) 
	   for(int k=j;k<n;k++) c(j,k)=sigma(j,t)*sigma(k,t)*rho(j,k);
		   
	   return c;
   }
   

   UTRMatrix<Real>& 
   LiborFactorLoading::nu(Real t, int p)
   {
	   UTRMatrix<Real> c(n-p,p);
	   for(int j=p;j<n;j++) 
	   for(int k=j;k<n;k++) c(j,k)=sigma(j,t)*sigma(k,t)*rho(j,k);
		   
	   return c.utrRoot();
   }
   
   

// TEST PROGRAMS 
		 

void LiborFactorLoading::selfTest()
{
	Real precision=0.001,       // maximum acceptable relative error in percent
		 epsilon=0.00000000001; // zero denominator reset to epsilon
	
    cout << "\nLIBOR FACTORLOADING SELFTEST:" << endl << toString();
	printFields();
	
	cout << "\nTesting the root L of the matrix C=logLiborCovariationMatrix(t):" << endl;
	for(int t=0;t<n-1;t++){
		
		UTRMatrix<Real>& C=logLiborCovariationMatrix(t);
		UTRMatrix<Real>& L=logLiborCovariationMatrixRoot(t);
		cout << "\nt = " << t << ": ";
		C.testEquals(L.aat(),precision,epsilon,"C=L*L'");
	}
	
	cout << "\n\n\nTesting the root nu of C(t):" << endl;
	for(int t=1;t<n-1;t++){            // not t=0, integrals!
		
		Real T_t=Tc[t];
		UTRMatrix<Real>& c=C(T_t,t);
		UTRMatrix<Real>& n=nu(T_t,t);
		cout << "\nt = " << t << ": ";
		c.testEquals(n.aat(),precision,epsilon,"C=nu*nu'");
	}
	
	
	cout << "\n\n\nTesting the matrix exponentials eB(t), eBInverse(t):" << endl;
	precision=0.001,       // maximum acceptable relative error in percent
    epsilon=0.0001;        // zero denominator reset to epsilon
	for(int t=1;t<n-1;t++){
		
	    Real T_t=Tc[t];
		UTRMatrix<Real>& eBt=eB(T_t,t);
		UTRMatrix<Real>& eBtinv=eBInverse(T_t,t);
		eBt*=eBtinv;     // now must be the identity matrix
		UTRMatrix<Real> I(n-t,t);
		for(int j=t;j<n;j++) I(j,j)=1.0;
		cout << "\nt = " << t << ": ";
		I.testEquals(eBt,precision,epsilon,"I=eBt*eBtinv");
	}

} // end test


void LiborFactorLoading::factorizationTest(int r)
{
	cerr << "\n\nApproximate factorization of all single time step log-Libor"
	     << "\ncovariation matrices C(t) as C(t)= R(t)R(t)' with R(t) of rank " << r
	     << "\n\nRelative errors (trace norm): " << endl;
	
	for(int t=0;t<n-2;t++){
		
		UTRMatrix<Real>& Ct=logLiborCovariationMatrix(t);
		Ct.testFactorization(r);
	}
}




/*******************************************************************************
 *
 *                       Class: CS_FactorLoading
 *
 ******************************************************************************/


// CONSTRUCTOR

    CS_FactorLoading::CS_FactorLoading
    (int dim, Real a, Real d, Real Alpha, Real Beta, Real Rho00, Real* C, Real* l, Real* deltas) :
	LiborFactorLoading(dim,l,deltas),
	A(a), D(d), alpha(Alpha), beta(Beta), rho00(Rho00), c(C),
	corr(dim-1,1)
	{
       // allocate and initialize the correlation base b[]=(b_1,...,b_{n-1})
	   int n=getDimension(); 
       Real b[n-1]; 
       for(int i=0;i<n-1;i++)
       { Real x=((Real)i)/(n-2); b[i]=exp(-f(x)); }
           
       // initialize the correlation matrix corr[i][j]=rho_{ij}
       for(int i=1;i<n;i++)
       for(int j=i;j<n;j++) corr(i,j)=b[i-1]/b[j-1]; 
  
    } // end constructor
	
	
                                                                                                        
       
// VOLATILITIES, CORRELATIONS, COVARIATION INTEGRALS  
  
 
   Real CS_FactorLoading::
   rho(int i, int j) const 
   { if(i<=j) return corr(i,j); return corr(j,i); } 


   Real CS_FactorLoading::
   sigma(int i, Real t) const 
   { 
	   Real T_i=getTenorStructure()[i];
	   return c[i]*g(1-t/T_i); 
   }

   
   Real CS_FactorLoading::
   integral_sgi_sgj_rhoij(int i, int j, Real t, Real T) const
   {
       Real* Tc=getTenorStructure();
	   Real T_i=Tc[i], T_j=Tc[j];
	   return c[i]*c[j]*rho(i,j)*integral_gg(t,T,T_i,T_j); 
   }


   
// SAMPLE FACTOR LOADING

   
    CS_FactorLoading* 
    CS_FactorLoading::sample(int n, Real delta=0.25) 
    {
        Real A=1.5, D=2, alpha=1.8, beta=0.1, rho00=0.4;
        
        Real* deltas=new Real[n];
		Real* c=new Real[n];
		Real* l=new Real[n];
        for(int i=0;i<n;i++){ deltas[i]=delta; c[i]=0.25; l[i]=0.04; }

        return new CS_FactorLoading(n,A,D,alpha,beta,rho00,c,l,deltas);
     
    } // end sample
   

	
// STRING REPRESENTATION
 
    string CS_FactorLoading::toString()  const 
    {
		int n=getDimension();
		Real* l=getInitialTermStructure();
		vector<Real> lvec(n,l);
		ostringstream os;
		os << endl << endl
		   << "\n\nCS_FactorLoading:\nn=" << n << ", A=" << A << ", D=" << D
		   << ", alpha=" << alpha << ", beta=" << beta << ", rho00=" << rho00
		   << endl << endl
		   << "Initial Libors: " << endl << lvec;
        return os.str();
    }
	
	
	void CS_FactorLoading::printFields()
    {
		int n=getDimension();
		Real* deltas=getDeltas();
		Real* Tc=getTenorStructure();
		Real* x=getInitialXLibors();
		vector<Real> delta_vec(n,deltas);
		vector<Real> Tc_vec(n+1,Tc);
		vector<Real> x_vec(n,x);
		vector<Real> c_vec(n,c);
		cout << "\ndeltas:" << endl << delta_vec
		     << "\nReset times T_j:" << endl << Tc_vec
		     << "\nInitial X-Libors:" << endl << x_vec
		     << "\nc_j values:" << endl << c_vec
		     << "\nCorrelation matrix: " << endl << corr;
	}
    

    Real CS_FactorLoading::
	integral_g_squared(Real T) const 
    { 
        if(T==0) return 0;
		Real R=T;
        R+=2*A*(G(-D,T)+D*D);
        R+=A*A*(H(-D/2,T)+D*D*D/4);
        //now R=int_0^T g^2(s)ds, see Tex document LMM.
        return R;
    } //end Integral_g_squared
    

	
    Real CS_FactorLoading::
	integral_gg(Real t, Real T, Real a, Real b) const 
    {
        if(t==T) return 0;
		
		Real R,d,da,db,Da,Db,Dab;
    
        d=A*exp(-1/D); da=d/a; db=d/b;       
        Da=D*a; Db=D*b; Dab=Da*Db/(Da+Db);
  
        R=T-t;
        R+=d*(F(Da,T)-F(Da,t));
        R+=d*(F(Db,T)-F(Db,t));
        R+=d*d*(F(Dab,T)-F(Dab,t));
  
        R-=da*(G(Da,T)-G(Da,t));
        R-=db*(G(Db,T)-G(Db,t));
        R-=d*(da+db)*(G(Dab,T)-G(Dab,t));
  
        R+=da*db*(H(Dab,T)-H(Dab,t));      
        //now R=int_t^Tg(1-u/a)g(1-u/b)du
  
       return R;
    
   } //end Integral_gg
        

               
/*******************************************************************************
 *
 *                     JR_FactorLoading
 *
 *******************************************************************************/        



// CONSTRUCTOR

    JR_FactorLoading::JR_FactorLoading
    (int dim, Real* L0, Real* deltas, Real a0, Real b0, Real c0, Real d0, Real Beta0, Real* k0) : 
    LiborFactorLoading(dim,L0,deltas),
    a(a0), b(b0), c(c0), d(d0), Beta(Beta0), k(k0), 
    corr(dim-1,1)
    {
       // initialize correlations rho_ij=exp(beta*(T_i-T_j)), i<=j.
	   Real* Tc=getTenorStructure();
       for(int i=1;i<dim;i++)
       for(int j=i;j<dim;j++) corr(i,j)=exp(Beta*(Tc[i]-Tc[j]));   
       
    } // end constructor
	

  
// STRING MESSAGE

    string JR_FactorLoading::toString() const
    {
        int n=getDimension();
		Real* l=getInitialTermStructure();
		vector<Real> lvec(n,l);
		ostringstream os;
		os << endl << endl
           << "\nJR_FactorLoading:" << endl
           << "n=" << n
		   <<", a="<< a
		   <<", b="<< b
		   <<", c="<< c
		   <<", d="<< d
		   <<", beta="<< Beta
		   << endl << endl
		   << "Initial Libors: " << endl << lvec;
        return os.str();
    }
	
	
	void JR_FactorLoading::printFields()
    {
		int n=getDimension();
		Real* deltas=getDeltas();
		Real* Tc=getTenorStructure();
		Real* x=getInitialXLibors();
		vector<Real> delta_vec(n,deltas);
		vector<Real> Tc_vec(n+1,Tc);
		vector<Real> x_vec(n,x);
		vector<Real> kvec(n,k);
		cout << "\ndeltas:" << endl << delta_vec
		     << "\nReset times T_j:" << endl << Tc_vec
		     << "\nInitial X-Libors:" << endl << x_vec
		     << "\nk_j values:" << endl << kvec
		     << "\nCorrelation matrix: " << endl << corr;
	}

	
// SAMPLE FACTOR LOADING

    JR_FactorLoading* JR_FactorLoading::
	sample(int n, Real delta=0.25)
    {
        Real a0=-0.05, b0=0.5, c0=1.5, d0=0.15, Beta0=0.1;
        
	    Real* L0=new Real[n];
        Real* k0=new Real[n];
        Real* deltas=new Real[n];       
        for(int i=0;i<n;i++){ k0[i]=2.6; L0[i]=0.04; deltas[i]=delta; }
        
        return new 
        JR_FactorLoading(n,L0,deltas,a0,b0,c0,d0,Beta0,k0);
     
    } // end sample
   

   Real JR_FactorLoading::
   Fij(int i, int j, Real t) const
   {
       Real* Tc=getTenorStructure();
	   Real f,A,B,C,
              ctmti=c*(t-Tc[i]),
              ctmtj=c*(t-Tc[j]),
              timtj=fabs(Tc[i]-Tc[j]),
              q=ctmti+ctmtj,                    // c(2t-ti-tj)
              ac=a*c,
              cd=c*d,
              bcd=b*c*d;
              
       f=k[i]*k[j]*exp(-Beta*timtj)/(c*c*c);
       A=ac*cd*(exp(ctmtj)+exp(ctmti))+c*cd*cd*t;
       B=bcd*(exp(ctmti)*(ctmti-1)+exp(ctmtj)*(ctmtj-1));
       C=exp(q)*(ac*(ac+b*(1-q))+b*b*(0.5*(1-q)+ctmti*ctmtj))/2;
       
       return f*(A-B+C);
   } // end Fij
   
   
   
               
/*******************************************************************************
 *
 *               Constant Volatility FactorLoading
 *
 *******************************************************************************/        


// CONSTRUCTOR
   
   
    // CS - correlations
    ConstVolLiborFactorLoading::
    ConstVolLiborFactorLoading
    (int dim, Real* L0, Real* deltas, Real _alpha, Real _beta, Real _r_oo, Real* _sg) : 
	LiborFactorLoading(dim,L0,deltas),
	alpha(_alpha), beta(_beta), r_oo(_r_oo), sg(_sg),
	corr(dim-1,1)
	{
       // allocate and initialize the correlation base b[]=(b_1,...,b_{n-1})
	   int n=getDimension(); 
       Real b[n-1]; 
       for(int i=0;i<n-1;i++)
       { Real x=((Real)i)/(n-2); b[i]=exp(-f(x)); }
           
       // initialize the correlation matrix corr[i][j]=rho_{ij}
       for(int i=1;i<n;i++)
       for(int j=i;j<n;j++) corr(i,j)=b[i-1]/b[j-1]; 
  
    } // end constructor

	
	// JR - correlations
    ConstVolLiborFactorLoading::
    ConstVolLiborFactorLoading
    (int dim, Real* L0, Real* deltas, Real _beta, Real* _sg) : 
    LiborFactorLoading(dim,L0,deltas),
    alpha(0.0), beta(_beta), r_oo(1.0), sg(_sg), 
	corr(dim-1,1)
    {
       // initialize correlations rho_ij=exp(beta*(T_i-T_j)), i<=j.
	   Real* Tc=getTenorStructure();
       for(int i=1;i<dim;i++)
       for(int j=i;j<dim;j++) corr(i,j)=exp(beta*(Tc[i]-Tc[j]));   
       
    } // end constructor
	

  
// STRING MESSAGE

    string ConstVolLiborFactorLoading::
	toString() const
    {
        int n=getDimension();
		Real* l=getInitialTermStructure();
		vector<Real> lvec(n,l);
		vector<Real> vols(n,sg);      // volatilities
		ostringstream os;
		os << endl << endl
           << "\nConstVolLiborFactorLoading:" 
           << "\nVolatilities: " << endl << vols
		   << "\nInitial Libors: " << endl << lvec;
        return os.str();
    }
	
	
	void ConstVolLiborFactorLoading::
	printFields()
    {
		int n=getDimension();
		Real* deltas=getDeltas();
		Real* Tc=getTenorStructure();
		Real* x=getInitialXLibors();
		vector<Real> delta_vec(n,deltas);
		vector<Real> Tc_vec(n+1,Tc);
		vector<Real> x_vec(n,x);
		vector<Real> vols(n,sg);
		cout << "\ndeltas:" << endl << delta_vec
		     << "\nReset times T_j:" << endl << Tc_vec
		     << "\nInitial X-Libors:" << endl << x_vec
		     << "\nVolatilities:" << endl << vols
		     << "\nCorrelation matrix: " << endl << corr;
	}

	
// SAMPLE FACTOR LOADING

    ConstVolLiborFactorLoading* ConstVolLiborFactorLoading::
	sample(int n, Real delta=0.25, int corrs=CS)
    {
	    Real* _L=new Real[n];
        Real* _sg=new Real[n];
        Real* _deltas=new Real[n];       
        for(int i=0;i<n;i++){ _sg[i]=0.4; _L[i]=0.04; _deltas[i]=delta; }		
		
		Real _alpha=1.8, _beta=0.1, _r_oo=0.4;
			
		if(corrs==CS)
		return new ConstVolLiborFactorLoading(n,_L,_deltas,_alpha,_beta,_r_oo,_sg);
        // else corrs == JR
        return new ConstVolLiborFactorLoading(n,_L,_deltas,_beta,_sg);
     
    } // end sample
        
