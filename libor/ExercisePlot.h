/* Created by Anjuta version 0.1.9 */
/*	This file will not be overwritten */

#include <plotter.h>
#include <fstream>
#include "BermudanOption.h"
#include "LmmLattice.h"


MTGL_BEGIN_NAMESPACE(Martingale)


struct Square { int e; int n; Square() : e(0), n(0) { } };
	
/** Allocates a lattice for driftless LMM with constant volatility surface
 *  and plots the exercise decision for a Bermudan swaption on \f$[T_p,T_n\f$
 *  at all nodes at time \f$T_t\f$ as a function of the 
 *  statistics \f$x=L_t(T_t)\f$ and \f$y=S_{t+1,n}(T_t)\f$ used in the Jaeckel exercise
 *  boundary parametrization (P. Jaeckel, <i>Monte Carlo Methods in Finance</i>, 
 *  pp. 171-181). We traverse all these nodes in the lattice
 *  and compute the pair (x,y) at each node. With k denoting the strike rate the range
 *  \f[0\leq x,y\leq 5k\f]
 *  is divided into 100^2 squares by dividing each axis into 100 intervals of equal length.
 *  Each square Q keeps count how often we have exercise (e) or no exercise (n) for a pair
 *  (x,y) falling into Q and is colored according to the exercise probability u=e/(e+n) as 
 *  follows:
 *  <p>
 *  \f$u\geq 0.99\f$: black (always exercised).<br>
 *  \f$0.8\leq u < 0.99\f$: brown.<br> 
 *  \f$0.6\leq u < 0.8\f$: green.<br> 
 *  \f$0.4\leq u < 0.6\f$: red (maximum ambiguity).<br> 
 *  \f$0.2\leq u < 0.4\f$: blue.<br> 
 *  \f$0.01\leq u < 0.2\f$: cyan.<br> 
 *  \f$u < 0.01\f$: yellow (never exercised).
 * </p>
 * <p>To make this meaningful we must have enough nodes to cover the significant squares
 * with sufficently many pairs (x,y). We use the default lattice which allocates close
 * to the maximum number of nodes that will fit into memory. We should then let t be 
 * 2n/3 or bigger.</p>
 *
 * <p>This needs the <strong>Gnu Plotutils</strong> library (&lt;plotter.h&gt;).
 */
void plotBermudanExercise()
{	
	cout << "Plotting exercise decision for Bermudan swaption on [T_p,T_n] "
	     << "at time T_t";
	cout << "\nEnter n = ";
	int n; cin >> n;
	cout << "\nEnter p = ";
	int p; cin>>p;
	cout << "\nEnter number of grid points on each axis: N = ";
	int N; cin>>N;
	
	// one training path (don't need these)
	BermudanSwaption* bswpn=BermudanSwaption::sample(p,n,1);
	// step size in x and y direction (square dimension)
	Real kappa = bswpn->getStrike();
	// coordinate system 0<=x,y<=a, 100 steps
	Real a=5*kappa, dx=a/N, dy=dx;
	

    // the lattice 	
	LmmLattice* lattice = bswpn->getDefaultLattice();
	typedef LmmLattice::NodeType Node;

	// price the swaption in the lattice, this sets the fields
	// pi=max{h_t,V_{t+1}} in each node
	bswpn->latticeForwardPrice(lattice);
	
	Array2D<Square*> squares(N,N,0,0);
	for(int i=0;i<N;i++)
	for(int j=0;j<N;j++) squares(i,j)=new Square();
	
    // loop over t
	for(int t=2*n/3-2;t<5*n/6;++t){
		
	    // iterate through the nodes and count exercise -- no exercise
	    // count the nodes
	    int nds=0,                            // number of nodes
	        steps = lattice->getSteps(),      // time steps per Libor accrual interval
	        s=t*steps;                        // time step in lattice at time T_t
	

        ofstream log_out("log.txt");
	
     	vector<Node*>* nodes = lattice->getNodeList(s);
	    vector<Node*>::iterator theNode=nodes->begin();
	    while(theNode!=nodes->end()) {
	
		    Node* node=*theNode;
		    Real h = bswpn->forwardPayoff(node,lattice,s);
		    Real V = node->getValue();
// sanity check
if(h>V) cerr << "h > V detected." << endl;
		    bool exercise = ((h>=V)&&(h>0.0));
		    // the statistics
		    Real x = lattice->L(t,node,s);
		    Real y = lattice->swapRate(t+1,n,node,s);

		    // the indices of the square Q containing the point (x,y)
		    // cutoff values out of range
		    //int i=(int)(x/dx), j=(int)(y/dy);
		    int i=floor(x/dx), j=floor(y/dy);
		    i=min(i,N-1); j=min(j,N-1);
// absurd
if((h>0.0)&&(y<kappa)) 
log_out << "Absurd: h = " << h << ", y-kappa = " << y-kappa 
        << ",   (x,y)=("<<x<<','<<y<<")." << endl;
	
	    	// count exercise / no exercise
		    if(exercise) (squares(i,j)->e)++;
		    else (squares(i,j)->n)++;
		
	    	theNode++; nds++;
			
    	} // end while(theNode)
        cout << "\n\n\nNodes: " << nds;	
	
	    // plot the squares
	    // set Plotter parameters
	    char* page="letter";
        Plotter::parampl ("PAGESIZE",page);
	
		// string conversion
	    Int_t tt(t);   
	    string outfile="Exercise"+tt.toString()+".ps";				
        ofstream fout(outfile.c_str());
     
    	//XPlotter plotter(cin, fout, cerr);  // declare X11-Plotter
        PSPlotter plotter(cin, fout, cerr);   // declare PostScript-Plotter
        if (plotter.openpl () < 0)            // open Plotter
        {
             cerr << "Couldn't open Plotter\n";
             exit(1);
         }

     	// x,y range: lower left (0,0) to upper right corner (a,a)
        plotter.fspace (0.0,0.0,a,a); 
	    plotter.erase ();                 // erase Plotter's graphics display
	    // draw axes y=kappa, x=kappa
	    plotter.flinewidth (dx/12);       // line thickness in user coordinates
	    plotter.fline(0.0,kappa,a,kappa);
	    plotter.fline(kappa,0.0,kappa,a);
	

	    // iterate through the squares
	    plotter.flinewidth (0.0); 
	    for(int i=0;i<N;i++)
	    for(int j=0;j<N;j++){
		
		    Square* theSquare=squares(i,j);
		    int e=theSquare->e,                // number of times exercised
		        n=theSquare->n;                // number of times not exercised                 
		
		    // if we have a sample in square(i,j)=Q
		    if(e+n>0){
			
		        Real ep = (1.0*e)/(e+n);             // exercise probability in Q
		
		       // objects (rectangles) will be filled 
		       plotter.filltype(1);              
		       if(ep>=0.99) plotter.fillcolorname("black");
			   if((0.8<=ep)&&(ep<0.99)) plotter.fillcolorname("brown");
		       if((0.6<=ep)&&(ep<0.8)) plotter.fillcolorname("green");
		       if((0.4<=ep)&&(ep<0.6)) plotter.fillcolorname("red");
		       if((0.2<=ep)&&(ep<0.4)) plotter.fillcolorname("blue");
		       if((0.01<=ep)&&(ep<0.2)) plotter.fillcolorname("cyan");
			   if((0.0<=ep)&&(ep<0.01)) plotter.fillcolorname("yellow");	
			  	
               plotter.fbox(i*dx,j*dy,(i+1)*dx,(j+1)*dy); 
		   } // end if	
	   } // end for i,j
	   
	   // close Plotter
       if (plotter.closepl () < 0)         
       {
            cerr << "Couldn't close Plotter\n";
       }
	   // clear squares
	   for(int i=0;i<N;i++)
       for(int j=0;j<N;j++){
		   
		   squares(i,j)->e=0;
		   squares(i,j)->n=0;
	   }
		   
    } // end for t
            
    for(int i=0;i<N;i++)
    for(int j=0;j<N;j++) delete squares(i,j);
 
} //end plotBermudanExercise()


	
	
	

MTGL_END_NAMESPACE(Martingale)