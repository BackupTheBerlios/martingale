/* Created by Anjuta version 0.1.9 */
/*	This file will not be overwritten */


#ifndef martingale_exerciseplot_h    
#define martingale_exerciseplot_h


#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)


struct Square { int e; int n; Square() : e(0), n(0) { } };
	
/** Allocates a lattice for driftless LMM with constant volatility surface
 *  and plots the exercise decision for a Bermudan swaption on \f$[T_p,T_n]\f$
 *  at all nodes at several times \f$T_t\f$ as a function of the 
 *  statistics 
 *  \f[x=L_t(T_t)\quad\hbox{and}\quad y=S_{t+1,n}(T_t)\f] 
 *  used in the Jaeckel exercise boundary parametrization 
 *  (P. Jaeckel, <i>Monte Carlo Methods in Finance</i>, pp. 171-181). 
 *  We traverse all these nodes in the lattice and compute the pair (x,y) at each node. 
 *  With k denoting the strike rate the range
 *  \f[0\leq x,y\leq 5k\f]
 *  is divided into N^2 squares by dividing each axis into N (user supplied) intervals of 
 *  equal length.
 *  Each square Q keeps count how often we have exercise (e) or no exercise (n) for a pair
 *  (x,y) falling into Q and is colored according to the exercise probability u=e/(e+n).
 *  The range of probabilities u is subdivided as follows:
 *  \f[0<0.01<0.2<0.4<0.6<0.8<0.99<1.\f]
 *  with respective colors
 *  <p>
 *  <table border="1" cellspacing="3" cellpadding="3">
 *  <tr>
 *  <td><font color="yellow">yellow:</font> </td><td>u in [0.0,0.01] </td><td>(never exercised). </td></tr>
 *  <tr>
 *  <td><font color="cyan">cyan:</font> </td><td>u in (0.01,0.2] </td><td></td></tr>
 *  <tr>
 *  <td><font color="blue">blue:</font> </td><td>u in (0.2,0.4] </td><td></td></tr>
 *  <tr>
 *  <td><font color="red">red:</font> </td><td>u in (0.4,0.6] </td><td>(maximum ambiguity). </td></tr>
 *  <tr>
 *  <td><font color="green">green:</font> </td><td>u in (0.6,0.8] </td><td></td></tr>
 *  <tr>
 *  <td><font color="brown">brown:</font> </td><td>u in (0.8,0.99] </td><td></td></tr>
 *  <tr>
 *  <td><font color="black">black:</font> </td><td>u in (0.99,1.0] </td><td>(always exercised). </td></tr>
 *  </table>
 *  </p>
 * <p>To make this meaningful we must have enough nodes to cover the significant squares
 * with sufficently many pairs (x,y). We use the default lattice and must see to it that
 * it allocates close to the maximum number of nodes that will fit into memory. 
 * A trivial change in the code of <code>LiborDerivative::getDefaultLattice()</code>
 * has this effect: count the variable <code>steps</code> down from 100 instead of 6.
 * Suggested values for N: 75-85.
 *
 * <p>This needs the <strong>Gnu Plotutils</strong> library (&lt;plotter.h&gt;).
 * Here is such a <a href="../Exercise44.ps">graph</a>. The file name here indicates that t=44.
 * The example was run with N=75, p=14, n=60. 
 *
 * <p>The vertical and horizontal axes shown are \f$x=k\f$ and \f$y=k\f$.
 * Note the occasional exercise when the first Libor is below the strike (black to the left
 * of x=k) which exercises into an immediate loss on the first swap leg if the
 * subsequent swaprate y is large enough to warrant this.
 * This phenmenon is more pronounced with other parameter values (try p=12, n=40).
 * The exercise boundary shifts quite a bit with the parameters p,n,t.
 */
void plotBermudanExercise();

	
MTGL_END_NAMESPACE(Martingale)

#endif