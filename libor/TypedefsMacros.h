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

#ifndef martingale_typedefsmacros_h
#define martingale_typedefsmacros_h

#include <iostream>                         // once and for all

#define MTGL_BEGIN_NAMESPACE(name) namespace name {
#define MTGL_END_NAMESPACE(name) }


#define SMALL       30               // matrix optimizations below this size.
#define SOBOL_MAX_DIM 300            // maximal dimension of the Sobol generator
#define LMM_MAX_DIM 300              // maximal dimension of a Libor process for
                                     // which a lattice can be built.
#define BASKET_MAX_DIM 300           // maximum number of assets in a basket for
                                     // which a lattice can be built.
#define LATTICE_MAX_STEPS 312        // maximum number of time steps in a two factor 
                                     // lattice (assumes 1GB of main memory)

// the basic scalar type
typedef double Real;




#endif

