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

#include <cassert>
#include <iostream>                         // once and for all
using std::cout;
using std::cerr;
using std::cin;
using std::endl;


#define MTGL_BEGIN_NAMESPACE(name) namespace name {
#define MTGL_END_NAMESPACE(name) }


#define TWO_PI 6.28318531
#define PI 3.14159266
#define SMALL 30

// the basic scalar type
typedef double Real;
// Real function of one real variable
typedef Real (*RealFunction)(Real);




#endif
