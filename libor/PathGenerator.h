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

#ifndef martingale_pathgenerator_h    
#define martingale_pathgenerator_h

#include "TypedefsMacros.h"


MTGL_BEGIN_NAMESPACE(Martingale)



/** <p> Interface to all entities generating paths.
 *  Forward a request for a new path to the concrete subclasses.
 *  Implicit assumption: the new path is stored and accessible as <i>the current
 *  path</i> of the PathGenerator.
 */
class PathGenerator {
	
public:
	
	virtual void newPath() = 0;
		
}; // end PathGenerator




MTGL_END_NAMESPACE(Martingale)

#endif
 
