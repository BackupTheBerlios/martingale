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
spyqqqdia@yahoo.com

*/


#include "PathGenerator.h"
#include "LiborMarketModel.h"
#include "Trigger.h"


MTGL_BEGIN_NAMESPACE(Martingale)

// not in header to avoid includes in the header
void 
LiborPathsToTriggerTime:: 
newPath()
{ 
    s=0;
	while(!(trigger->isTriggered(0,s))){ LMM->timeStep(s,i); ++s; }
}


MTGL_END_NAMESPACE(Martingale)

 
