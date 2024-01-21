/*
 *  utility.c
 *  twoyukawa
 *
 *  Created by Marcus Hennig on 5/13/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#if 0
#include "TY_utility.h"
#endif
double chop( double x )
{
	if ( fabs(x) < 1E-6 )
		return 0;
	else 
		return x;
}