#ifndef GLOBALS_H_
#define GLOBALS_H_

// IO includes
#include <iostream>
#include <fstream>
// C++ STL includes
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>

// C functions
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <string>
// useful Utilities
#include "Util.h"

// Type definitions
typedef double Double;
typedef long double Double1;
typedef enum {
	BAYES, MARKOV
} NETWORK_TYPE;
typedef enum {
	POSITIVE, DET
} NETWORK_MODE;

#define INVALID_VALUE -1
static bool is_equal(Double d1, Double d2)
{
	if (d1==d2)
		return true;
	if (d1==0.0 || d2==0.0)
		return false;
    if(fabs(log(d1)-log(d2))<0.00000001)
       return true;
    return false;
}


#endif
