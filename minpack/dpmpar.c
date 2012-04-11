#include <float.h>

double dpmpar_ (int k)
{
    if (k == 1) {
	return DBL_EPSILON;
    } else if (k == 2) {
	return DBL_MIN;
    } else if (k == 3) {
	return DBL_MAX;
    } else {
	return 0;
    }
}


