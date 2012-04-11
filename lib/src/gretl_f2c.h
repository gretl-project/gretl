#ifndef GRETL_F2C_H
#define GRETL_F2C_H

/* minimal support for f2c'd source files: we assume that
   most fortran-isms will be purged from such files
*/

#include <stdint.h>

typedef int32_t integer;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef int logical;

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
#else
typedef int flag;
typedef int ftnlen;
#endif

#ifndef min
# define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
# define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif

#endif /* GRETL_F2C_H */
