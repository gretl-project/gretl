# AC_C_AVX
# -----------
# Run a double check for avx: first check that it's supported by the
# compiler, and if so then also check that it's supported by the CPU.
# Some systems may have a cc that's "more advanced" than the CPU, and
# and in that case one could end up with a non-functional binary.
#
# If we are cross-compiling, however, we assume that if the (cross-)
# compiler supports avx then so will the target CPU. If that 
# assumption is invalid, it will be necessary to pass the option
# --disable-avx to the configure script explicitly.
# -----------
AC_DEFUN([AC_C_AVX],
[
  AC_MSG_CHECKING([whether to use AVX])
  AC_ARG_ENABLE(avx,
    [AS_HELP_STRING([--enable-avx], [use avx if available [default=auto]])],
    [enable_avx=$enableval]
  )
  avx_result=no
  save_CFLAGS=$CFLAGS
  AVX_CFLAGS=
  if test "$enable_avx" != "no" ; then
    if test "x$AVX_CFLAGS" = "x" ; then
      if test "x$SUNCC" = "xyes" && test "x$AMD64_ABI" = "xno" ; then
        AVX_CFLAGS="-xarch=avx"
      else
        AVX_CFLAGS="-mavx -Winline"
      fi
    fi
  
    have_avx_intrinsics=no
    CFLAGS="$AVX_CFLAGS $CFLAGS"  
  
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 5))
#   if !defined(__amd64__) && !defined(__x86_64__)
#      error "Need GCC >= 4.5 for AVX intrinsics on x86"
#   endif
#endif
#include <immintrin.h>
int main () {
    __m256i a = _mm256_set1_epi32 (0), b = _mm256_set1_epi32 (0), c;
	c = _mm256_permute2f128_si256 (a, b, 0);
    return 0;
}      
    ])], have_avx_intrinsics=yes,AVX_CFLAGS="")
    if test "$have_avx_intrinsics" = "yes" ; then
      AC_RUN_IFELSE([AC_LANG_SOURCE([
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#define GET_FEATURE_FLAGS 1

int main (void)
{
    int32_t ecx, edx;
    __asm__("cpuid"
	    : "=c" (ecx),
	      "=d" (edx)
	    : "a" (GET_FEATURE_FLAGS)
	    );
    if ((ecx >> 28) & 1) {
	return 0;
    } else {
        return 1;
    }
}
      ])], avx_result=yes, avx_result=no, avx_result=yes)
    fi
 
    if test "$avx_result" = "yes" ; then      
       AC_DEFINE(USE_AVX)
    else
       AVX_CFLAGS=""
    fi 
  fi
  CFLAGS="$save_CFLAGS"
  AC_MSG_RESULT([$avx_result])
  AC_SUBST([AVX_CFLAGS])
])

