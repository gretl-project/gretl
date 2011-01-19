# AC_C_SSE2
# -----------
# Run a double check for sse2: first check that it's supported by the
# compiler, and if so then also check that it's supported by the CPU.
# Some systems may have a cc that's "more advanced" than the CPU, and
# and in that case one could end up with a non-functional binary.
# -----------
AC_DEFUN([AC_C_SSE2],
[
  AC_MSG_CHECKING([whether to use SSE2])
  AC_ARG_ENABLE(sse2,
    [AS_HELP_STRING([--enable-sse2], [use sse2 if available [default=auto]])],
    [enable_sse2=$enableval]
  )
  sse2_result=no
  if test ! "$enable_sse2" = no; then
    if test "x$SSE2_CFLAGS" = "x" ; then
      if test "x$SUNCC" = "xyes"; then
        # SSE2 is enabled by default in the Sun Studio 64-bit environment
        if test "$AMD64_ABI" = "no" ; then
          SSE2_CFLAGS="-xarch=sse2"
        fi
      else
        SSE2_CFLAGS="-msse2"
      fi
    fi
  
    have_sse2_intrinsics=no
    save_CFLAGS=$CFLAGS
    CFLAGS="$SSE2_CFLAGS $CFLAGS"  
  
    AC_COMPILE_IFELSE([
#if defined(__GNUC__) && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
#   if !defined(__amd64__) && !defined(__x86_64__)
#      error "Need GCC >= 4.2 for SSE2 intrinsics on x86"
#   endif
#endif
#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
int main () {
   __m128i a = _mm_set1_epi32 (0), b = _mm_set1_epi32 (0), c;
   c = _mm_xor_si128 (a, b);
   return 0;
}      
    ], have_sse2_intrinsics=yes,)
    if test "$have_sse2_intrinsics" = "yes" ; then
      AC_RUN_IFELSE([
#include <stdio.h>
#define cpuid(func,ax,bx,cx,dx)\
	__asm__ __volatile__ ("cpuid":\
	"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func));
int main (void)
{
    int SSE2 = 0x01 << 26;
    int a, b, c, d;
    cpuid(0x01, a, b, c, d);
    if (d & SSE2) {
	return 0;
    } else {
	return 1;
    }
}
      ], sse2_result=yes,)
    fi
 
    if test "$sse2_result" = "yes" ; then      
       AC_DEFINE(USE_SSE2)
    else
       CFLAGS="$save_CFLAGS"  
    fi      
  fi
  AC_MSG_RESULT([$sse2_result])
])

