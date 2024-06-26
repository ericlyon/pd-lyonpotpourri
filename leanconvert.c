#include "fftease.h"

void leanconvert( t_float *S, t_float *C, int N2 )
{
  int    real, imag, amp, phase;
  t_float  a, b;
  int    i;
  double hypot(), atan2();

  for ( i = 0; i <= N2; i++ ) {
    imag = phase = ( real = amp = i<<1 ) + 1;
    a = ( i == N2 ? S[1] : S[real] );
    b = ( i == 0 || i == N2 ? 0. : S[imag] );
    C[amp] = hypot( a, b );
    C[phase] = -atan2( b, a );
  }
}
