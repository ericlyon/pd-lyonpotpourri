#include "fftease.h"



void unconvert( t_float *C, t_float *S, int N2, t_float *lastphase, t_float fundamental, t_float factor )

{
  int i, real, imag, amp, freq;
  t_float mag, phase;
  double sin(), cos();

  for ( i = 0; i <= N2; i++ ) {
    imag = freq = ( real = amp = i<<1 ) + 1;

    if ( i == N2 )
      real = 1;

    mag = C[amp];
    lastphase[i] += C[freq] - i*fundamental;
    phase = lastphase[i]*factor;
    S[real] = mag*cos( phase );

    if ( i != N2 )
      S[imag] = -mag*sin( phase );

  }
}
