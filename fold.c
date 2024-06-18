#include "fftease.h"

/*
 * multiply current input I by window W (both of length Nw);
 * using modulus arithmetic, fold and rotate windowed input
 * into output array O of (FFT) length N according to current
 * input time n
 */
void lpp_fold( t_float *I, t_float *W, int Nw, t_float *O, int N, int n )
{
  int i;

  for ( i = 0; i < N; i++ )
    O[i] = 0.;

  while ( n < 0 )
    n += N;
  n %= N;
  for ( i = 0; i < Nw; i++ ) {
    O[n] += I[i]*W[i];
    if ( ++n == N )
      n = 0;
  }
}
