/*
 * input I is a folded spectrum of length N; output O and
 * synthesis window W are of length Nw--overlap-add windowed,
 * unrotated, unfolded input data into output O
 */

#include "fftease.h"

void overlapadd( t_float *I, int N, t_float *W, t_float *O, int Nw, int n )

{
  int i ;
  while ( n < 0 )
    n += N ;
  n %= N ;
  for ( i = 0 ; i < Nw ; i++ ) {
    O[i] += I[n]*W[i] ;
    if ( ++n == N )
      n = 0 ;
  }
}
