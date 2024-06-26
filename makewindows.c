#include "fftease.h"

void lpp_makewindows( t_float *H, t_float *A, t_float *S, int Nw, int N, int I )

{
  int i ;
  t_float sum ;

  for ( i = 0 ; i < Nw ; i++ )
    H[i] = A[i] = S[i] = 0.54 - 0.46*cos( TWOPI*i/(Nw - 1) ) ;

  if ( Nw > N ) {
    t_float x ;

    x = -(Nw - 1)/2. ;
    for ( i = 0 ; i < Nw ; i++, x += 1. )
      if ( x != 0. ) {
        A[i] *= N*sin( PI*x/N )/(PI*x) ;
        if ( I )
          S[i] *= I*sin( PI*x/I )/(PI*x) ;
      }
  }

  for ( sum = i = 0 ; i < Nw ; i++ )
    sum += A[i] ;

  for ( i = 0 ; i < Nw ; i++ ) {
    t_float afac = 2./sum ;
    t_float sfac = Nw > N ? 1./afac : afac ;
    A[i] *= afac ;
    S[i] *= sfac ;
  }

  if ( Nw <= N && I ) {
    for ( sum = i = 0 ; i < Nw ; i += I )
      sum += S[i]*S[i] ;
    for ( sum = 1./sum, i = 0 ; i < Nw ; i++ )
      S[i] *= sum ;
  }
}

void lpp_makehamming( t_float *H, t_float *A, t_float *S, int Nw, int N, int I, int odd )

{
  int i;
  t_float sum ;



  if (odd) {
    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = sqrt(0.54 - 0.46*cos( TWOPI*i/(Nw - 1) ));
  }

  else {

    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = 0.54 - 0.46*cos( TWOPI*i/(Nw - 1) );

  }

  if ( Nw > N ) {
    t_float x ;

    x = -(Nw - 1)/2. ;
    for ( i = 0 ; i < Nw ; i++, x += 1. )
      if ( x != 0. ) {
        A[i] *= N*sin( PI*x/N )/(PI*x) ;
        if ( I )
          S[i] *= I*sin( PI*x/I )/(PI*x) ;
      }
  }
  for ( sum = i = 0 ; i < Nw ; i++ )
    sum += A[i] ;

  for ( i = 0 ; i < Nw ; i++ ) {
    t_float afac = 2./sum ;
    t_float sfac = Nw > N ? 1./afac : afac ;
    A[i] *= afac ;
    S[i] *= sfac ;
  }

  if ( Nw <= N && I ) {
    for ( sum = i = 0 ; i < Nw ; i += I )
      sum += S[i]*S[i] ;
    for ( sum = 1./sum, i = 0 ; i < Nw ; i++ )
      S[i] *= sum ;
  }
}


void lpp_makehanning( t_float *H, t_float *A, t_float *S, int Nw, int N, int I, int odd )
{
  int i;
  t_float sum ;


  if (odd) {
    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = sqrt(0.5 * (1. + cos(PI + TWOPI * i / (Nw - 1))));
  }

  else {

    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = 0.5 * (1. + cos(PI + TWOPI * i / (Nw - 1)));

  }

  if ( Nw > N ) {
    t_float x ;

    x = -(Nw - 1)/2. ;
    for ( i = 0 ; i < Nw ; i++, x += 1. )
      if ( x != 0. ) {
        A[i] *= N*sin( PI*x/N )/(PI*x) ;
        if ( I )
          S[i] *= I*sin( PI*x/I )/(PI*x) ;
      }
  }
  for ( sum = i = 0 ; i < Nw ; i++ )
    sum += A[i] ;

  for ( i = 0 ; i < Nw ; i++ ) {
    t_float afac = 2./sum ;
    t_float sfac = Nw > N ? 1./afac : afac ;
    A[i] *= afac ;
    S[i] *= sfac ;
  }

  if ( Nw <= N && I ) {
    for ( sum = i = 0 ; i < Nw ; i += I )
      sum += S[i]*S[i] ;
    for ( sum = 1./sum, i = 0 ; i < Nw ; i++ )
      S[i] *= sum ;
  }
}
