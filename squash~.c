#include "MSPd.h"
#include "fftease.h"

#define OBJECT_NAME "squash~"

/* Pd version of squash~ */

static t_class *squash_class;


typedef struct _squash
{

  t_object x_obj;
  t_float x_f;
  int D;
  int N;
  int Nw;
  int N2;
  int incnt;
  int outcnt;
  float *Wanal;
  float *Wsyn;
  float *Hwin;
  float *buffer;
  float *input;
  float *output;
  float thresh;
  float ratio;
  float nt;
  float nmult;
  short mute;
} t_squash;


//static float boundrand(float min, float max);
static void *squash_new(t_symbol *msg, int argc, t_atom *argv);
static void squash_mute(t_squash *x, t_floatarg toggle);
//static void squash_assist (t_squash *x, void *b, long msg, long arg, char *dst);
//static void squash_dsp_free(t_squash *x);
static double squash_squat( float *buffer, float thresh, float ratio, float nt, float nmult, int N );
static void squash_mute(t_squash *x, t_floatarg f);
static void squash_thresh(t_squash *x, t_floatarg f);
static void squash_nt(t_squash *x, t_floatarg f);
static void squash_ratio(t_squash *x, t_floatarg f);
static void squash_nmult(t_squash *x, t_floatarg f);
static void squash_mute(t_squash *x, t_floatarg f);
static void squash_free(t_squash *x);
static void squash_dsp(t_squash *x, t_signal **sp);


static void overlapadd( float *I, int N, float *W, float *O, int Nw, int n );
static void fold( float *I, float *W, int Nw, float *O, int N, int n );
/*
void makehanning( float *H, float *A, float *S, int Nw, int N, int I, int odd );
void makehamming( float *H, float *A, float *S, int Nw, int N, int I, int odd );
void makewindows( float *H, float *A, float *S, int Nw, int N, int I );
*/

void squash_tilde_setup(void) {
  squash_class = class_new(gensym("squash~"), (t_newmethod)squash_new,
                           (t_method)squash_free, sizeof(t_squash),0,A_GIMME,0);
  CLASS_MAINSIGNALIN(squash_class, t_squash, x_f);
  class_addmethod(squash_class, (t_method)squash_dsp, gensym("dsp"), A_CANT, 0);
  class_addmethod(squash_class,(t_method)squash_thresh, gensym("thresh"), A_FLOAT,  0);
  class_addmethod(squash_class,(t_method)squash_nt, gensym("nt"), A_FLOAT,  0);
  class_addmethod(squash_class,(t_method)squash_nmult, gensym("nmult"), A_FLOAT,  0);
  class_addmethod(squash_class,(t_method)squash_ratio, gensym("ratio"), A_FLOAT,  0);
  class_addmethod(squash_class,(t_method)squash_mute, gensym("mute"), A_FLOAT,  0);
  potpourri_announce(OBJECT_NAME);
}

void *squash_new(t_symbol *msg, int argc, t_atom *argv)
{

  t_squash *x = (t_squash *)pd_new(squash_class);
  outlet_new(&x->x_obj, gensym("signal"));

  x->D = sys_getblksize();
  x->N = x->D * 4;
  x->Nw = x->N;
  x->N2 = x->N / 2;
  x->incnt = - x->Nw;
  x->Wanal = (float *) getbytes(x->Nw * sizeof(float));
  x->Wsyn = (float *) getbytes(x->Nw *  sizeof(float));
  x->Hwin = (float *) getbytes(x->Nw * sizeof(float));
  x->input = (float *) getbytes(x->Nw * sizeof(float));
  x->output = (float *) getbytes(x->Nw * sizeof(float));
  x->buffer = (float *) getbytes(x->N * sizeof(float));
  lpp_makehanning(x->Hwin, x->Wanal, x->Wsyn, x->Nw, x->N, x->D, 0);
  x->thresh = 0.1;
  x->ratio = 1.;
  x->nt = .0000001;
  x->nmult = 0.1;
  x->mute = 0;
  return x;
}


void squash_free(t_squash *x)
{
    freebytes(x->Wanal, x->Nw * sizeof(float));
    freebytes(x->Wsyn, x->Nw *  sizeof(float));
    freebytes(x->Hwin, x->Nw * sizeof(float));
    freebytes(x->input, x->Nw * sizeof(float));
    freebytes(x->output, x->Nw * sizeof(float));
    freebytes(x->buffer, x->N * sizeof(float));
}
/*
void makewindows( float *H, float *A, float *S, int Nw, int N, int I )

{
  int i ;
  float sum ;

  for ( i = 0 ; i < Nw ; i++ )
    H[i] = A[i] = S[i] = 0.54 - 0.46*cos( TWOPI*i/(Nw - 1) ) ;

  if ( Nw > N ) {
    float x ;

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
    float afac = 2./sum ;
    float sfac = Nw > N ? 1./afac : afac ;
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




void makehamming( float *H, float *A, float *S, int Nw, int N, int I, int odd )

{
  int i;
  float sum ;



  if (odd) {
    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = sqrt(0.54 - 0.46*cos( TWOPI*i/(Nw - 1) ));
  }

  else {

    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = 0.54 - 0.46*cos( TWOPI*i/(Nw - 1) );

  }

  if ( Nw > N ) {
    float x ;

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
    float afac = 2./sum ;
    float sfac = Nw > N ? 1./afac : afac ;
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



void makehanning( float *H, float *A, float *S, int Nw, int N, int I, int odd )
{
  int i;
  float sum ;


  if (odd) {
    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = sqrt(0.5 * (1. + cos(PI + TWOPI * i / (Nw - 1))));
  }

  else {

    for ( i = 0 ; i < Nw ; i++ )
      H[i] = A[i] = S[i] = 0.5 * (1. + cos(PI + TWOPI * i / (Nw - 1)));

  }

  if ( Nw > N ) {
    float x ;

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
    float afac = 2./sum ;
    float sfac = Nw > N ? 1./afac : afac ;
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




*/

void overlapadd( float *I, int N, float *W, float *O, int Nw, int n )
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

void fold( float *I, float *W, int Nw, float *O, int N, int n )

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

void squash_mute(t_squash *x, t_floatarg f) {
  x->mute = (short)f;
}

void squash_thresh(t_squash *x, t_floatarg f) {
  x->thresh = (float)f;
}

void squash_nt(t_squash *x, t_floatarg f) {
  x->nt = (float)f;
}

void squash_ratio(t_squash *x, t_floatarg f) {
  x->ratio = (float)f;
}

void squash_nmult(t_squash *x, t_floatarg f) {
  x->nmult = (float)f;
}

t_int *squash_perform(t_int *w)
{
  t_squash *x = (t_squash *) (w[1]);
  float *in = (t_float *)(w[2]);
  float *out = (t_float *)(w[3]);
  int n = (int) w[4];

  int j;
  float *input = x->input;
  float *output = x->output;
  int D = x->D;
  int Nw = x->Nw;
  int N = x->N;
  float *buffer = x->buffer;
  float *Wanal = x->Wanal;
  float *Wsyn = x->Wsyn;
  float thresh = x->thresh;
  float ratio = x->ratio;
  float nt = x->nt;
  float nmult = x->nmult;

  if(x->mute) {
    memset((void *)out, 0, n * sizeof(float) );
    return w + 5;
  }

  x->incnt += D;

  for ( j = 0 ; j < Nw - D ; j++ )
    input[j] = input[j+D];

  for ( j = Nw-D; j < Nw; j++ ) {
    input[j] = *in++;
  }
  fold( input, Wanal, Nw, buffer, N, x->incnt );
  squash_squat( buffer, thresh, ratio, nt, nmult, Nw );
  overlapadd( buffer, N, Wsyn, output, Nw, x->incnt );

  for ( j = 0; j < D; j++ )
    *out++ = output[j];

  for ( j = 0; j < Nw - D; j++ )
    output[j] = output[j+D];

  for ( j = Nw - D; j < Nw; j++ )
    output[j] = 0.;

  x->incnt = x->incnt % Nw;
  return w + 5;
}

void squash_dsp(t_squash *x, t_signal **sp)
{
  if(sp[0]->s_n != x->D ) {
    x->D = sp[0]->s_n;
    pd_error(0, "blocksize change not implemented yet!");
  } else {
    dsp_add(squash_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, (t_int)sp[0]->s_n);
  }
}

#define DIAG 0

/*  compression/expansion routine! */

double squash_squat( float *buffer, float thresh, float ratio, float nt, float nmult, int N )
{
  register int    i;
  double    rms = 0.;
  float   dbthr;
  register float  mult;

  dbthr = 10. * log10(thresh);


  for ( i=0; i < N; i++ )
    rms += ( *(buffer+i) * *(buffer+i) );

  rms = sqrt(rms/(float)N);
  if (rms < nt && ratio < 1.)
    mult = nmult;
  else
    mult = pow( 10., (( dbthr - (( dbthr - (10. * log10(rms)) ) * ratio) ) / 10.) ) / rms;

  if (DIAG) {
    if (rms <= thresh)
      fprintf(stderr,"below  dbthr: %f  dbrms: %f rms: %f  mult: %f\n", dbthr, (10. * log10(rms)), rms, mult);
    else
      fprintf(stderr,"above  dbthr: %f  dbrms: %f rms: %f  mult: %f\n", dbthr, (10. * log10(rms)), rms, mult);
  }

  for ( i=0; i < N; i++ )
    *(buffer+i) *= mult;

  return rms;
}
