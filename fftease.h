#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*
  The new improved fftease.h
*/

#define getbytes t_getbytes
#define freebytes t_freebytes
#define resizebytes t_resizebytes

#define FFTEASE_ANNOUNCEMENT "- a member of FFTease 2.52"

#ifndef PI
#define PI 3.141592653589793115997963468544185161590576171875
#endif

#ifndef TWOPI
#define TWOPI 6.28318530717958623199592693708837032318115234375
#endif

#define MAX_N (1073741824)
#define MAX_N2 (MAX_N/2)
#define MAX_Nw (MAX_N * 4)

void lpp_convert(float *S, float *C, int N2, float *lastphase, float fundamental, float factor );
void lpp_unconvert( float *C, float *S, int N2, float *lastphase, float fundamental,  float factor );
void lpp_rfft( float *x, int N, int forward );
void lpp_cfft( float *x, int NC, int forward );
void lpp_bitreverse( float *x, int N );
void lpp_fold( float *I, float *W, int Nw, float *O, int N, int n );
void lpp_init_rdft(int n, int *ip, float *w);
void lpp_rdft(int n, int isgn, float *a, int *ip, float *w);
void lpp_bitrv2(int n, int *ip, float *a);
void lpp_cftsub(int n, float *a, float *w);
void lpp_rftsub(int n, float *a, int nc, float *c);
void lpp_makewt(int nw, int *ip, float *w);
void lpp_makect(int nc, int *ip, float *c);
void lpp_leanconvert( float *S, float *C, int N2 );
void lpp_leanunconvert( float *C, float *S, int N2 );
void lpp_makewindows( float *H, float *A, float *S, int Nw, int N, int I );
void lpp_makehamming( float *H, float *A, float *S, int Nw, int N, int I,int odd );
void lpp_makehanning( float *H, float *A, float *S, int Nw, int N, int I,int odd );
void lpp_overlapadd( float *I, int N, float *W, float *O, int Nw, int n );
void lpp_bloscbank( float *S, float *O, int D, float iD, float *lf, float *la,
                float *bindex, float *tab, int len, float synt, int lo, int hi );

float lpp_randf( float min, float max );
int lpp_randi( int min, int max );
int lpp_power_of_two(int test);


//void freebytes2(void *fatso, size_t nbytes);
//void *getbytes2(size_t nbytes);
//void *resizebytes2(void *old, size_t oldsize, size_t newsize);
void lpp_limit_fftsize(int *N, int *Nw, char *OBJECT_NAME);

/* THE FUNCTIONS */
