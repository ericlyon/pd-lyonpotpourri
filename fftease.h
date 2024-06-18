#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef PD
# include <m_pd.h>
#else
# define t_float float
#endif


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

void lpp_convert(t_float *S, t_float *C, int N2, t_float *lastphase, t_float fundamental, t_float factor );
void lpp_unconvert( t_float *C, t_float *S, int N2, t_float *lastphase, t_float fundamental,  t_float factor );
void lpp_rfft( t_float *x, int N, int forward );
void lpp_cfft( t_float *x, int NC, int forward );
void lpp_bitreverse( t_float *x, int N );
void lpp_fold( t_float *I, t_float *W, int Nw, t_float *O, int N, int n );
void lpp_init_rdft(int n, int *ip, t_float *w);
void lpp_rdft(int n, int isgn, t_float *a, int *ip, t_float *w);
void lpp_bitrv2(int n, int *ip, t_float *a);
void lpp_cftsub(int n, t_float *a, t_float *w);
void lpp_rftsub(int n, t_float *a, int nc, t_float *c);
void lpp_makewt(int nw, int *ip, t_float *w);
void lpp_makect(int nc, int *ip, t_float *c);
void lpp_leanconvert( t_float *S, t_float *C, int N2 );
void lpp_leanunconvert( t_float *C, t_float *S, int N2 );
void lpp_makewindows( t_float *H, t_float *A, t_float *S, int Nw, int N, int I );
void lpp_makehamming( t_float *H, t_float *A, t_float *S, int Nw, int N, int I,int odd );
void lpp_makehanning( t_float *H, t_float *A, t_float *S, int Nw, int N, int I,int odd );
void lpp_overlapadd( t_float *I, int N, t_float *W, t_float *O, int Nw, int n );
void lpp_bloscbank( t_float *S, t_float *O, int D, t_float iD, t_float *lf, t_float *la,
                t_float *bindex, t_float *tab, int len, t_float synt, int lo, int hi );

t_float lpp_randf( t_float min, t_float max );
int lpp_randi( int min, int max );
int lpp_power_of_two(int test);


//void freebytes2(void *fatso, size_t nbytes);
//void *getbytes2(size_t nbytes);
//void *resizebytes2(void *old, size_t oldsize, size_t newsize);
void lpp_limit_fftsize(int *N, int *Nw, char *OBJECT_NAME);

/* THE FUNCTIONS */
