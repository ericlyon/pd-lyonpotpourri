#include "MSPd.h"

/* Pd 3.0 Version */


#define MAX_N 16384
#define MAX_N2 8192
#define MAX_Nw 16384

static void convert(float *S, float *C, int N2, float *lastphase, float fundamental, float factor );
static void fold( float *I, float *W, int Nw, float *O, int N, int n );
static void init_rdft(int n, int *ip, float *w);
static void rdft(int n, int isgn, float *a, int *ip, float *w);
static void bitrv2(int n, int *ip, float *a);
static void cftsub(int n, float *a, float *w);
static void rftsub(int n, float *a, int nc, float *c);
static void makewt(int nw, int *ip, float *w);
static void makect(int nc, int *ip, float *c);
static void makehanning( float *H, float *A, float *S, int Nw, int N, int I,int odd );
static int power_of_two(int test);

static t_class *magfreq_analysis_class;

#define OBJECT_NAME "magfreq_analysis~"

typedef struct _magfreq_analysis
{
    t_object x_obj;
    float x_f;
    float R;
    int N;
    int N2;
    int Nw;
    int Nw2;
    int D;
    int i;
    int inCount;
    float *Wanal;
    float *Wsyn;
    float *input;
    float *Hwin;
    float *buffer;
    float *channel;
    float *output;
    // for convert
    float *c_lastphase_in;
    float *c_lastphase_out;
    float c_fundamental;
    float c_factor_in;
    float c_factor_out;
    // for oscbank
    int NP;
    float P;
    int L;
    int first;
    float Iinv;
    float *lastamp;
    float *lastfreq;
    float *index;
    float *table;
    float myPInc;
    float ffac;
    //
    float lofreq;
    float hifreq;
    int lo_bin;
    int hi_bin;
    float topfreq;
    float synt;
    // for fast fft
    float mult;
    float *trigland;
    int *bitshuffle;
    //
    int bypass_state;
    int pitch_connected;
    int synt_connected;
    int overlap;
    int winfac;
    short mute;
} t_magfreq_analysis;

static void *magfreq_analysis_new(t_symbol *s, int argc, t_atom *argv);
static t_int *magfreq_analysis_perform(t_int *w);
static void magfreq_analysis_dsp(t_magfreq_analysis *x, t_signal **sp);
// static void magfreq_analysis_bypass(t_magfreq_analysis *x, t_floatarg state);
// static void magfreq_analysis_float(t_magfreq_analysis *x, double f);
static void magfreq_analysis_free(t_magfreq_analysis *x);
static void magfreq_analysis_mute(t_magfreq_analysis *x, t_floatarg tog);
static void magfreq_analysis_init(t_magfreq_analysis *x, short initialized);
//static void magfreq_analysis_lowfreq(t_magfreq_analysis *x, t_floatarg f);
//static void magfreq_analysis_highfreq(t_magfreq_analysis *x, t_floatarg f);
//static void magfreq_analysis_overlap(t_magfreq_analysis *x, t_floatarg o);
//static void magfreq_analysis_winfac(t_magfreq_analysis *x, t_floatarg f);
static void magfreq_analysis_fftinfo(t_magfreq_analysis *x);;


void magfreq_analysis_tilde_setup(void)
{
    magfreq_analysis_class = class_new(gensym("magfreq_analysis~"), (t_newmethod)magfreq_analysis_new,
                                       (t_method)magfreq_analysis_free ,sizeof(t_magfreq_analysis), 0,A_GIMME,0);
    CLASS_MAINSIGNALIN(magfreq_analysis_class, t_magfreq_analysis, x_f);
    class_addmethod(magfreq_analysis_class, (t_method)magfreq_analysis_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(magfreq_analysis_class, (t_method)magfreq_analysis_mute, gensym("mute"), A_DEFFLOAT,0);
//    class_addmethod(magfreq_analysis_class, (t_method)magfreq_analysis_highfreq, gensym("highfreq"), A_DEFFLOAT,0);
//    class_addmethod(magfreq_analysis_class, (t_method)magfreq_analysis_lowfreq, gensym("lowfreq"), A_DEFFLOAT,0);
    class_addmethod(magfreq_analysis_class, (t_method)magfreq_analysis_fftinfo, gensym("fftinfo"),0);
    potpourri_announce(OBJECT_NAME);
}


void magfreq_analysis_mute(t_magfreq_analysis *x, t_floatarg tog)
{
    x->mute = (short)tog;
}

void magfreq_analysis_overlap(t_magfreq_analysis *x, t_floatarg f)
{
    int i = (int) f;
    if(!power_of_two(i)) {
        pd_error(0, "%f is not a power of two",f);
        return;
    }
    x->overlap = i;
    magfreq_analysis_init(x,1);
}

void magfreq_analysis_winfac(t_magfreq_analysis *x, t_floatarg f)
{
    int i = (int)f;
    
    if(!power_of_two(i)) {
        pd_error(0, "%f is not a power of two",f);
        return;
    }
    x->winfac = i;
    magfreq_analysis_init(x,2);
}

void magfreq_analysis_fftinfo(t_magfreq_analysis *x)
{
    if( ! x->overlap ) {
        post("zero overlap!");
        return;
    }
    post("%s: FFT size %d, hopsize %d, windowsize %d", OBJECT_NAME, x->N, x->N/x->overlap, x->Nw);
    post("sample rate: %f", x->R);
    post("fundamental analysis frequency %f", x->c_fundamental);
}

void magfreq_analysis_free(t_magfreq_analysis *x ) {
    freebytes(x->c_lastphase_in,(MAX_N2+1) * sizeof(float));
    freebytes(x->c_lastphase_out,(MAX_N2+1) * sizeof(float));
    freebytes(x->trigland,MAX_N * 2 * sizeof( float ));
    freebytes(x->bitshuffle,MAX_N * 2 * sizeof( int ));
    freebytes(x->Wanal,(MAX_Nw) * sizeof(float));
    freebytes(x->Wsyn,(MAX_Nw) * sizeof(float));
    freebytes(x->input,MAX_Nw * sizeof(float));
    freebytes(x->Hwin,(MAX_Nw) * sizeof(float));
    freebytes(x->buffer,MAX_N * sizeof(float));
    freebytes(x->channel,(MAX_N+2) * sizeof(float));
    freebytes(x->output,MAX_Nw * sizeof(float));
    freebytes(x->lastamp,(MAX_N+1) * sizeof(float));
    freebytes(x->lastfreq,(MAX_N+1) * sizeof(float));
    freebytes(x->index,(MAX_N+1) * sizeof(float));
    freebytes(x->table,x->L * sizeof(float));
}

void magfreq_analysis_highfreq(t_magfreq_analysis *x, t_floatarg f)
{
    float curfreq;
    
    if(f < x->lofreq) {
        pd_error(0, "current minimum is %f",x->lofreq);
        return;
    }
    if(f > x->R/2 ) {
        f = x->R/2;
    }
    x->hifreq = f;
    x->hi_bin = 1;
    curfreq = 0;
    while(curfreq < x->hifreq) {
        ++(x->hi_bin);
        curfreq += x->c_fundamental;
    }
}

void magfreq_analysis_lowfreq(t_magfreq_analysis *x, t_floatarg f)
{
    float curfreq;
    
    if(f > x->hifreq) {
        pd_error(0, "current maximum is %f",x->lofreq);
        return;
    }
    if(f < 0 ) {
        f = 0;
    }
    x->lofreq = f;
    x->lo_bin = 0;
    curfreq = 0;
    while( curfreq < x->lofreq ) {
        ++(x->lo_bin);
        curfreq += x->c_fundamental ;
    }
}


void magfreq_analysis_init(t_magfreq_analysis *x, short initialized)
{
    int i;
    float curfreq;
    x->R = sys_getsr();
    x->D = sys_getblksize();
   
    if(!x->R){//temp init if MSP functions returned zero
        x->R = 48000.0;
    }
    if(!x->D){
        x->D = 64;
    }
    if(x->P <= 0){
        x->P = 1.0;
    }
    if(!power_of_two(x->overlap)){
        x->overlap = 2;
    }
    if(!power_of_two(x->winfac)){
        x->winfac = 2;
    }
    x->N = x->D * x->overlap;
    x->Nw = x->N * x->winfac;
    x->N2 = x->N / 2;
    x->Nw2 = x->Nw / 2;
    x->inCount = -(x->Nw);
    x->bypass_state = 0;
    x->mult = 1. / (float) x->N;
    x->pitch_connected = 0;
    x->synt_connected = 0;
    x->L = 8192 ;
    x->c_fundamental =  x->R/(float)( (x->N2)<<1 );
    x->c_factor_in =  x->R/((float)x->D * TWOPI);
    x->c_factor_out = TWOPI * (float)  x->D / (float) x->R;
    x->Iinv = 1./(float)x->D;
    x->myPInc = x->P*x->L/x->R;
    x->ffac = x->P * PI/(float)x->N;

    if(!initialized) {
        x->Wanal = (float *) getbytes( (MAX_Nw) * sizeof(float));
        x->Wsyn = (float *) getbytes( (MAX_Nw) * sizeof(float));
        x->Hwin = (float *) getbytes( (MAX_Nw) * sizeof(float));
        x->input = (float *) getbytes(MAX_Nw * sizeof(float) );
        x->output = (float *) getbytes(MAX_Nw * sizeof(float) );
        x->buffer = (float *) getbytes(MAX_N * sizeof(float) );
        x->channel = (float *) getbytes( (MAX_N+2) * sizeof(float) );
        x->bitshuffle = (int *) getbytes(MAX_N * 2 * sizeof( int ) );
        x->trigland = (float *) getbytes(MAX_N * 2 * sizeof( float ) );
        x->c_lastphase_in = (float *) getbytes( (MAX_N2+1) * sizeof(float) );
        x->c_lastphase_out = (float *) getbytes( (MAX_N2+1) * sizeof(float) );
        x->lastamp = (float *) getbytes( (MAX_N+1) * sizeof(float) );
        x->lastfreq = (float *) getbytes( (MAX_N+1) * sizeof(float) );
        x->index = (float *) getbytes( (MAX_N+1) * sizeof(float) );
        x->table = (float *) getbytes( x->L * sizeof(float) );
        x->P = 1.0;
        x->ffac = x->P * PI/(float)MAX_N;
        x->mute = 0;
        //    x->threshgen = .0001;
    }
/*
    memset((char *)x->input,0,x->Nw * sizeof(float));
    memset((char *)x->output,0,x->Nw * sizeof(float));
    memset((char *)x->c_lastphase_in,0,(x->N2+1) * sizeof(float));
    memset((char *)x->c_lastphase_out,0,(x->N2+1) * sizeof(float));
    memset((char *)x->lastamp,0,(x->N+1) * sizeof(float));
    memset((char *)x->lastfreq,0,(x->N+1) * sizeof(float));
    memset((char *)x->index,0,(x->N+1) * sizeof(float));
*/
    for ( i = 0; i < x->L; i++ ) {
        x->table[i] = (float) x->N * cos((float)i * TWOPI / (float)x->L);
    }
    init_rdft( x->N, x->bitshuffle, x->trigland);
    
    makehanning( x->Hwin, x->Wanal, x->Wsyn, x->Nw, x->N, x->D, 0);

    if( x->hifreq < x->c_fundamental ) {
        x->hifreq = 3000.0 ;
    }
    x->hi_bin = 1;
    curfreq = 0;

    while( curfreq < x->hifreq ) {
        ++(x->hi_bin);
        curfreq += x->c_fundamental ;
    }
    
    
    x->lo_bin = 0;
    curfreq = 0;
    while( curfreq < x->lofreq ) {
        ++(x->lo_bin);
        curfreq += x->c_fundamental ;
    }

}

void *magfreq_analysis_new(t_symbol *s, int argc, t_atom *argv)
{
    
    t_magfreq_analysis *x = (t_magfreq_analysis *)pd_new(magfreq_analysis_class);
    outlet_new(&x->x_obj, gensym("signal"));
    outlet_new(&x->x_obj, gensym("signal"));
    outlet_new(&x->x_obj, gensym("signal"));
    
 //   x->lofreq = atom_getfloatarg(0,argc,argv);
 //   x->hifreq = atom_getfloatarg(1,argc,argv);
    x->overlap = atom_getfloatarg(0,argc,argv);
    x->winfac = atom_getfloatarg(1,argc,argv);

    // these might not actually have any effect on analysis:
    x->lofreq = 0;
    x->hifreq = 4000;
    /*
    if(x->lofreq <0 || x->lofreq> 22050)
        x->lofreq = 0;
    if(x->hifreq <50 || x->hifreq> 22050)
        x->hifreq = 4000;
    */
    x->P = 1.0;
    if(!power_of_two(x->overlap)) {
        x->overlap = 4;
    }
    if(!power_of_two(x->winfac)) {
        x->winfac = 2;
    }
    x->R = sys_getsr();
    x->D = sys_getblksize();
    
    magfreq_analysis_init(x,0);
    
    return x;
}

t_int *magfreq_analysis_perform(t_int *w)
{
    int   j, in,on;
    int    amp,freq,chan;
    
    t_magfreq_analysis *x = (t_magfreq_analysis *) (w[1]);
    t_float *inbuf = (t_float *)(w[2]);
    t_float *magnitude_vec = (t_float *)(w[3]);
    t_float *frequency_vec = (t_float *)(w[4]);
    t_float *index_vec = (t_float *)(w[5]);
    int n = (int) w[6];
    
    int D = x->D;
    int Nw = x->Nw;
    int N = x->N ;
    int N2 = x-> N2;
    float fundamental = x->c_fundamental;
    float factor_in =  x->c_factor_in;
    int *bitshuffle = x->bitshuffle;
    float *trigland = x->trigland;
    float *lastphase_in = x->c_lastphase_in;
    
    
    float *Wanal = x->Wanal;
    float *input = x->input;;
    float *buffer = x->buffer;
    float *channel = x->channel;
    in = on = x->inCount ;
    
    
    if(x->mute) {
        for( j = 0; j < n; j++ ) {
            *magnitude_vec++ = 0;
            *frequency_vec++ = 0;
            *index_vec++ = j;
        }
        return w+7;
    }
    
    
    if (x->bypass_state) {
        for( j = 0; j < n; j++ ) {
            *magnitude_vec++ = 0;
            *frequency_vec++ = 0;
            *index_vec++ = 0;
        }
        return w+7;
    }
    
    in = on = x->inCount ;
    
    in += D;
    //    on += I;
    
    for ( j = 0 ; j < (Nw - D) ; j++ ) {
        input[j] = input[j+D];
    }
    for ( j = (Nw-D); j < Nw; j++) {
        input[j] = *inbuf++;
    }
    
    fold( input, Wanal, Nw, buffer, N, in );
    rdft( N, 1, buffer, bitshuffle, trigland );
    convert( buffer, channel, N2, lastphase_in, fundamental, factor_in );
    
    
    // start osc bank
    
    for ( chan = 0; chan < n; chan++ ) {
        
        freq = ( amp = ( chan << 1 ) ) + 1;
        frequency_vec[chan] = channel[freq];
        magnitude_vec[chan] = channel[amp];
        
        index_vec[chan] = chan;
    }
    
    
    
    // restore state variables
    x->inCount = in % Nw;
    return w+7;
}

void magfreq_analysis_bypass(t_magfreq_analysis *x, t_floatarg state)
{
    x->bypass_state = state;
}


void magfreq_analysis_dsp(t_magfreq_analysis *x, t_signal **sp)
{
    
    if(x->D != sp[0]->s_n || x->R != sp[0]->s_sr ) {
        x->D = sp[0]->s_n;
        x->R = sp[0]->s_sr;
        magfreq_analysis_init(x,1);
    }
 
    dsp_add(magfreq_analysis_perform, 6, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
            (t_int)sp[0]->s_n);

}

// IMPORTED FUNCTIONS

int power_of_two(int test)
{
  int limit = 8192;
  int compare = 1;
  //  post("testing what we thing is an int:%d",test);
  do {
    if(test == compare){
      //      post("good power of 2 found!");
      return 1;
    }
    compare *= 2;
  } while (compare <= limit);
  
  return 0;
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


void convert(float *S, float *C, int N2, float *lastphase, float fundamental, float factor )
{
  float     phase,
        phasediff;
  int         real,
        imag,
        amp,
        freq;
  float     a,
        b;
  int         i;
    for ( i = 0; i <= N2; i++ ) {
      imag = freq = ( real = amp = i<<1 ) + 1;
      a = ( i == N2 ? S[1] : S[real] );
      b = ( i == 0 || i == N2 ? 0. : S[imag] );

      C[amp] = hypot( a, b );
      if ( C[amp] == 0. )
    phasediff = 0.;
      else {
    phasediff = ( phase = -atan2( b, a ) ) - lastphase[i];
    lastphase[i] = phase;
    
    while ( phasediff > PI )
      phasediff -= TWOPI;
    while ( phasediff < -PI )
      phasediff += TWOPI;
      }
      C[freq] = phasediff*factor + i*fundamental;
    }
}
// more libraries
void init_rdft(int n, int *ip, float *w)
{

  int    nw,
    nc;

  void    makewt(int nw, int *ip, float *w);
  void    makect(int nc, int *ip, float *c);

  nw = n >> 2;
  makewt(nw, ip, w);

  nc = n >> 2;
  makect(nc, ip, w + nw);

  return;
}


void rdft(int n, int isgn, float *a, int *ip, float *w)
{

  int        j,
        nw,
        nc;

  float        xi;

  void        bitrv2(int n, int *ip, float *a),
        cftsub(int n, float *a, float *w),
        rftsub(int n, float *a, int nc, float *c);

    
  nw = ip[0];
  nc = ip[1];
  
  if (isgn < 0) {
    a[1] = 0.5 * (a[1] - a[0]);
    a[0] += a[1];

    for (j = 3; j <= n - 1; j += 2) {
      a[j] = -a[j];
    }

    if (n > 4) {
      rftsub(n, a, nc, w + nw);
      bitrv2(n, ip + 2, a);
    }

    cftsub(n, a, w);

    for (j = 1; j <= n - 1; j += 2) {
      a[j] = -a[j];
    }
  }

  else {

    if (n > 4) {
      bitrv2(n, ip + 2, a);
    }

    cftsub(n, a, w);

    if (n > 4) {
      rftsub(n, a, nc, w + nw);
    }

    xi = a[0] - a[1];
    a[0] += a[1];
    a[1] = xi;
  }
}


void bitrv2(int n, int *ip, float *a)
{
  int j, jj1, k, k1, l, m, m2;
  float xr, xi;
    
  ip[0] = 0;
  l = n;
  m = 1;

  while ((m << 2) < l) {
    l >>= 1;
    for (j = 0; j <= m - 1; j++) {
      ip[m + j] = ip[j] + l;
    }
    m <<= 1;
  }

  if ((m << 2) > l) {

    for (k = 1; k <= m - 1; k++) {

      for (j = 0; j <= k - 1; j++) {
    jj1 = (j << 1) + ip[k];
    k1 = (k << 1) + ip[j];
    xr = a[jj1];
    xi = a[jj1 + 1];
    a[jj1] = a[k1];
    a[jj1 + 1] = a[k1 + 1];
    a[k1] = xr;
    a[k1 + 1] = xi;
      }
    }
  }

  else {
    m2 = m << 1;

    for (k = 1; k <= m - 1; k++) {

      for (j = 0; j <= k - 1; j++) {
    jj1 = (j << 1) + ip[k];
    k1 = (k << 1) + ip[j];
    xr = a[jj1];
    xi = a[jj1 + 1];
    a[jj1] = a[k1];
    a[jj1 + 1] = a[k1 + 1];
    a[k1] = xr;
    a[k1 + 1] = xi;
    jj1 += m2;
    k1 += m2;
    xr = a[jj1];
    xi = a[jj1 + 1];
    a[jj1] = a[k1];
    a[jj1 + 1] = a[k1 + 1];
    a[k1] = xr;
    a[k1 + 1] = xi;
      }
    }
  }
}


void cftsub(int n, float *a, float *w)
{
  int j, jj1, j2, j3, k, k1, ks, l, m;
  float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
  l = 2;

  while ((l << 1) < n) {
    m = l << 2;

    for (j = 0; j <= l - 2; j += 2) {
      jj1 = j + l;
      j2 = jj1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[jj1];
      x0i = a[j + 1] + a[jj1 + 1];
      x1r = a[j] - a[jj1];
      x1i = a[j + 1] - a[jj1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      a[j2] = x0r - x2r;
      a[j2 + 1] = x0i - x2i;
      a[jj1] = x1r - x3i;
      a[jj1 + 1] = x1i + x3r;
      a[j3] = x1r + x3i;
      a[j3 + 1] = x1i - x3r;
    }

    if (m < n) {
      wk1r = w[2];

      for (j = m; j <= l + m - 2; j += 2) {
    jj1 = j + l;
    j2 = jj1 + l;
    j3 = j2 + l;
    x0r = a[j] + a[jj1];
    x0i = a[j + 1] + a[jj1 + 1];
    x1r = a[j] - a[jj1];
    x1i = a[j + 1] - a[jj1 + 1];
    x2r = a[j2] + a[j3];
    x2i = a[j2 + 1] + a[j3 + 1];
    x3r = a[j2] - a[j3];
    x3i = a[j2 + 1] - a[j3 + 1];
    a[j] = x0r + x2r;
    a[j + 1] = x0i + x2i;
    a[j2] = x2i - x0i;
    a[j2 + 1] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[jj1] = wk1r * (x0r - x0i);
    a[jj1 + 1] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[j3] = wk1r * (x0i - x0r);
    a[j3 + 1] = wk1r * (x0i + x0r);
      }

      k1 = 1;
      ks = -1;

      for (k = (m << 1); k <= n - m; k += m) {
    k1++;
    ks = -ks;
    wk1r = w[k1 << 1];
    wk1i = w[(k1 << 1) + 1];
    wk2r = ks * w[k1];
    wk2i = w[k1 + ks];
    wk3r = wk1r - 2 * wk2i * wk1i;
    wk3i = 2 * wk2i * wk1r - wk1i;

    for (j = k; j <= l + k - 2; j += 2) {
      jj1 = j + l;
      j2 = jj1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[jj1];
      x0i = a[j + 1] + a[jj1 + 1];
      x1r = a[j] - a[jj1];
      x1i = a[j + 1] - a[jj1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      x0r -= x2r;
      x0i -= x2i;
      a[j2] = wk2r * x0r - wk2i * x0i;
      a[j2 + 1] = wk2r * x0i + wk2i * x0r;
      x0r = x1r - x3i;
      x0i = x1i + x3r;
      a[jj1] = wk1r * x0r - wk1i * x0i;
      a[jj1 + 1] = wk1r * x0i + wk1i * x0r;
      x0r = x1r + x3i;
      x0i = x1i - x3r;
      a[j3] = wk3r * x0r - wk3i * x0i;
      a[j3 + 1] = wk3r * x0i + wk3i * x0r;
    }
      }
    }

    l = m;
  }

  if (l < n) {

    for (j = 0; j <= l - 2; j += 2) {
      jj1 = j + l;
      x0r = a[j] - a[jj1];
      x0i = a[j + 1] - a[jj1 + 1];
      a[j] += a[jj1];
      a[j + 1] += a[jj1 + 1];
      a[jj1] = x0r;
      a[jj1 + 1] = x0i;
    }
  }
}


void rftsub(int n, float *a, int nc, float *c)
{
  int j, k, kk, ks;
  float wkr, wki, xr, xi, yr, yi;
    
  ks = (nc << 2) / n;
  kk = 0;

  for (k = (n >> 1) - 2; k >= 2; k -= 2) {
    j = n - k;
    kk += ks;
    wkr = 0.5 - c[kk];
    wki = c[nc - kk];
    xr = a[k] - a[j];
    xi = a[k + 1] + a[j + 1];
    yr = wkr * xr - wki * xi;
    yi = wkr * xi + wki * xr;
    a[k] -= yr;
    a[k + 1] -= yi;
    a[j] += yr;
    a[j + 1] -= yi;
  }
}


void makewt(int nw, int *ip, float *w)
{
    void bitrv2(int n, int *ip, float *a);
    int nwh, j;
    float delta, x, y;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j <= nwh - 2; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}


void makect(int nc, int *ip, float *c)
{
    int nch, j;
    float delta;
    
    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = 0.5;
        c[nch] = 0.5 * cos(delta * nch);
        for (j = 1; j <= nch - 1; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}
