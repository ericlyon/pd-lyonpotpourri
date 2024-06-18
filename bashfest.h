#include "MSPd.h"
#include "ugens.h"
#include <string.h>
/* calling codes for DSP modules */
#define TRANSPOSE 0
#define RINGMOD 1
#define FLANGE 2
#define BUTTER 3
#define TRUNCATE 4
#define SWEEPRESON 5
#define COMB 6
#define SLIDECOMB 7
#define REVERB1 8
#define ELLIPSE 9
#define COMPDIST 10
#define FEED1 11
#define RETRO 12
#define FLAM1 13
#define FLAM2 14
#define EXPFLAM 15
#define COMB4 16
#define RINGFEED 17
#define RESONADSR 18
#define STV 19
//////
#define ROOT2 (1.4142135623730950488)
#define PI2 (6.2831853071795862319959)
// #define BUFFER_SIZE (1<<15)
#define LOPASS 0
#define HIPASS 1
#define BANDPASS 2
#define COMBFADE (.04 )
#define MAXFILTER 12 /*must be at least as big as see below */
#define ELLIPSE_FILTER_COUNT 11 /*actual number of filters*/
#define MAX_COEF 48
#define MY_MAX 2147483647.0 /* for rand() */
/*data types */

typedef struct
{
  t_float *data;//contains cycle data
  int len;//length of array
  int p;//position pointer
} t_cycle;

typedef struct
{
  long phase; // current phase in frames
  double phasef; // current phase in frames
  t_float gain; // gain for this note
  t_float gainL;// left gain
  t_float gainR;// right gain
  short status;// status of this event slot
  t_float *workbuffer;//sample processing space (both input and output)
  t_float *inbuf;//pointer to input part of workbuffer
  t_float *outbuf;//pointer to output part of workbuffer
  int in_start;// location in workbuffer to read from input
  int out_start;// location in workbuffer to write output
  int sample_frames;//actual size in frames of sample, which changes if it gets bigger
  int countdown;//latency counter before we actually start sending out samples
  int out_channels; //number of channels per frame of output
  short completed;//did the defer call do its thing?


} t_event;

typedef struct _bashfest
{
  t_object x_obj;
  t_float x_f;
  t_float sr; // sampling rate
  t_symbol *wavename; // name of waveform buffer
  short hosed; // buffers are bad
  t_float fadeout; // fadeout time in sample frames (if truncation)
  t_float sync; // input from groove sync signal
  t_float increment; // read increment
  int most_recent_event; // position in array where last note was initiated
  long b_nchans; // channels of buffer
  long b_valid; // state of buffer
  long b_frames; // number of frames in sample buffer
  t_word *b_samples; // pointer samples in buffer
  int overlap_max; // max number of simultaneous plays
  t_event *events; //note attacks
  int active_events; // how many currently activated notes?
  int buf_samps;//total samples in workbuffer
  int halfbuffer;//buf_samps / 2
  int buf_frames;// number of sample frames in workbuffer
  int latency_samples;// amount of samples to count down before playing sample
  t_float *params; // parameter list
  t_float *odds;// odds for each process happening
  int max_process_per_note;//what it says
  int min_process_per_note;//ditto
  int new_slot;//position for newest note
  t_float new_gain;//recently assigned gain
  short verbose;//toggle Max window error reporting
  t_float work_buffer_size;// size in ms of work buffers
  t_cycle tcycle;//contains an optional transposition cycle
  short block_dsp;//flag to turn off all dsp and play straight from MSP buffer
  short sound_lock;//keep current processed sound in buffer
  short grab;//flag to copy immediate processed buffer into MSP buffer
  char sound_name[256];
  t_float *trigger_vec;//stores incoming trigger vectors
  int vs;//Max/MSP vector size

  /* stuff for bashfest DSP */
  t_float *sinewave;
  int sinelen;
  short mute;
  t_float maxdelay;
  t_float *delayline1;
  t_float *delayline2;
  LSTRUCT *eel; // for ellipse processor
  t_float *mini_delay[4]; // small delay lines for allpass filter
  t_float max_mini_delay ;
  t_float *transfer_function;
  int tf_len; // length of transfer function
  t_float *feedfunc1;
  t_float *feedfunc2;
  t_float *feedfunc3;
  t_float *feedfunc4;
  int feedfunclen;
  int flamfunc1len;
  t_float *flamfunc1;
  CMIXCOMB *combies;
  CMIXADSR *adsr;
  t_float max_comb_lpt;
  t_float *reverb_ellipse_data;
  t_float **ellipse_data;
  t_float *dcflt;
  CMIXOSC oscar;
  CMIXRESON resies[2];

} t_bashfest;



/*function prototypes*/
void lpp_putsine (t_float *arr, int len);
t_float lpp_boundrand(t_float min, t_float max);
void lpp_mycombset(t_float loopt,t_float rvt,int init,t_float *a,t_float srate);
t_float lpp_mycomb(t_float samp,t_float *a);
void lpp_setweights(t_float *a, int len);
void lpp_delset2(t_float *a,int *l,t_float xmax, t_float srate);
void lpp_delput2(t_float x,t_float *a,int *l);
t_float lpp_dliget2(t_float *a,t_float dwait,int *l,t_float srate);
void lpp_butterLopass( t_float *in, t_float *out, t_float cutoff, int frames, int channels, t_float srate);
void lpp_butterBandpass(t_float *in, t_float *out,  t_float center, t_float bandwidth, int frames,int  channels, t_float srate);
void lpp_butterHipass(t_float *in, t_float *out,  t_float cutoff, int frames,int channels, t_float srate);
void lpp_butset(t_float *a);
void lpp_lobut(t_float *a, t_float cutoff,t_float SR);
void lpp_hibut(t_float *a, t_float cutoff, t_float SR);
void lpp_bpbut(t_float *a, t_float formant, t_float bandwidth, t_float SR);
void lpp_butter_filter(t_float *in,t_float *out,t_float *a, int frames, int channels, int channel);
void lpp_rsnset2(t_float cf,t_float bw,t_float scl,t_float xinit,t_float *a,t_float srate);
t_float lpp_reson(t_float x,t_float *a);

void lpp_ellipset(t_float *list, LSTRUCT *eel, int  *nsects, t_float *xnorm);
t_float lpp_ellipse(t_float x, LSTRUCT *eel, int nsects, t_float xnorm);
t_float lpp_allpass(t_float samp,t_float *a);
void lpp_init_reverb_data(t_float *a);
void lpp_init_ellipse_data(t_float **a);

void lpp_setExpFlamFunc(t_float *arr, int flen, t_float v1,t_float v2,t_float alpha);
void lpp_setflamfunc1(t_float *arr, int flen);
void lpp_funcgen1(t_float *outArray, int outlen, t_float duration, t_float outMin, t_float outMax,
              t_float speed1, t_float speed2, t_float gain1, t_float gain2, t_float *phs1, t_float *phs2,
              t_float *sine, int sinelen);
void lpp_normtab(t_float *inarr,t_float *outarr, t_float min, t_float max, int len);
t_float lpp_mapp(t_float in,t_float imin,t_float imax,t_float omin,t_float omax);
t_float lpp_oscil(t_float amp,t_float si,t_float *farray,int len,t_float *phs);
void lpp_set_dcflt(t_float *a);

void lpp_set_distortion_table(t_float *arr, t_float cut, t_float max, int len);
t_float lpp_dlookup(t_float samp,t_float *arr,int len);
void lpp_do_compdist(t_float *in,t_float *out,int sampFrames,int nchans,int channel,
                 t_float cutoff,t_float maxmult,int lookupflag,t_float *table,int range,t_float bufMaxamp);
t_float lpp_getmaxamp(t_float *arr, int len);
void lpp_buildadsr(CMIXADSR *a);
/*bashfest dsp functions */
void lpp_feed1(t_float *inbuf, t_float *outbuf, int in_frames, int out_frames,int channels, t_float *functab1,
           t_float *functab2,t_float *functab3,t_float *functab4,int funclen,
           t_float duration, t_float maxDelay, t_bashfest *x);
void lpp_reverb1me(t_float *in, t_float *out, int inFrames, int out_frames, int nchans,
               int channel, t_float revtime, t_float dry, t_bashfest *x);
void lpp_killdc( t_float *inbuf, int in_frames, int channels, t_bashfest *x);
