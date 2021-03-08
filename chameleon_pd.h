#include "MSPd.h"
#include "ugens_pd.h"
/* calling codes for DSP modules */
#define BENDY 0
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
#define BITCRUSH 12
#define FLAM1 13
#define SLIDEFLAM 14
#define RINGMOD4 15
#define COMB4 16
#define RINGFEED 17
#define RESONADSR 18
#define STV 19	
//////
#define ROOT2 (1.4142135623730950488)
#define PI2 (6.2831853071795862319959)
#define LOPASS 0
#define HIPASS 1
#define BANDPASS 2
#define COMBFADE (.04 )
#define MAXFILTER 12 /*must be at least as big as see below */
#define ELLIPSE_FILTER_COUNT 11 /*actual number of filters*/
#define MAX_COEF 48
#define MY_MAX 2147483647.0 /* for rand() */
/* maximum number of each DSP unit possible for each pattern (costs significant MBs of memory) */
// #define MAX_DSP_UNITS 12
#define MAX_SLIDECOMB_DELAY 0.05
#define MAX_SLIDEFLAM_DELAY 0.5
#define MAX_MINI_DELAY 0.25
#define FEEDFUNCLEN 8192
#define TONECOMB_MAX_DELAY 0.02
#define BENDY_MINSEG 0.025
#define BENDY_MAXSEG 0.2
#define BENDY_MAXDEL 0.25
#define FLAM1_MAX_DELAY 2.0
#define STV_MAX_DELAY 0.011

/*data types */

typedef struct _flange_unit {
    double *flange_dl1;
    double *flange_dl2;
    int *dv1; // cmix bookkeeping
    int *dv2;
    double phase;
} t_flange_unit;

typedef struct _butterworth_unit {
    int ftype;
    double cutoff;
    double cf;
    double bw;
    double *data1; // two filters
    double *data2;
} t_butterworth_unit;

typedef struct _truncate_unit {
    long counter;
    long state;
    long segsamples;
} t_truncate_unit;

typedef struct _sweepreson_unit {
    double minfreq;
    double maxfreq;
    double bwfac;
    double speed;
    double phase;
    double *q1; // 5 elements
    double *q2;
} t_sweepreson_unit;

typedef struct _slidecomb_unit {
    double *delayline1;
    double *delayline2;
    double start_delay;
    double end_delay;
    double feedback;
    long sample_length;
    long counter;
    int *dv1; // 2 elements
    int *dv2;
} t_slidecomb_unit;

typedef struct _reverb1_unit {
    double *dels; // allocate 4x
    double **alpo1;
    double **alpo2;
    double xnorm;
    LSTRUCT *eel1, *eel2;
    double wet;
    double dry;
    double revtime;
    int nsects;
} t_reverb1_unit;

typedef struct _ellipseme_unit {
    LSTRUCT *eel1, *eel2;
    double xnorm;
    int nsects;
    int filtercode;
} t_ellipseme_unit;

typedef struct _feed1_unit {
    double *func1;
    double *func2;
    double *func3;
    double *func4;
    double mindelay;
    double maxdelay;
    double speed1;
    double speed2;
    double *delayLine1a;
    double *delayLine2a;
    double *delayLine1b;
    double *delayLine2b;
    int *dv1a, *dv2a;
    int *dv1b, *dv2b;
    double funcPhs;
    double duration;
} t_feed1_unit;

typedef struct _flam1_unit {
    double *delayline1;
    double *delayline2;
    double dt;
    long sample_length;
    long counter;
    int *dv1;
    int *dv2;
} t_flam1_unit;

typedef struct _slideflam_unit {
    double *delayline1;
    double *delayline2;
    double dt1,dt2;
    double feedback;
    long sample_length;
    long counter;
    int *dv1;
    int *dv2;
} t_slideflam_unit;


typedef struct _comb4_unit {
    double **combs1;
    double **combs2;
    double revtime;
} t_comb4_unit;

typedef struct _resonfeed_unit {
    double *res1q;
    double *res2q;
    double osc1phs;
    double osc2phs;
    double osc1si;
    double osc2si;
    double *comb1arr;
    double *comb2arr;
} t_resonfeed_unit;

typedef struct _resonadsr_unit {
    double *q1;
    double *q2;
    double phs;
    double si;
    double bwfac;
    CMIXADSR *adsr;
} t_resonadsr_unit;

typedef struct _stv_unit {
    double osc1phs;
    double osc2phs;
    double osc1si;
    double osc2si;
    double *delayline1;
    double *delayline2;
    int *dv1;
    int *dv2;
    double osc1amp;
    double osc2amp;
    double fac1;
    double fac2;
} t_stv_unit;

typedef struct _bendy_unit {
    double val1;
    double val2;
    long counter;
    long segment_samples;
    double min_delay;
    double max_delay;
    double min_speed;
    double max_speed;
    double *delayline1;
    double *delayline2;
    int *dv1;
    int *dv2;
} t_bendy_unit;

typedef struct _slot {
    long pcount;
    double *params;
} t_slot;

typedef struct _chameleon {
    t_object x_obj;
    float x_f;
    float sr; // sampling rate
    long vs;//Max/MSP vector size
    long vecsize; // Max signal vector size
    long pcount; // number of parameters
    double *params; // parameter list
    double fadeout; // fadeout time in sample frames (if truncation)
    t_double *chan1buf; // work vector buffer inside of perform routine
    t_double *chan2buf; // work vector buffer inside of perform routine
    float *odds;// odds for each process happening
    int max_process_per_note;//what it says
    int min_process_per_note;//ditto
    int new_slot;//position for newest note
    double new_gain;//recently assigned gain
    short verbose;//toggle Max window error reporting
    void *listo;
    t_atom *data; // data array to send to list outlet
    char sound_name[256];
    double *distortion_function;
    int distortion_length;
    int set_parameters_flag;
    int recall_parameters_flag;
    t_slot *slots;
    long recall_slot; // slot to be recalled
    /* stuff for chameleon DSP */
    double *sinewave;
    int sinelen;
    short mute;
    double maxdelay;
    double **comb_delay_pool1;
    double **comb_delay_pool2;
    double *delayline1; // delay lines for comber
    double *delayline2;
    double *ringmod_phases; // store phase of ringmod units
    double *ringmod4_phases; // store phase of ringmod4 units
    long max_dsp_units; // set once on initialization to determine maximum number of per-process units
    long membytes; // store the additional memory allocated for the current object
    t_flange_unit *flange_units; // all DSP elements for flanging
    t_butterworth_unit *butterworth_units;// all DSP elements for butterworth filters
    t_truncate_unit *truncate_units;// all DSP elements for truncate
    t_sweepreson_unit *sweepreson_units;
    t_slidecomb_unit *slidecomb_units;
    t_reverb1_unit *reverb1_units;
    t_ellipseme_unit *ellipseme_units;
    t_feed1_unit *feed1_units;
    t_flam1_unit *flam1_units;
    t_comb4_unit *comb4_units;
    t_resonfeed_unit *resonfeed_units;
    t_resonadsr_unit *resonadsr_units;
    t_stv_unit *stv_units;
    t_bendy_unit *bendy_units;
    t_slideflam_unit *slideflam_units;
    
    int *flange_dv1; // CMIX bookkeeping
    int *flange_dv2;
    double flange_phase; // maintain phase for flange oscillator
    double max_flangedelay; // constant to allocate memory for flange DDLs
    
    LSTRUCT *eel; // for ellipse processor
    double *mini_delay[4]; // small delay lines for allpass filter
    double max_mini_delay ;
    double *transfer_function;
    int tf_len; // length of transfer function
    double *feedfunc1;
    double *feedfunc2;
    double *feedfunc3;
    double *feedfunc4;
    int feedfunclen;
    int flamfunc1len;
    double *flamfunc1;
    CMIXADSR *adsr;
    double max_comb_lpt;
    double *reverb_ellipse_data;
    double **ellipse_data;
    double *dcflt;
    CMIXOSC oscar;
    CMIXRESON resies[2];
    double *ratios; // = {9.0/8.0,4.0/3.0, 5.0/4.0, 6.0/5.0, 1.5};
    double *bitcrush_factors;
} t_chameleon;
