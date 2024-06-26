#include "MSPd.h"

#define MAXGRAINS (512) // just for present to get lower overhead

#define MAXSCALE (8192)
#define OBJECT_NAME "granule~"


static t_class *granule_class;


typedef struct {
  t_float amplitude;
  t_float panL;
  t_float panR;
  long delay; // samples to wait until event starts
  long duration;// length in samples of event
  t_float phase; // phase for frequency oscillator
  t_float ephase; // phase for envelope
  t_float si; // sampling increment for frequency
  t_float esi; // sampling increment for envelope
} t_grain;

typedef struct {
  t_word *b_samples;
  long b_frames;
  long b_nchans;
} t_pdbuffer;


typedef struct _granule
{

  t_object x_obj;
  t_float x_f;
  t_pdbuffer *wavebuf; // holds waveform samples
  t_pdbuffer *windowbuf; // holds window samples
  t_symbol *wavename; // name of waveform buffer
  t_symbol *windowname; // name of window buffer

  t_float sr; // sampling rate
  short mute;
  short hosed; // buffers are bad
  /* Global grain data*/
  long events; // number of events in a block
  long horizon; // length of block for random events
  t_float minfreq; // minimum frequency for a grain
  t_float maxfreq; // maximum frequency for a grain
  t_float minpan; // minimum pan for a grain
  t_float maxpan; // maximum pan for a grain
  t_float minamp; // minimum amplitude for a grain
  t_float maxamp; // maximum amplitude for a grain
  t_float mindur; // minimum duration for a grain
  t_float maxdur; // maximum duration for a grain
  t_grain *grains; // stores grain data
  t_float *pitchscale; // contains a frequency grid for pitch constraint
  int pitchsteps; // number of members in scale
  t_float transpose; // factor for scaling all pitches
  t_float pitch_deviation; // factor to adjust scaled pitches
  short steady; // toggles pulsed rhythmic activity
  t_float lowblock_freq; //lowest allowed frequency
  t_float highblock_freq;// highest allowed frequency
  t_float mindur_ms;//store duration in ms
  t_float maxdur_ms;//ditto
  t_float horizon_ms;//ditto
  short constrain_scale;//flag to only use bounded portion of scale rather than all of it
} t_granule;

static void granule_setbuf(t_granule *x, t_symbol *wavename, t_symbol *windowname);
static void *granule_new(t_symbol *msg, int argc, t_atom *argv);
static t_int *granule_perform(t_int *w);
static t_int *granule_performhose(t_int *w);
static void granule_dsp(t_granule *x, t_signal **sp);
static void granule_reload(t_granule *x);
static void granule_spray(t_granule *x);
static void granule_pitchspray(t_granule *x);
static void granule_transpose(t_granule *x, t_floatarg t);
static void granule_pitchdev(t_granule *x, t_floatarg d);
static void granule_lowblock(t_granule *x, t_floatarg f);
static void granule_highblock(t_granule *x, t_floatarg f);
static void granule_events(t_granule *x, t_floatarg e);
static t_float granule_boundrand(t_float min, t_float max);
static void *granule_grist(t_granule *x, t_symbol *msg, int argc, t_atom *argv);
static void *granule_grain(t_granule *x, t_symbol *msg, int argc, t_atom *argv);
static void *granule_setscale(t_granule *x, t_symbol *msg, int argc, t_atom *argv);
static void granule_info(t_granule *x);
static void granule_mute(t_granule *x, t_floatarg toggle);
static void granule_steady(t_granule *x, t_floatarg toggle);
static void granule_constrain_scale(t_granule *x, t_floatarg toggle);
static void granule_dsp_free(t_granule *x);
static void granule_init(t_granule *x,short initialized);
static void granule_constrain(int *index_min, int *index_max, t_float minfreq, t_float maxfreq, t_float *scale, int steps);

void granule_tilde_setup(void) {
  granule_class = class_new(gensym("granule~"), (t_newmethod)granule_new,
                            (t_method)granule_dsp_free,sizeof(t_granule), 0,A_GIMME,0);
  CLASS_MAINSIGNALIN(granule_class, t_granule, x_f);
  class_addmethod(granule_class,(t_method)granule_dsp, gensym("dsp"), A_CANT, 0);
  class_addmethod(granule_class,(t_method)granule_mute,gensym("mute"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_setbuf,gensym("setbuf"),A_DEFSYM,A_DEFSYM,0);
  class_addmethod(granule_class,(t_method)granule_spray,gensym("spray"),0);
  class_addmethod(granule_class,(t_method)granule_info,gensym("info"),0);
  class_addmethod(granule_class,(t_method)granule_pitchspray,gensym("pitchspray"),0);
  class_addmethod(granule_class,(t_method)granule_transpose,gensym("transpose"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_events,gensym("events"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_pitchdev,gensym("pitchdev"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_lowblock,gensym("lowblock"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_highblock,gensym("highblock"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_steady,gensym("steady"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_constrain_scale,gensym("constrain_scale"),A_FLOAT,0);
  class_addmethod(granule_class,(t_method)granule_grist,gensym("grist"),A_GIMME,0);
  class_addmethod(granule_class,(t_method)granule_grain,gensym("grain"),A_GIMME,0);
  class_addmethod(granule_class,(t_method)granule_setscale,gensym("setscale"),A_GIMME,0);
  potpourri_announce(OBJECT_NAME);
}

void granule_constrain_scale(t_granule *x, t_floatarg toggle)
{
  x->constrain_scale = toggle;
}
void granule_lowblock(t_granule *x, t_floatarg f)
{
  if(f > 0) {
    x->lowblock_freq = f;
  }
}

void granule_highblock(t_granule *x, t_floatarg f)
{
  if(f > 0) {
    x->highblock_freq = f;
  }
}

void granule_pitchdev(t_granule *x, t_floatarg d)
{
  if(d < 0 ) {
    pd_error(0, "pitch deviation must be positive");
    return;
  }
  x->pitch_deviation = d;
}

void granule_mute(t_granule *x, t_floatarg toggle)
{
  x->mute = toggle;
}

void granule_steady(t_granule *x, t_floatarg toggle)
{
  x->steady = toggle;
}

void granule_events(t_granule *x, t_floatarg e)
{
  if( e <= 0 ) {
    post("events must be positive!");
    return;
  }
  x->events = e;
  //  x->steady_dur = x->horizon / (t_float) x->events;
}

void granule_transpose(t_granule *x, t_floatarg t)
{
  if( t <= 0 ) {
    pd_error(0, "transpose factor must be greater than zero!");
    return;
  }
  x->transpose = t;
}

void *granule_setscale(t_granule *x, t_symbol *msg, int argc, t_atom *argv)
{
  int i;
  t_float *pitchscale = x->pitchscale;
  if( argc >= MAXSCALE ) {
    pd_error(0, "%d is the maximum size scale", MAXSCALE);
    return 0;
  }
  if( argc < 2 ) {
    pd_error(0, "there must be at least 2 members in scale");
    return 0;
  }
  for(i=0; i < argc; i++) {
    pitchscale[i] = atom_getfloatarg(i,argc,argv);
  }
  x->pitchsteps = argc;
  //  post("read %d values into scale", x->pitchsteps);
  return 0;
}

void granule_constrain(int *index_min, int *index_max, t_float minfreq, t_float maxfreq, t_float *scale, int steps)
{
  int imax = steps - 1;
  int imin = 0;
  while(scale[imin] < minfreq && imin < imax) {
    ++imin;
  }
  if(imin == imax) {
    //    post("could not constrain minimum index  - your grist parameters are out of range for this scale");
    *index_min = 0;
    *index_max = steps - 1;
    return;
  }
  while(scale[imax] > maxfreq && imax > 0) {
    --imax;
  }
  if(imax < 1 || imax <= imin) {
    //    post("could not constrain maximum index - your grist parameters are out of range for this scale");
    *index_min = 0;
    *index_max = steps - 1;
    return;
  }
  *index_min = imin;
  *index_max = imax;
}

void granule_pitchspray(t_granule *x)
{
  int i,j;

  long eframes = x->windowbuf->b_frames;
  long frames = x->wavebuf->b_frames;
  t_float sr = x->sr;
  long horizon = x->horizon; // length of block for random events
  t_float mindur = x->mindur;
  t_float maxdur = x->maxdur;
  t_float minfreq = x->minfreq; // minimum frequency for a grain
  t_float maxfreq = x->maxfreq; // maximum frequency for a grain
  t_float minpan = x->minpan; // minimum pan for a grain
  t_float maxpan = x->maxpan; // maximum pan for a grain
  t_float minamp = x->minamp; // minimum amplitude for a grain
  t_float maxamp = x->maxamp; // maximum amplitude for a grain
  t_float transpose = x->transpose; // pitch scalar
  t_float lowblock_freq = x->lowblock_freq;
  t_float highblock_freq = x->highblock_freq;
  short steady = x->steady;
  t_float pitch_deviation = x->pitch_deviation;
  t_float pdev = 0;
  t_float pdev_invert = 0;
  //  t_float pscale;
  t_float pan;
  int index_min, index_max;
  int steps = x->pitchsteps;
  t_float *scale = x->pitchscale;
  int windex;
  short inserted = 0;
  short constrain_scale = x->constrain_scale;
  t_grain *grains = x->grains;


  if( steps < 2 ) {
    pd_error(0, "scale is undefined");
    return;
  }
  if( pitch_deviation ) {
    pdev = 1.0 + pitch_deviation;
    pdev_invert = 1.0 / pdev;
  }
  for( i = 0; i < x->events; i++ ) {
    inserted = 0;
    for(j = 0; j < MAXGRAINS; j++ ) {
      if( grains[j].ephase >= eframes ) {
        if(steady) {
          grains[j].delay = (t_float)(i * horizon) / (t_float) x->events ;
        } else {
          grains[j].delay = granule_boundrand(0.0,(t_float) horizon);
        }
        grains[j].duration = (long) granule_boundrand(mindur, maxdur);
        grains[j].phase = 0.0;
        grains[j].ephase = 0.0;
        pan = granule_boundrand(minpan, maxpan);
        grains[j].panL = cos(pan * PIOVERTWO);
        grains[j].panR = sin(pan * PIOVERTWO);
        grains[j].amplitude = granule_boundrand(minamp, maxamp);
        grains[j].esi =  (t_float) eframes / (t_float) grains[j].duration ;
        if(constrain_scale) {
          granule_constrain(&index_min,&index_max,minfreq, maxfreq, scale, steps);
          windex = (int) granule_boundrand((t_float)index_min, (t_float)index_max);
        } else {
          windex = (int) granule_boundrand(0.0, (t_float)(steps-1));
        }
        grains[j].si = transpose * scale[windex] * (t_float) frames / sr;
        if( pitch_deviation ) {
          grains[j].si *= granule_boundrand(pdev_invert,pdev);
        }
        /* must add this code to spray, and also do for high frequencies
         */
        if(lowblock_freq > 0.0) {
          if(grains[j].si * (sr/frames) < lowblock_freq) {
            post("lowblock: aborted grain with %f frequency",grains[j].si * (sr/frames));
            grains[j].ephase = eframes; // abort grain
          }
        }
        if(highblock_freq > 0.0) {
          if(grains[j].si * (sr/frames) > highblock_freq) {
            post("highblock: aborted grain with %f frequency, greater than %f",
                 grains[j].si * (sr/frames), highblock_freq);
            grains[j].ephase = eframes; // abort grain
          }
        }
        inserted = 1;
        goto nextgrain;
      }
    }
    if(!inserted) {
      pd_error(0, "could not insert grain");
      return;
    }
  nextgrain: ;
  }
}

void granule_spray(t_granule *x)
{
  int i,j;
  long eframes = x->windowbuf->b_frames;
  long frames = x->wavebuf->b_frames;
  t_float sr = x->sr;
  long horizon = x->horizon; // length of block for random events
  t_float mindur = x->mindur;
  t_float maxdur = x->maxdur;
  t_float minfreq = x->minfreq; // minimum frequency for a grain
  t_float maxfreq = x->maxfreq; // maximum frequency for a grain
  t_float minpan = x->minpan; // minimum pan for a grain
  t_float maxpan = x->maxpan; // maximum pan for a grain
  t_float minamp = x->minamp; // minimum amplitude for a grain
  t_float maxamp = x->maxamp; // maximum amplitude for a grain
  t_float transpose = x->transpose; // pitch scalar
  //  t_float steady_dur = x->steady_dur;
  short steady = x->steady;
  t_float pan;
  t_grain *grains = x->grains;
  short inserted;

  for( i = 0; i < x->events; i++ ) {
    inserted = 0;
    for(j = 0; j < MAXGRAINS; j++ ) {
      if( grains[j].ephase >= eframes ) {
        if(steady) {
          grains[j].delay = (t_float)(i * horizon) / (t_float) x->events ;
        } else {
          grains[j].delay = granule_boundrand(0.0,(t_float) horizon);
        }
        grains[j].duration = (long) granule_boundrand(mindur, maxdur);
        grains[j].phase = 0.0;
        grains[j].ephase = 0.0;
        pan = granule_boundrand(minpan, maxpan);
        grains[j].panL = cos(pan * PIOVERTWO);
        grains[j].panR = sin(pan * PIOVERTWO);
        grains[j].amplitude = granule_boundrand(minamp, maxamp);
        grains[j].esi =  (t_float) eframes / (t_float) grains[j].duration ;
        grains[j].si = transpose * granule_boundrand(minfreq, maxfreq) * (t_float) frames / sr;
        inserted = 1;
        goto nextgrain;
      }
    }
    if(! inserted) {
      pd_error(0, "could not insert grain");
      return;
    }
  nextgrain: ;
  }
}

void *granule_grain(t_granule *x, t_symbol *msg, int argc, t_atom *argv)
{
  short inserted;
  int j;
  t_float duration, frequency, amplitude, pan;
  t_grain *grains;
  long eframes;
  long frames;
  t_float sr;

  grains = x->grains;
  eframes = x->windowbuf->b_frames;
  frames = x->wavebuf->b_frames;
  sr = x->sr;

  if(argc < 4) {
    pd_error(0, "grain takes 4 arguments, not %d",argc);
    post("duration frequency amplitude pan");
    return 0;
  }
  duration = atom_getintarg(0,argc,argv);
  frequency = atom_getfloatarg(1,argc,argv); // in ms
  amplitude = atom_getfloatarg(2,argc,argv);
  pan = atom_getfloatarg(3,argc,argv);
  if(duration <= 0.0) {
    pd_error(0, "illegal duration:%f",duration);
    return 0;
  }
  if(frequency <= 0.0) {
    pd_error(0, "illegal frequency:%f",frequency);
    return 0;
  }
  if(pan < 0.0 || pan > 1.0) {
    pd_error(0, "illegal pan:%f",pan);
    return 0;
  }
  inserted = 0;
  for(j = 0; j < MAXGRAINS; j++ ) {
    if( grains[j].ephase >= eframes ) {
      grains[j].delay = 0.0;// immediate deployment
      grains[j].duration = (long) (.001 * x->sr * duration);
      grains[j].phase = 0.0;
      grains[j].ephase = 0.0;
      grains[j].panL = cos(pan * PIOVERTWO);
      grains[j].panR = sin(pan * PIOVERTWO);
      grains[j].amplitude = amplitude;
      grains[j].esi =  (t_float) eframes / (t_float) grains[j].duration ;
      grains[j].si = frequency * (t_float) frames / sr;
      return 0;
    }
  }

  pd_error(0, "could not insert grain");
  return 0;

}

t_float granule_boundrand(t_float min, t_float max)
{
  return min + (max-min) * ((t_float) (rand() % RAND_MAX)/ (t_float) RAND_MAX);
}


void *granule_new(t_symbol *msg, int argc, t_atom *argv)
{

  t_granule *x = (t_granule *)pd_new(granule_class);
  outlet_new(&x->x_obj, gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));
  x->wavebuf = (t_pdbuffer*)getbytes(sizeof(t_pdbuffer));
  x->windowbuf = (t_pdbuffer*)getbytes(sizeof(t_pdbuffer));
  srand(time(0));

  x->pitchscale = (t_float *) getbytes(MAXSCALE * sizeof(t_float));
  x->grains = (t_grain *) getbytes(MAXGRAINS * sizeof(t_grain));


  // default names
  x->wavename = gensym("waveform");
  x->windowname = gensym("window");

  /* MaxMSP bug that may soon be fixed, this does not work:
     x->wavename = atom_getsymarg(0,argc,argv);
     x->windowname = atom_getsymarg(1,argc,argv); */

  // apparently Pd lacks this Max/MSP bug
  x->wavename = atom_getsymbolarg(0,argc,argv);
  x->windowname = atom_getsymbolarg(1,argc,argv);


  x->sr = sys_getsr();
  if(! x->sr )
    x->sr = 44100;

  granule_init(x,0);


  return (x);
}

void granule_init(t_granule *x,short initialized)
{
  int i;

  if(!initialized) {
    x->pitchsteps = 0; // we could predefine a 12t scale
    x->mute = 0;
    x->steady = 0;
    x->events = 10;
    x->horizon_ms = 1000;
    x->minfreq = 220.0;
    x->maxfreq = 880.0;
    x->minpan = .1;
    x->maxpan = .9;
    x->minamp = .1;
    x->maxamp = 1.0;
    x->mindur_ms = 150;
    x->maxdur_ms = 750;
    x->transpose = 1.0;
    x->pitch_deviation = 0.0;
    x->lowblock_freq = 0.0; // by default we do not block any frequencies
    x->highblock_freq = 0.0; // ditto
    x->constrain_scale = 0;
  }
  x->horizon = x->horizon_ms * .001 * x->sr;
  x->mindur = x->mindur_ms * .001 * x->sr;
  x->maxdur = x->maxdur_ms * .001 * x->sr;
  for( i = 0; i < MAXGRAINS; i++ ){ // this is what we test for a legal place to insert grain
    x->grains[i].ephase = 9999999999.0;
  }
}

void granule_info(t_granule *x)
{
  int tcount = 0;
  t_grain *grains = x->grains;
  long eframes = x->windowbuf->b_frames;
  int i;

  for(i = 0; i < MAXGRAINS; i++ ) {
    if( grains[i].ephase < eframes )
      ++tcount;
  }
  post("%d active grains", tcount);
  post("wavename %s", x->wavename->s_name);
  post("windowname %s", x->windowname->s_name);
}


void *granule_grist(t_granule *x, t_symbol *msg, int argc, t_atom *argv)
{
  if(argc < 10 ) {
    pd_error(0, "grist takes 10 arguments:");
    post("events horizon minfreq maxfreq minpan maxpan minamp maxamp mindur maxdur");
    return 0;
  }
  x->events = atom_getintarg(0,argc,argv);
  x->horizon_ms = atom_getfloatarg(1,argc,argv);
  x->minfreq = atom_getfloatarg(2,argc,argv);
  x->maxfreq = atom_getfloatarg(3,argc,argv);
  x->minpan = atom_getfloatarg(4,argc,argv);
  x->maxpan = atom_getfloatarg(5,argc,argv);
  x->minamp = atom_getfloatarg(6,argc,argv);
  x->maxamp = atom_getfloatarg(7,argc,argv);
  x->mindur_ms = atom_getfloatarg(8,argc,argv);
  x->maxdur_ms = atom_getfloatarg(9,argc,argv);

  x->mindur = .001 * x->sr * x->mindur_ms ;
  x->maxdur = .001 * x->sr * x->maxdur_ms;
  x->horizon = .001 * x->sr * x->horizon_ms;

  if(x->minfreq < 0) {
    x->minfreq *= -1.0;
  }
  if(x->maxfreq < 0) {
    x->maxfreq *= -1.0;
  }
  if(x->minpan < 0.0) {
    x->minpan = 0.0;
  }
  if(x->maxpan > 1.0) {
    x->maxpan = 1.0;
  }
  if(x->events < 0) {
    x->events = 0;
  }
  return 0;
}


void granule_reload(t_granule *x)
{
  granule_setbuf(x, x->wavename, x->windowname);
}


void granule_setbuf(t_granule *x, t_symbol *wavename, t_symbol *windowname)
{
  t_garray *a;
  int frames;

  x->hosed = 0;
  x->wavebuf->b_frames = 0;
  x->windowbuf->b_frames = 0;
  x->wavebuf->b_nchans = 1;
  x->windowbuf->b_nchans = 1;
  if (!(a = (t_garray *)pd_findbyclass(wavename, garray_class))) {
    if (*wavename->s_name) pd_error(x, "granule~: %s: no such array", wavename->s_name);
    x->hosed = 1;
  }
  else if (!garray_getfloatwords(a, &frames, &x->wavebuf->b_samples)) {
    pd_error(x, "%s: bad template for granule~", wavename->s_name);
    x->hosed = 1;
  }
  else  {
    x->wavebuf->b_frames = frames;
    garray_usedindsp(a);
  }

  if (!(a = (t_garray *)pd_findbyclass(windowname, garray_class))) {
    if (*wavename->s_name) pd_error(x, "granule~: %s: no such array", windowname->s_name);
    x->hosed = 1;
  }
  else if (!garray_getfloatwords(a, &frames, &x->windowbuf->b_samples)) {
    pd_error(x, "%s: bad template for granule~", windowname->s_name);
    x->hosed = 1;
  }
  else  {
    x->windowbuf->b_frames = frames;
    garray_usedindsp(a);
  }
}


t_int *granule_performhose(t_int *w)
{
  //  t_granule *x = (t_granule *) (w[1]);
  t_float *outputL = (t_float *)(w[3]);
  t_float *outputR = (t_float *)(w[4]);
  int n = (int) w[5];
  while(n--) *outputL++ = *outputR++ = 0;
  return (w+6);
}

t_int *granule_perform(t_int *w)
{
  t_granule *x = (t_granule *) (w[1]);
  //  t_float *in = (t_float *)(w[2]); // ignoring input
  t_float *outputL = (t_float *)(w[3]);
  t_float *outputR = (t_float *)(w[4]);
  int n = (int) w[5];

  t_pdbuffer *wavebuf = x->wavebuf;
  t_pdbuffer *windowbuf = x->windowbuf;
  t_word *wavetable = wavebuf->b_samples;
  t_word *window = windowbuf->b_samples;
  t_grain *grains = x->grains;
  t_float sample;
  t_float envelope;
  t_float amplitude;
  t_float panL, panR;
  t_float si;
  t_float esi;
  t_float phase;
  t_float ephase;
  long delay;
  long frames = wavebuf->b_frames;
  long eframes = windowbuf->b_frames;
  int i,j;



  /* grain parameters */


  if( x->mute ) {
    while(n--) *outputL++ = *outputR++ = 0;
    return (w+6);
  }

  // pre-clean buffer
  for( i = 0; i < n; i++ ) {
    outputL[i] = outputR[i] = 0;
  }

  for (j=0; j<MAXGRAINS; j++) {

    if(grains[j].ephase >= eframes) {
      goto nextgrain;
    }
    amplitude = grains[j].amplitude;
    si =  grains[j].si;
    esi = grains[j].esi;
    phase =  grains[j].phase;
    ephase = grains[j].ephase;
    delay =  grains[j].delay;
    panL = grains[j].panL;
    panR = grains[j].panR;


    for(i = 0; i < n; i++ ) {
      // ++(x->sampcount); // not really needed
      if( delay > 0 ) {
        --delay;
      }
      if( delay <= 0 && ephase < eframes) {
        sample = wavetable[(int)phase].w_float;

        envelope = amplitude * window[(int)ephase].w_float;
        sample *= envelope;
        outputL[i] += panL * sample;
        outputR[i] += panR * sample;
        phase += si;
        ephase += esi;
        while( phase >= frames )
          phase -= frames;

        if( ephase >= eframes ) {
          grains[j].ephase = ephase;
          goto nextgrain; // must escape loop now
        }

      }
    }
    grains[j].phase = phase;
    grains[j].ephase = ephase;
    grains[j].delay = delay;

  nextgrain: ;
  }

  return (w+6);


}

void granule_dsp_free(t_granule *x)
{

  freebytes(x->grains, MAXGRAINS * sizeof(t_grain));
  freebytes(x->pitchscale, MAXSCALE * sizeof(t_float));
}

void granule_dsp(t_granule *x, t_signal **sp)
{

  granule_reload(x);

  if( x->hosed ) {
    post("You need some valid buffers");
    dsp_add(granule_performhose, 5, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
    return;
  }
  if( x->sr != sp[0]->s_sr) {
    x->sr = sp[0]->s_sr;
    if( !x->sr ) {
      post("warning: zero sampling rate!");
      x->sr = 44100;
    }
    granule_init(x,1);
  }
  dsp_add(granule_perform, 5, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, (t_int)sp[0]->s_n);
}
