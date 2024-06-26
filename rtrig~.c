#include "MSPd.h"

static t_class *rtrig_class;

#define OBJECT_NAME "rtrig~"

#define MAXFLAMS (16)
#define MAXATTACKS (128)
#define STOPGAIN (.001)

#define RAND_MAX2 (0x7fffffff)

typedef struct
{
  int attack_count; // number of triggers per flam event
  t_float *attack_times; // trigger times in seconds
  int *attack_points; // trigger times in samples
  int fdex; // current flam
  t_float gainatten; // attenuation factor
  t_float amp; // current amp
  int atks;// number of attacks per flam
  long counter; // internal clock
  short active; // flag that flam is turned on


} t_flam;

typedef struct _rtrig
{

  t_object x_obj;
  t_float x_f;
  short mute;
  t_float min;
  t_float max;
  t_float odds;

} t_rtrig;


static void *rtrig_new(t_symbol *s, int argc, t_atom *argv);
static t_int *rtrig_perform(t_int *w);
static void rtrig_dsp(t_rtrig *x, t_signal **sp);
static void rtrig_free(t_rtrig *x);
static void rtrig_mute(t_rtrig *x, t_floatarg tog);
static void rtrig_min(t_rtrig *x, t_floatarg v);
static void rtrig_max(t_rtrig *x, t_floatarg v);
static void rtrig_odds(t_rtrig *x, t_floatarg v);


void rtrig_tilde_setup(void)
{
  rtrig_class = class_new(gensym("rtrig~"),(t_newmethod)rtrig_new,(t_method)rtrig_free, sizeof(t_rtrig), 0, A_GIMME,0);
  CLASS_MAINSIGNALIN(rtrig_class,t_rtrig, x_f );

  class_addmethod(rtrig_class,(t_method)rtrig_dsp,gensym("dsp"),A_CANT,0);
  class_addmethod(rtrig_class,(t_method)rtrig_mute,gensym("mute"),A_FLOAT,0);
  class_addmethod(rtrig_class,(t_method)rtrig_min,gensym("min"),A_FLOAT,0);
  class_addmethod(rtrig_class,(t_method)rtrig_max,gensym("max"),A_FLOAT,0);
  class_addmethod(rtrig_class,(t_method)rtrig_odds,gensym("odds"),A_FLOAT,0);
  potpourri_announce(OBJECT_NAME);
}

void rtrig_mute(t_rtrig *x, t_floatarg tog)
{
  x->mute = (short) tog;
}

void rtrig_min(t_rtrig *x, t_floatarg v)
{
  x->min = (t_float) v;
}

void rtrig_max(t_rtrig *x, t_floatarg v)
{
  x->max = (t_float) v;
}

void rtrig_odds(t_rtrig *x, t_floatarg v)
{
  x->odds = (t_float) v;
}


void *rtrig_new(t_symbol *s, int argc, t_atom *argv)
{

  t_rtrig *x = (t_rtrig *)pd_new(rtrig_class);
  outlet_new(&x->x_obj, gensym("signal"));
  x->mute = 0;
  x->min = atom_getfloatarg(0,argc,argv);
  x->max = atom_getfloatarg(1,argc,argv);
  x->odds = atom_getfloatarg(2,argc,argv);

  srand(time(0));
  return (x);
}



void rtrig_free(t_rtrig *x)
{
}

t_int *rtrig_perform(t_int *w)
{

  t_rtrig *x = (t_rtrig *) (w[1]);
  t_float *out_vec = (t_float *)(w[2]);
  int n = (int) w[3];

  t_float rval;
  t_float min = x->min;
  t_float max = x->max;
  t_float odds = x->odds;

  if(x->mute) {
    memset( (void *)out_vec, 0, n * sizeof(t_float) );
    return (w+4);
  }

  while( n-- ) {
    rval = (t_float) rand() / (t_float) RAND_MAX2;
    if(rval < odds) {
      rval = min + (max-min) * ((t_float) rand() / (t_float) RAND_MAX2);
      *out_vec++ = rval;
    } else {
      *out_vec++ = 0.0;
    }
  }


  return (w+4);
}

void rtrig_dsp(t_rtrig *x, t_signal **sp)
{
  dsp_add(rtrig_perform, 3, x, sp[1]->s_vec, (t_int)sp[0]->s_n);
}
