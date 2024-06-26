#include "MSPd.h"
#define FUNC_LEN (32768)
#define OBJECT_NAME "phasemod~"


static t_class *phasemod_class;


typedef struct _phasemod
{

  t_object x_obj;
  t_float x_f;
  t_float x_val;
  t_float mygain;
  t_float *wavetab;
  t_float phs;
  t_float bendphs;
  t_float frequency;
  t_float alpha;
  short mute;
  short connections[4];
  t_float si_fac;
  t_float sr;
} t_phasemod;

static void *phasemod_new(t_symbol *s, int argc, t_atom *argv);
static t_int *phasemod_perform(t_int *w);
static void phasemod_mute(t_phasemod *x, t_floatarg toggle);
static void phasemod_dsp(t_phasemod *x, t_signal **sp);
static void phasemod_dsp_free(t_phasemod *x);

void phasemod_tilde_setup(void) {
  phasemod_class = class_new(gensym("phasemod~"), (t_newmethod)phasemod_new,
                             (t_method)phasemod_dsp_free,sizeof(t_phasemod), 0,A_GIMME,0);
  CLASS_MAINSIGNALIN(phasemod_class, t_phasemod, x_f);
  class_addmethod(phasemod_class,(t_method)phasemod_dsp, gensym("dsp"), A_CANT, 0);
  class_addmethod(phasemod_class,(t_method)phasemod_mute,gensym("mute"),A_FLOAT,0);
  potpourri_announce(OBJECT_NAME);
}

void phasemod_dsp_free( t_phasemod *x )
{
  freebytes(x->wavetab, FUNC_LEN * sizeof(t_float));
}

void phasemod_mute(t_phasemod *x, t_floatarg toggle)
{
  x->mute = toggle;
}
void phasemod_assist (t_phasemod *x, void *b, long msg, long arg, char *dst)
{
  if (msg==1) {
    switch (arg) {
    case 0:
      sprintf(dst,"(signal/t_float) Frequency ");
      break;
    case 1:
      sprintf(dst,"(signal/t_float) Slope Factor ");
      break;
    }
  } else if (msg==2) {
    sprintf(dst,"(signal) Output ");
  }
}

void *phasemod_new(t_symbol *s, int argc, t_atom *argv)
{
  int i;
  t_phasemod *x = (t_phasemod *)pd_new(phasemod_class);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"), gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));
  x->phs = 0;
  x->mute = 0;
  x->frequency = 440.0;

  x->wavetab = (t_float *) getbytes(FUNC_LEN * sizeof(t_float));
  for( i = 0 ; i < FUNC_LEN; i++ ) {
    x->wavetab[i] = sin( TWOPI * ((t_float)i/(t_float) FUNC_LEN)) ;
  }
  x->bendphs = 0;
  x->sr = sys_getsr();
  if(!x->sr)
    x->sr = 44100.0;
  x->si_fac = 1.0/x->sr;
  return (x);
}

t_int *phasemod_perform(t_int *w)
{

  t_float phs;

  t_phasemod *x = (t_phasemod *) (w[1]);
  t_float *frequency_vec = (t_float *)(w[2]);
  t_float *alpha_vec = (t_float *)(w[3]);
  t_float *out = (t_float *)(w[4]);
  int n = (int) w[5];

  short *connections = x->connections;
  t_float bendphs = x->bendphs;
  t_float *wavetab = x->wavetab;
  t_float si_fac = x->si_fac;

  t_float incr = x->frequency * si_fac ;
  t_float alpha = x->alpha;

  if( x->mute ) {
    while(n--) {
      *out++ = 0.0;
    }
    return (w + 6);
  }

  while (n--) {
    if( connections[1] ) {
      alpha = *alpha_vec++;
    }
    if( alpha == 0 ) {
      alpha = .000001;
    }

    if( connections[0] ) {
      incr = *frequency_vec++ * si_fac ;
    }
    // NO NEGATIVE FREQUENCIES
    if( incr < 0 )
      incr = -incr;

    bendphs += incr ;
    while( bendphs > 1.0 )
      bendphs -= 1.0 ;
    phs =   FUNC_LEN * ( (1 - exp(bendphs * alpha))/(1 - exp(alpha))  );

    while( phs < 0.0 ) {
      phs += FUNC_LEN;
    }
    while( phs >= FUNC_LEN ) {
      phs -= FUNC_LEN;
    }
    *out++ =  wavetab[(int) phs] ;
  }

  x->bendphs = bendphs;
  return (w+6);
}

void phasemod_dsp(t_phasemod *x, t_signal **sp)
{

  x->connections[0] = 1;
  x->connections[1] = 1;

  if(x->sr != sp[0]->s_sr) {
    if(!sp[0]->s_sr) {
      pd_error(0, "zero sampling rate");
      return;
    }
    x->sr = sp[0]->s_sr;
    x->si_fac = 1.0/x->sr;
  }
  dsp_add(phasemod_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,  (t_int)sp[0]->s_n);
}
