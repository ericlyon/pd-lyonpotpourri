#include "MSPd.h"


static t_class *expflam_class;

#define OBJECT_NAME "expflam~"

#define MAXFLAMS (64)
#define MAXATTACKS (128)
#define STOPGAIN (.001)


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

typedef struct _expflam
{

  t_object x_obj;
  t_float x_f;
  t_flam *flams; // contain flams
  t_float start_delay; // initial flam delay
  t_float end_delay;// end delay
  t_float atten; // attenuation factor
  t_float slope;// slope of curve
  int atks;// number of attacks per flam
  t_float sr;
  t_float *trigvec; // hold input vector (to protect from memory sharing)
  t_float *bypvec; // ditto for flamgate vector
  short flamall; // flag to put a flam on everything
  short bypass; // flag to copy input to output without flam
  short flamgate_connected; // flag that a flamgate logical signal is connected to inlet 2

} t_expflam;


static void *expflam_new(void);
static t_int *expflam_perform(t_int *w);
static void expflam_dsp(t_expflam *x, t_signal **sp);
static void expflam_setflam(t_expflam *x, t_symbol *msg, int argc, t_atom *argv);
static void expflam_free(t_expflam *x);
// static void expflam_assist(t_expflam *x, void *b, long msg, long arg, char *dst);
static void expflam_flamall(t_expflam *x, t_floatarg tog);
static void expflam_bypass(t_expflam *x, t_floatarg tog);


void expflam_tilde_setup(void)
{
  expflam_class = class_new(gensym("expflam~"),(t_newmethod)expflam_new,
                            (t_method)expflam_free, sizeof(t_expflam), 0, 0);
  CLASS_MAINSIGNALIN(expflam_class,t_expflam, x_f );
  class_addmethod(expflam_class,(t_method)expflam_dsp,gensym("dsp"),A_CANT,0);
  class_addmethod(expflam_class,(t_method)expflam_setflam,gensym("setflam"),A_GIMME,0);
  class_addmethod(expflam_class,(t_method)expflam_flamall,gensym("flamall"),A_FLOAT,0);
  class_addmethod(expflam_class,(t_method)expflam_bypass,gensym("bypass"),A_FLOAT,0);
  potpourri_announce(OBJECT_NAME);
}

void expflam_flamall(t_expflam *x, t_floatarg tog)
{
  x->flamall = (short) tog;
}

void expflam_bypass(t_expflam *x, t_floatarg tog)
{
  x->bypass = (short) tog;
}

/*
void expflam_assist(t_expflam *x, void *b, long msg, long arg, char *dst)
{
  if (msg==1) {
    switch (arg) {
    case 0: sprintf(dst,"(signal) Trigger Click"); break;
    case 1: sprintf(dst,"(signal) Flam Gate"); break;
    }
  } else if (msg==2) {
    sprintf(dst,"(signal) Flam Clicks");
  }
}
*/

void *expflam_new(void)
{
  int i;
  t_expflam *x = (t_expflam *)pd_new(expflam_class);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"), gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));
  x->flams = (t_flam *) getbytes(MAXFLAMS * sizeof(t_flam));
  for(i = 0; i < MAXFLAMS; i++) {
    x->flams[i].attack_times = (t_float *) getbytes(MAXATTACKS * sizeof(t_float));
    x->flams[i].attack_points = (int *) getbytes(MAXATTACKS * sizeof(int));
  }

  x->trigvec = (t_float *) getbytes(8192 * sizeof(t_float)); // maximum vector size
  x->bypvec = (t_float *) getbytes(8192 * sizeof(t_float)); // maximum vector size
  x->sr = sys_getsr();
  x->start_delay = .025;
  x->end_delay = 0.1;
  x->slope = -3.0;
  x->atks = 8;
  x->atten = 0.8;
  x->bypass = 0;
  x->flamall = 0;

  return x;
}

void expflam_setflam(t_expflam *x, t_symbol *msg, int argc, t_atom *argv)
{
  if( argc != 5 ) {
    pd_error(0, "%s: setflam format: startdelay enddelay attacks slope gainatten",OBJECT_NAME);
    return;
  }
  x->start_delay = atom_getfloatarg(0,argc,argv) * 0.001;
  x->end_delay = atom_getfloatarg(1,argc,argv) * 0.001;
  x->atks = (int) atom_getfloatarg(2,argc,argv);
  x->slope = atom_getfloatarg(3,argc,argv);
  x->atten = atom_getfloatarg(4,argc,argv);
  if(x->slope == 0)
    x->slope = .0001;
  if(x->start_delay <= 0)
    x->start_delay = .00001;
  if(x->end_delay <= 0)
    x->end_delay = .00001;
  if(x->atks < 2)
    x->atks = 2;
  if(x->atks > MAXATTACKS) {
    post("%s: exceeded maximum of %d attacks",OBJECT_NAME, MAXATTACKS);
    x->atks = MAXATTACKS;
  }
}

void expflam_free(t_expflam *x)
{
  int i;

  freebytes(x->trigvec, 8192 * sizeof(t_float));
  freebytes(x->bypvec, 8192 * sizeof(t_float));
  for(i = 0; i < MAXFLAMS; i++) {
    freebytes(x->flams[i].attack_times, MAXATTACKS * sizeof(t_float));
    freebytes(x->flams[i].attack_points, MAXATTACKS * sizeof(int));
  }
  freebytes(x->flams, MAXFLAMS * sizeof(t_flam));

}

t_int *expflam_perform(t_int *w)
{
  int i,j,k;
  t_expflam *x = (t_expflam *) (w[1]);
  t_float *in_vec = (t_float *)(w[2]);
  t_float *in2_vec = (t_float *)(w[3]);
  t_float *out_vec = (t_float *)(w[4]);
  int n = (int) w[5];

  t_float *trigvec = x->trigvec;
  t_float *flamgate_vec = x->bypvec;
  t_flam *flams = x->flams;
  int atks = x->atks;
  t_float atten = x->atten;
  t_float slope = x->slope;
  t_float start_delay = x->start_delay;
  t_float end_delay = x->end_delay;
  t_float sr = x->sr;
  short flamgate_connected = x->flamgate_connected;
  short flamall = x->flamall;

  /* in flamgate mode copy input to output and return */
  if(x->bypass) {
    memcpy( (void *)out_vec, (void *)in_vec, n * sizeof(t_float) );
    return (w+6);
  }
  /* copy input vectors */
  memcpy( (void *)flamgate_vec, (void *)in2_vec, n * sizeof(t_float) );// the order of these mcopies matters
  memcpy( (void *)trigvec, (void *)in_vec, n * sizeof(t_float) );
  memcpy( (void *)out_vec, (void *)in_vec, n * sizeof(t_float) );// copy triggers to output for a start

  /* look for activation triggers */
  for(i = 0; i < n; i++) {
    if(trigvec[i] && (flamgate_vec[i] || ! flamgate_connected || flamall ) ) {
//    post("triggered with t %f and flamgate %f",trigvec[i],flamgate_vec[i]);
      j = 0;
      while(flams[j].active && j < MAXFLAMS) {
        ++j;
      }
      if(j >= MAXFLAMS) {
        post("too many flams");
      }
      else {
//        post("inserting flam at location %d",j);
        flams[j].active = 1;
        flams[j].attack_times[0] = 0.0;
        flams[j].attack_points[0] = i;
        flams[j].gainatten = atten;
        flams[j].amp = trigvec[i];
        flams[j].counter = 0;
        flams[j].fdex = 0;
        flams[j].atks = atks;

        for(k = 1; k < atks; k++) {
          flams[j].attack_times[k] = start_delay + (end_delay - start_delay) * ((1.0 - exp((t_float)k * slope/((t_float)atks-1.0)))/(1.0-exp(slope)));
          flams[j].attack_times[k] += flams[j].attack_times[k - 1];
          flams[j].attack_points[k] = flams[j].attack_times[k] * sr + i;
        }
      }
    }
  }
  /* now iterate through active flams */
  for( i = 0; i < n; i++) {
    for(j = 0; j < MAXFLAMS; j++) {
      if(flams[j].active) {
        if(flams[j].counter >= flams[j].attack_points[flams[j].fdex]) {
          out_vec[i] += flams[j].amp;
          flams[j].amp *= flams[j].gainatten;
          if( flams[j].amp <= STOPGAIN ) {
            flams[j].active = 0;
          }
          flams[j].fdex++;
          if(flams[j].fdex >= flams[j].atks) {
            flams[j].active = 0;
          }
        }
        flams[j].counter++;
      }
    }
  }

  return w+6;
}

void expflam_dsp(t_expflam *x, t_signal **sp)
{

  x->flamgate_connected = 1;
  dsp_add(expflam_perform, 5, x,
          sp[0]->s_vec,
          sp[1]->s_vec,
          sp[2]->s_vec,
          (t_int)sp[0]->s_n
    );
}
