#include "MSPd.h"
#define OBJECT_NAME "npan~"

#define SPEAKERMAX (1024)

static t_class *npan_class;

typedef struct _npan
{
  t_object x_obj;
  double pi_over_two;
  double twopi;
  t_float x_f;
  int outcount;
  t_float *input_locvec;
  t_float *panner_locvec;
} t_npan;

static void *npan_new(t_symbol *s, int argc, t_atom *argv);
static t_int *npan_perform(t_int *w);
static void npan_dsp(t_npan *x, t_signal **sp);
static void npan_free(t_npan *x);

void npan_tilde_setup(void) {
  npan_class = class_new(gensym("npan~"), (t_newmethod)npan_new,
                         (t_method)npan_free, sizeof(t_npan),0,A_GIMME,0);
  CLASS_MAINSIGNALIN(npan_class, t_npan, x_f);
  class_addmethod(npan_class, (t_method)npan_dsp, gensym("dsp"), A_CANT, 0);
  potpourri_announce(OBJECT_NAME);
}

void *npan_new(t_symbol *s, int argc, t_atom *argv)
{
  t_npan *x;
  x = (t_npan *)pd_new(npan_class);
  int i;
  x->outcount = (int) atom_getfloatarg(0, argc, argv);
  if( x->outcount < 2 || x->outcount > SPEAKERMAX ) {

    pd_error(0, "npan~: output count %d exceeded range limits of 2 to %d",x->outcount, SPEAKERMAX);
    x->outcount = SPEAKERMAX;
  }
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("signal"),gensym("signal")); // create 2nd inlet
  for(i = 0; i < x->outcount; i++) {
    outlet_new(&x->x_obj, gensym("signal"));
  }
  x->input_locvec = (t_float *) getbytes(8192 * sizeof(t_float));
  x->panner_locvec = (t_float *) getbytes(8192 * sizeof(t_float));
  x->pi_over_two = 1.5707963267948965;
  x->twopi = 6.283185307179586;
  return x;
}

void npan_free(t_npan *x)
{
  freebytes(x->panner_locvec, 8192 * sizeof(t_float));
  freebytes(x->input_locvec, 8192 * sizeof(t_float));
}

t_int *npan_perform(t_int *w)
{
  t_npan *x = (t_npan*) w[1];
  int outcount = x->outcount;
  t_float *input = (t_float*) w[2];
  t_float *panner = (t_float*) w[3];
  t_float *outlet1, *outlet2, *cleanoutlet;
  t_float *input_locvec = x->input_locvec;
  t_float *panner_locvec = x->panner_locvec;

  double gain1, gain2;
  double insamp;
  int chan1, chan2;
  double panloc, frak;
  int i,j;

  double pi_over_two = x->pi_over_two;
  int n = (int) w[4 + outcount];
  // copy inputs
  for(i = 0; i < n; i++) {
    input_locvec[i] = input[i];
    panner_locvec[i] = panner[i];
  }
  // clean all outlets
  for( i = 0; i < outcount; i++ ) {
    cleanoutlet = (t_float*) w[4 + i];
    for(j = 0; j < n; j++) {
      cleanoutlet[j] = 0.0;
    }
  }
  for(i = 0; i < n; i++) {
    insamp = input_locvec[i];
    panloc = panner_locvec[i];
    if( panloc < 0 )
      panloc = 0.0;
    if( panloc >= 1 ) // wrap around (otherwise crash on outlet out of range)
      panloc = 0.0;
    panloc *= (double) outcount;
    chan1 = floor( panloc );
    chan2 = (chan1 + 1) % outcount;
    frak = ( panloc - chan1 ) * pi_over_two;
    gain1 = cos( frak );
    gain2 = sin( frak );
    outlet1 = (t_float*) w[chan1 + 4]; // add offset
    outlet2 = (t_float*) w[chan2 + 4];

    outlet1[i] = insamp * gain1;
    outlet2[i] = insamp * gain2;
  }
  return (w + outcount + 5);
}

void npan_dsp(t_npan *x, t_signal **sp)
{
  long i;
  t_int **sigvec;
  int pointer_count = x->outcount + 4;
  sigvec  = (t_int **) getbytes(pointer_count * sizeof(t_int *));
  for(i = 0; i < pointer_count; i++) {
    sigvec[i] = (t_int *) getbytes(sizeof(t_int) * 1);
  }
  sigvec[0] = (t_int *)x; // first pointer is to the object
  sigvec[pointer_count - 1] = (t_int *)sp[0]->s_n; // last pointer is to vector size (N)
  for(i = 1; i < pointer_count - 1; i++){ // now attach the inlet and all outlets
    sigvec[i] = (t_int *)sp[i-1]->s_vec;
  }
  dsp_addv(npan_perform, pointer_count, (t_int *) sigvec);
  freebytes(sigvec,sizeof(t_int) * 1);
}
