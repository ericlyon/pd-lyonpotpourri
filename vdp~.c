
#include "MSPd.h"

#define F_LEN 16384
#define MAX_DELAY_TIME 3600000.0 // in seconds
#define OBJECT_NAME "vdp~"

static t_class *vdp_class;

typedef struct
{
  t_float coef;
  t_float cutoff;
  t_float x1;
} t_lpf;

typedef struct {
  t_word *b_samples;
  long b_valid;
  long b_nchans;
  long b_frames;
  t_symbol *wavename;
} t_guffer; // stuff we care about from garrays and buffers

typedef struct _vdp
{

  t_object x_obj;
  t_float x_f;
  t_float sr;

  t_lpf lpf;
  short filter;

  t_float speed;
  t_float feedback;
  t_float delay_time;
  t_float delay_samps;
  t_float maxdel;

  t_float *delay_line ;
  t_float *write_ptr; // location to write current input
  t_float *startmem; // first address in delay line
  t_float *endmem; // last address to read in delay line
  int len;
  int phs;
  t_float tap;
  short connections[4];

  short feedback_protect;
  short mute;
  short interpolate;
  short inf_hold;
  /* copy to buffer */
  t_guffer *destbuf; /* for copying to another buffer */
  /* tapering */
  long taper_count;
  t_float taper_feedback;
} t_vdp;

static t_int *vdp_perform(t_int *w);
static void vdp_protect(t_vdp *x, double state);
static void vdp_dsp(t_vdp *x, t_signal **sp);
static void *vdp_new(t_symbol *s, int argc, t_atom *argv);
//static void vdp_float(t_vdp *x, double f);
static void vdp_mute(t_vdp *x, t_floatarg t);
static void vdp_interpolate(t_vdp *x, t_floatarg t);
static void vdp_show(t_vdp *x);
static void vdp_coef(t_vdp *x, t_floatarg f);
static void vdp_filter(t_vdp *x, t_floatarg t);
static void vdp_init(t_vdp *x,short initialized);
static void vdp_clear(t_vdp *x);
static void vdp_inf_hold(t_vdp *x,t_floatarg state);
static void vdp_copy_to_buffer(t_vdp *x, t_symbol *msg, int argc, t_atom *argv);
static int vdp_setdestbuf(t_vdp *x, t_symbol *wavename);
static void vdp_redraw(t_vdp *x);
static void vdp_free(t_vdp *x);


void vdp_tilde_setup(void)
{
  vdp_class = class_new(gensym("vdp~"),(t_newmethod)vdp_new,(t_method)vdp_free, sizeof(t_vdp), 0, A_GIMME,0);
  CLASS_MAINSIGNALIN(vdp_class,t_vdp, x_f );

  class_addmethod(vdp_class,(t_method)vdp_dsp,gensym("dsp"),A_CANT,0);
  class_addmethod(vdp_class,(t_method)vdp_protect,gensym("protect"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_mute,gensym("mute"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_filter,gensym("filter"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_coef,gensym("coef"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_show,gensym("show"),0);
  class_addmethod(vdp_class,(t_method)vdp_clear,gensym("clear"),0);
  class_addmethod(vdp_class,(t_method)vdp_inf_hold,gensym("inf_hold"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_interpolate,gensym("interpolate"),A_FLOAT,0);
  class_addmethod(vdp_class,(t_method)vdp_copy_to_buffer,gensym("copy_to_buffer"),A_GIMME,0);
  potpourri_announce(OBJECT_NAME);
}

void vdp_mute(t_vdp *x, t_floatarg t)
{
  x->mute = (short)t;
}

void vdp_inf_hold(t_vdp *x, t_floatarg t)
{
  x->inf_hold = (short)t;
  x->taper_feedback = 1.0;
  x->taper_count = 0;
}


void vdp_filter(t_vdp *x, t_floatarg t)
{
  x->filter = (short)t;
}

void vdp_coef(t_vdp *x, t_floatarg f)
{
  x->lpf.coef = (t_float)f;
}

void vdp_show(t_vdp *x)
{
  post("feedback %f delay %f",x->feedback, x->delay_time);
}

void vdp_interpolate(t_vdp *x, t_floatarg t)
{
  x->interpolate = (short)t;
}

t_int *vdp_perform(t_int *w)
{
  // DSP config
  t_vdp *x = (t_vdp *)(w[1]);
  t_float *input = (t_float *)(w[2]);
  t_float *delay_vec = (t_float *)(w[3]);
  t_float *feedback_vec = (t_float *)(w[4]);
  t_float *output = (t_float *)(w[5]);
  int n = (int) w[6];

  t_float fdelay;
  t_float insamp;
  t_float outsamp = 0.0;
  t_float frac;
  t_float *write_ptr = x->write_ptr;
  t_float *read_ptr;
  t_float *startmem = x->startmem;
  t_float *endmem = x->endmem;
  int len = x->len;
  t_float tap = x->tap;
  t_float feedback = x->feedback;
  t_float delay_samps = x->delay_samps;
  short *connections = x->connections;
  t_float sr = x->sr;
  t_float msr = sr * 0.001;
  short feedback_protect = x->feedback_protect;
  short interpolate = x->interpolate;
  t_lpf lpf = x->lpf;
  short filter = x->filter;
  t_float x1,x2;
  int idelay;

  short inf_hold = x->inf_hold;


  /**********************/

  if( x->mute ) {
    /* while(n--) {
     *output++ = 0.0;
     } */
    memset( (char *)output, 0, n * sizeof(t_float) );
    return (w+7);
  }

  fdelay = delay_samps;
  idelay = floor(fdelay);

  /* loop only for infinite hold */

  if(inf_hold) {
    while( n-- ) {
      read_ptr = write_ptr;
      outsamp = *read_ptr;
      write_ptr++;
      if( write_ptr >= endmem ) {
        write_ptr = startmem;
      }
      *output++ = outsamp;
    }
    x->write_ptr = write_ptr;
    x->delay_samps = fdelay;
    return (w+7);
  }


  /* normal main loop*/
  while( n-- ) {

    // Pull Data off Signal buffers
    insamp = *input++;


    if ( connections[1]) {
      fdelay = *delay_vec++ * msr; // convert delay from milliseconds to samples
      if (fdelay < 0.0 )
        fdelay = 0.0;
      if( fdelay >= len )
        fdelay = len - 1;
      idelay = floor(fdelay);
    }

    if(connections[2]) {
      feedback = *feedback_vec++;
      if( feedback_protect ) {
        if( feedback > 0.99)
          feedback = 0.99;
        if( feedback < -0.99 )
          feedback = -0.99;
      }

    }


    /* make fdelay behave */
    if(fdelay < 0.0) {
      fdelay = 0.0;
    }
    else if(fdelay >= len) {
      fdelay = len - 1;
    }
    idelay = floor(fdelay);

    if(interpolate) {
      frac = (fdelay - idelay);
      read_ptr = write_ptr - idelay;
      if( read_ptr < startmem ) {
        read_ptr += len;
      }
      x1 = *read_ptr--;
      if( read_ptr < startmem ) {
        read_ptr += endmem - startmem;
      }
      x2 = *read_ptr;
      outsamp = x1 + frac * (x2 - x1);

    }
    else { // no interpolation case
      read_ptr = write_ptr - idelay;
      if( read_ptr < startmem ) {
        read_ptr += len;
      }
      outsamp = *read_ptr;

    }
    if(filter) {
      outsamp += lpf.x1 * lpf.coef;
      outsamp /= (1.0+lpf.coef);
      lpf.x1 = outsamp;
    }

    *write_ptr++ = insamp + outsamp * feedback;

    if( write_ptr >= endmem ) {
      write_ptr = startmem;
    }

    *output++ = outsamp;
  }

  x->tap = tap;
  x->feedback = feedback;
  x->delay_time = fdelay;
  x->delay_samps = fdelay;
  x->write_ptr = write_ptr;

  return (w+7);


}

void *vdp_new(t_symbol *s, int argc, t_atom *argv)
{
  int i;

  t_vdp *x = (t_vdp *)pd_new(vdp_class);
  for(i = 0; i < 2; i++) {
    inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"), gensym("signal"));
  }
  outlet_new(&x->x_obj, gensym("signal") );
  x->sr = sys_getsr();
  if(!x->sr) {
    pd_error(0, "zero sampling rate - set to 44100");
    x->sr = 44100;
  }
  // DSP CONFIG


  // SET DEFAULTS
  x->maxdel = 50.0; // milliseconds
  x->feedback = 0.5;
  x->delay_time  = 0.0;

  /*
    atom_arg_getfloat(&x->maxdel,0,argc,argv);
    atom_arg_getfloat(&x->delay_time,1,argc,argv);
    atom_arg_getfloat(&x->feedback,2,argc,argv);
  */
  x->maxdel = atom_getfloatarg(0,argc,argv);
  x->delay_time = atom_getfloatarg(1,argc,argv);
  x->feedback = atom_getfloatarg(2,argc,argv);
  x->interpolate = atom_getfloatarg(3,argc,argv);
  if(!x->maxdel)
    x->maxdel = 50.0;

  vdp_init(x,0);
  return (x);
}

void vdp_free(t_vdp *x)
{
  freebytes(x->delay_line, (x->len + 2) *  sizeof(t_float));
}

void vdp_clear(t_vdp *x)
{
  memset((char*)x->delay_line,0,(x->len + 2) * sizeof(t_float));
}


void vdp_init(t_vdp *x,short initialized)
{
  //int i;

  if(!initialized) {
    x->feedback_protect = 0;
    x->interpolate = 1;
    x->filter = 0;
    x->inf_hold = 0;
    if( x->maxdel < .00001 ) {
      x->maxdel = .00001;
    }
    if( x->maxdel > MAX_DELAY_TIME ) {
      pd_error(0, "%s: %f is too long, delay time set to max of %f",OBJECT_NAME,x->maxdel, MAX_DELAY_TIME);
      x->maxdel = MAX_DELAY_TIME;
    }
    x->len = x->maxdel * .001 * x->sr;
    x->lpf.coef = 0.5;
    x->lpf.x1 = 0.0;
    x->delay_line = (t_float *) getbytes((x->len + 2) *  sizeof(t_float));
    x->destbuf = (t_guffer *) getbytes(1 * sizeof(t_guffer));
    x->phs = 0;
    x->mute = 0;
    x->tap = 0;
  } else {
    x->len = x->maxdel * .001 * x->sr;
    x->delay_line = (t_float *) realloc(x->delay_line, (x->len + 2) * sizeof(t_float));
    memset((char*)x->delay_line,0,(x->len + 2) * sizeof(t_float));
  }
  x->startmem = x->delay_line;
  x->endmem = x->startmem + x->len;
  // x->endmem = x->startmem + (x->len - 1);

  x->write_ptr = x->startmem;
  /*
    post("startmem %d endmem %d len %d diff %d diffback %d", x->startmem, x->endmem, x->len, x->endmem - x->startmem, ( x->endmem - x->startmem) /sizeof(t_float));
  */
}

void vdp_dsp(t_vdp *x, t_signal **sp)
{
  // DSP CONFIG
  x->connections[1] = 1;
  x->connections[2] = 1;

  if(x->sr != sp[0]->s_sr) {
    x->sr = sp[0]->s_sr;
    vdp_init(x,1);
  }
  dsp_add(vdp_perform, 6, x,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec,
          (t_int)sp[0]->s_n);

}



void vdp_protect(t_vdp *x, double state)
{
  x->feedback_protect = state;
}

void vdp_copy_to_buffer(t_vdp *x, t_symbol *msg, int argc, t_atom *argv)
{
  t_symbol *destname;

  t_float *b_samples = x->delay_line;
  long b_nchans = 1;
  long b_frames = x->len + 2;
  long i;

  t_word *b_dest_samples;
  long b_dest_nchans;
  long b_dest_frames;


  destname = atom_getsymarg(0,argc,argv);

  if(! vdp_setdestbuf(x, destname)) {
    post("could not find buffer");
    return;
  }
  b_dest_samples = x->destbuf->b_samples;
  b_dest_nchans = x->destbuf->b_nchans;
  b_dest_frames = x->destbuf->b_frames;


  if(b_nchans != 1) {
    pd_error(0, "%s: buffer must be mono",OBJECT_NAME);
    return;
  }
  if(b_dest_frames < b_frames ) {
    // post("%s: destination buffer %s is too small, truncating",OBJECT_NAME,destname->s_name);
    b_frames = b_dest_frames; // local copy only
  }
  // post("cleaning out %d frames",b_dest_frames);
  /* first clean out destination */
  memset((char *)b_dest_samples, 0, b_dest_frames * 1 * sizeof(t_word));

  // post("copying %d frames",b_frames);

  /* now copy segment */
  for (i=0; i<b_frames; i++) {
    b_dest_samples[i].w_float = b_samples[i];
  }
  vdp_redraw(x);
}

int vdp_setdestbuf(t_vdp *x, t_symbol *wavename)
{

  t_garray *a;
  //  t_symbol *wavename = x->buffername;
  int b_frames;
  t_word *b_samples;
  if (!(a = (t_garray *)pd_findbyclass(wavename, garray_class))) {
    if (*wavename->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME,wavename->s_name);

    x->destbuf->b_valid = 0;
    return 0;
  }
  else if (!garray_getfloatwords(a, &b_frames, &b_samples)) {
    pd_error(x, "%s: bad array for %s", wavename->s_name,OBJECT_NAME);
    x->destbuf->b_valid = 0;
    return 0;
  }
  else  {
    x->destbuf->b_nchans = 1;
    x->destbuf->b_frames = b_frames;
    x->destbuf->b_samples = b_samples;
    x->destbuf->b_valid = 1;
    x->destbuf->wavename = wavename;
    garray_usedindsp(a);
    return(1);
  }
}

void vdp_redraw(t_vdp *x)
{
  t_garray *a;
  t_symbol *wavename = x->destbuf->wavename;
  if (!(a = (t_garray *)pd_findbyclass(wavename, garray_class))) {
    if (*wavename->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME, wavename->s_name);
    x->destbuf->b_valid = 0;
  }
  else  {
    garray_redraw(a);
  }
}
