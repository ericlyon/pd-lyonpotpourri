#include "MSPd.h"

#define OBJECT_NAME "waveshape~"

#define ws_MAXHARMS (256)


static t_class *waveshape_class;

typedef struct _waveshape
{
    
    t_object x_obj;
    t_float x_f;
    int flen;
    t_float *wavetab;
    t_float *tempeh; // work function
    int hcount;
    t_float *harms;
    short mute;
} t_waveshape;

static void *waveshape_new(void);
static t_int *waveshape_perform(t_int *w);
static void waveshape_dsp(t_waveshape *x, t_signal **sp);
static void waveshape_list (t_waveshape *x, t_symbol *msg, int argc, t_atom *argv);
static void update_waveshape_function( t_waveshape *x );
//t_float mapp();
static void waveshape_mute(t_waveshape *x, t_floatarg tog);
static void waveshape_free(t_waveshape *x);


void waveshape_tilde_setup(void) {
    waveshape_class = class_new(gensym("waveshape~"), (t_newmethod)waveshape_new,
                                (t_method)waveshape_free,sizeof(t_waveshape), 0,0);
    CLASS_MAINSIGNALIN(waveshape_class, t_waveshape, x_f);
    class_addmethod(waveshape_class,(t_method)waveshape_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(waveshape_class,(t_method)waveshape_mute,gensym("mute"),A_FLOAT,0);
    class_addmethod(waveshape_class,(t_method)waveshape_list,gensym("list"),A_GIMME,0);
    potpourri_announce(OBJECT_NAME);
}


void waveshape_free(t_waveshape *x)
{
    freebytes(x->wavetab, x->flen * sizeof(t_float));
    freebytes(x->tempeh, x->flen * sizeof(t_float));
    freebytes(x->harms, ws_MAXHARMS * sizeof(t_float));
}


void waveshape_list (t_waveshape *x, t_symbol *msg, int argc, t_atom *argv)
{
    short i;
 //   t_symbol *fraud;
 //   fraud = msg;
    x->hcount = 0;
    for (i=0; i < argc; i++) {
        if (argv[i].a_type == A_FLOAT) {
            x->harms[ x->hcount ] = argv[i].a_w.w_float;
            ++(x->hcount);
        }
    }
    update_waveshape_function( x );
    
}

void waveshape_mute(t_waveshape *x, t_floatarg tog)
{
    x->mute = tog;
}

void *waveshape_new(void)
{
    t_waveshape *x = (t_waveshape *)pd_new(waveshape_class);
    long i;
    outlet_new(&x->x_obj, gensym("signal"));
    
    x->flen = 65536;
    x->wavetab = (t_float *) getbytes(x->flen * sizeof(t_float));
    x->tempeh = (t_float *) getbytes(x->flen * sizeof(t_float));
    x->harms = (t_float *) getbytes(ws_MAXHARMS * sizeof(t_float));
    
    x->hcount = 4;
    x->harms[0] = 0;
    x->harms[1] = .33;
    x->harms[2] = .33;
    x->harms[3] = .33;
    x->mute = 0;
    update_waveshape_function(x);
    return x;
}

void update_waveshape_function( t_waveshape *x ) {
    t_float point;
    int i, j;
    t_float min, max;
    for( i = 0; i < x->flen; i++ ) {
        x->tempeh[i] = 0;
    }
    for( i = 0 ; i < x->hcount; i++ ) {
        if( x->harms[i] > 0.0 ) {
            for( j = 0; j < x->flen; j++ ) {
                point = -1.0 + 2.0 * ( (t_float) j / (t_float) x->flen) ;
                x->tempeh[j] += x->harms[i] * cos( (t_float) i * acos( point ) );
            }
        }
    }
    min = 1; max = -1;
    for( j = 0; j < x->flen; j++ ) {
        if( min > x->tempeh[j] )
            min = x->tempeh[j];
        if( max < x->tempeh[j] )
            max = x->tempeh[j];
    }
    //  post("min:%f, max:%f",min,max);
    // normalize from -1 to +1
    if( (max - min) == 0 ) {
        post("all zero function - watch out!");
        return;
    }
    for( j = 0; j < x->flen; j++ ) {
        x->tempeh[j] = -1.0 + ( (x->tempeh[j] - min) / (max - min) ) * 2.0 ;
    }
    // put tempeh into waveshape function
    for( j = 0; j < x->flen; j++ ) {
        x->wavetab[j] = x->tempeh[j];
    }
}

t_int *waveshape_perform(t_int *w)
{
    t_float insamp; // , waveshape, ingain ;
    int windex ;
    
    t_waveshape *x = (t_waveshape *) (w[1]);
    t_float *in = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int) w[4];
    int flenm1 = x->flen - 1;
    t_float *wavetab = x->wavetab;
    
    if(x->mute) {
        while(n--) {
            *out++ = 0.0;
        }
        return w+5;
    }
    
    while (n--) {
        insamp = *in++;
        if(insamp > 1.0) {
            insamp = 1.0;
        }
        else if(insamp < -1.0) {
            insamp = -1.0;
        }
        windex = ((insamp + 1.0)/2.0) * (t_float)flenm1 ;
        *out++ = wavetab[windex] ;
    }
    
    return (w+5);
}

void waveshape_dsp(t_waveshape *x, t_signal **sp)
{
    dsp_add(waveshape_perform, 4, x, sp[0]->s_vec,sp[1]->s_vec,(t_int)sp[0]->s_n);
}
