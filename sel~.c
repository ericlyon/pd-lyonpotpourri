#include "MSPd.h"

static t_class *sel_class;

#define MAXBEATS (256)
#define OBJECT_NAME "sel~"
#define COMPILE_DATE "9.02.07"
#define OBJECT_VERSION "2.01"
// #define DATE "prerelease"

/* Pd version of sel~ */

typedef struct _sel
{
    t_object x_obj;
    t_float x_f;
    t_float *matches; // store numbers to match against
    t_float *trigger_vec; // copy of input vector
    int length; // number of matches to check
    t_float **ins; // array of input signal vectors
    t_float **outs; // array of output signal vectors
} t_sel;

static void *sel_new(t_symbol *msg, int argc, t_atom *argv);
static void sel_free(t_sel *x);
static void sel_dsp(t_sel *x, t_signal **sp);


void sel_tilde_setup(void) {
    sel_class = class_new(gensym("sel~"), (t_newmethod)sel_new,
                          (t_method)sel_free, sizeof(t_sel),0,A_GIMME,0);
    CLASS_MAINSIGNALIN(sel_class, t_sel, x_f);
    class_addmethod(sel_class, (t_method)sel_dsp, gensym("dsp"), A_CANT, 0);
    potpourri_announce(OBJECT_NAME);
}

void *sel_new(t_symbol *msg, int argc, t_atom *argv)
{
    int i;
    
    t_sel *x = (t_sel *)pd_new(sel_class);
    x->length = (int)argc;
    
    for(i=0;i< x->length ;i++) {
        outlet_new(&x->x_obj, gensym("signal"));
    }
    
    x->matches = (t_float *) getbytes(x->length * sizeof(t_float));
    
    for(i = 0; i < argc; i++) {
        x->matches[i] = (double)atom_getfloatarg(i,argc,argv);
    }
    
    x->ins = (t_float **) getbytes(1 * sizeof(t_float *));
    x->outs = (t_float **) getbytes(x->length * sizeof(t_float *));
    // only 1 inlet
    for(i = 0; i < 1; i++) {
        x->ins[i] = (t_float *) getbytes(8192 * sizeof(t_float));
    }
    return x;
}

void sel_free(t_sel *x)
{
    freebytes(x->matches, x->length * sizeof(t_float));
    freebytes(x->outs, x->length * sizeof(t_float *));
    freebytes(x->ins[0], 8192 * sizeof(t_float));
    freebytes(x->ins, 1 * sizeof(t_float *));
}

t_int *sel_perform(t_int *w)
{
    int i, j;
    t_sel *x = (t_sel *) w[1];
    t_float **ins = x->ins;
    t_float **outs = x->outs;
    t_float *invec;
    t_float *inlet;
    t_float *match_outlet;
    t_float *matches = x->matches;
    int length = x->length;
    
    int n = (int) w[length + 3]; // obj, func, 1 inlet
    
    // copy input vectors (just 1 here)
    for(i = 0; i < 1; i++) {
        invec = (t_float *) w[2 + i];
        for(j = 0; j < n; j++) {
            ins[i][j] = invec[j];
        }
    }
    inlet = ins[0];
    // assign output vector pointers
    for(i = 0; i < length; i++) {
        outs[i] = (t_float *) w[3 + i];
    }
    
    // clean each outlet
    for(j = 0; j < length; j++) {
        match_outlet = (t_float *) outs[j];
        for(i = 0; i < n; i++) {
            match_outlet[i] = 0.0;
        }
    }
    // now match and route any clicks in the input
    for(i = 0; i < n; i++) {
        if(inlet[i]) {
            for(j = 0; j < length; j++) {
                if( inlet[i] == matches[j]) {
                    match_outlet = (t_float *) outs[j];
                    match_outlet[i] = 1.0; // always send a unity click
                }
            }
        }
    }
    return (w + length + 4);
}

/*
 void sel_perform64(t_sel *x, t_object *dsp64, double **ins,
 long numins, double **outs,long numouts, long n,
 long flags, void *userparam)
 {
 int i, j;
 t_double *inlet = ins[0];
 t_double *match_outlet;
 t_double *matches = x->matches;
 int length = x->length;
 
 // clean each outlet
 for(j = 0; j < length; j++) {
 match_outlet = (t_double *) outs[j];
 for(i = 0; i < n; i++) {
 match_outlet[i] = 0.0;
 }
 }
 // now match and route any clicks in the input
 for(i = 0; i < n; i++) {
 if(inlet[i]) {
 for(j = 0; j < length; j++) {
 if( inlet[i] == matches[j]) {
 match_outlet = (t_double *) outs[j];
 match_outlet[i] = 1.0; // always send a unity click
 }
 }
 }
 }
 }
 
 
 t_int *sel_dsp64(t_sel *x, t_object *dsp64, short *count, double sr, long n, long flags)
 {
 if(!sp[0]->s_sr)
 return;
 object_method(dsp64, gensym("dsp_add64"),x,sel_perform64,0,NULL);
 }
 */

void sel_dsp(t_sel *x, t_signal **sp)
{
    long i;
    t_int **sigvec;
    int pointer_count = x->length + 3; // output chans + object + inchan + vectorsize
    if(!sp[0]->s_sr) {
        return;
    }
    sigvec  = (t_int **) getbytes(pointer_count * sizeof(t_int *));
    for(i = 0; i < pointer_count; i++) {
        sigvec[i] = (t_int *) getbytes(sizeof(t_int) * 1);
    }
    sigvec[0] = (t_int *)x; // first pointer is to the object
    sigvec[pointer_count - 1] = (t_int *)sp[0]->s_n; // last pointer is to vector size (N)
    for(i = 1; i < pointer_count - 1; i++){ // now attach the inlet and all outlets
        sigvec[i] = (t_int *)sp[i-1]->s_vec;
    }
    dsp_addv(sel_perform, pointer_count, (t_int *)sigvec);
    freebytes(sigvec, pointer_count * sizeof(t_int *));
}
