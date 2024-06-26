#include "MSPd.h"

static t_class *markov_class;

#define OBJECT_NAME "markov~"

typedef struct _markov
{
    t_object x_obj;
    t_float x_f;
    // for markov
    int event_count;
    int maximum_length;
    t_float **event_weights;
    int current_event;
    t_float *values;
    t_float current_value;
    // for rhythm
    int count;
    int event_samples;
    int subdiv;
    t_float tempo;
    
    t_float sr;
    short manual_override;
    short trigger;
} t_markov;

//void main(void)
static void *markov_new(t_floatarg event_count);
static t_int *markov_perform(t_int *w);
static void markov_dsp(t_markov *x, t_signal **sp);
static int markov_domarkov( int current_event, t_float **event_weights, int event_count );
static void markov_subdiv(t_markov *x, t_floatarg subdiv);
static void markov_tempo(t_markov *x, t_floatarg tempo);
static void markov_set_length(t_markov *x, t_floatarg length);
static void markov_manual_override(t_markov *x, t_floatarg toggle);
static void markov_values(t_markov *x, t_symbol *msg, int argc, t_atom *argv);
static void markov_event_odds(t_markov *x, t_symbol *msg, int argc, t_atom *argv);
static void markov_free( t_markov *x);
static void markov_bang( t_markov *x);

void markov_tilde_setup(void) {
    markov_class = class_new(gensym("markov~"), (t_newmethod)markov_new,
                             (t_method)markov_free,sizeof(t_markov), 0,A_DEFFLOAT,0);
    CLASS_MAINSIGNALIN(markov_class, t_markov, x_f);
    class_addmethod(markov_class,(t_method)markov_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(markov_class,(t_method)markov_subdiv,gensym("subdiv"),A_FLOAT,0);
    class_addmethod(markov_class,(t_method)markov_tempo,gensym("tempo"),A_FLOAT,0);
    class_addmethod(markov_class,(t_method)markov_set_length,gensym("set_length"),A_FLOAT,0);
    class_addmethod(markov_class,(t_method)markov_manual_override,gensym("manual_override"),A_FLOAT,0);
    class_addmethod(markov_class,(t_method)markov_values,gensym("values"),A_GIMME,0);
    class_addmethod(markov_class,(t_method)markov_event_odds,gensym("event_odds"),A_GIMME,0);
    class_addbang(markov_class,(t_method)markov_bang);
    potpourri_announce(OBJECT_NAME);
}


void markov_free( t_markov *x)
{
    int i;
    for( i = 0; i < 10; i++ ) {
        freebytes(x->event_weights[i], x->maximum_length * sizeof(t_float));
    }
    freebytes(x->values, x->maximum_length * sizeof(t_float));
    freebytes(x->event_weights, x->maximum_length * sizeof(t_float *));
}

void markov_manual_override(t_markov *x, t_floatarg toggle)
{
    x->manual_override = toggle;
}

void markov_bang(t_markov *x)
{
    x->trigger = 1;
}

void markov_values(t_markov *x, t_symbol *msg, int argc, t_atom *argv)
{
    int i;
    
    if( argc != x->event_count ) {
        pd_error(0, "there must be %d values in this list", x->event_count);
        return;
    }
    for( i = 0; i < x->event_count ; i++) {
        x->values[i] = atom_getfloatarg(i, argc, argv);
    }
}

void markov_event_odds(t_markov *x, t_symbol *msg, int argc, t_atom *argv)
{
    int i;
    int event;
    t_float sum = 0.0;
    
    t_float **event_weights = x->event_weights;
    
    if( argc != x->event_count + 1) {
        pd_error(0, "there must be %d values in this list", x->event_count + 1);
        return;
    }
    event = atom_getfloatarg(0, argc, argv);
    if( event < 0 || event > x->event_count - 1 ) {
        pd_error(0, "attempt to set event outside range of 0 to %d",x->event_count - 1);
        return;
    }
    for( i = 0; i < x->event_count; i++) {
        event_weights[event][i] = atom_getfloatarg( (i+1), argc, argv);
        sum += event_weights[event][i];
    }
    if( sum == 0.0 ) {
        pd_error(0, "zero sum for odds - this is a very bad thing");
        return;
    } else if( sum != 1.0 ) {
        // post("sum was %f, rescaling to 1.0", sum);
        for( i = 0; i < x->event_count; i++ ) {
            event_weights[event][i] /= sum;
        }
    }
}

void markov_set_length(t_markov *x, t_floatarg length)
{
    if( length < 1 || length > x->maximum_length ) {
        pd_error(0, "%d is an illegal length", (int) length);
        return;
    }
    x->event_count = length;
    
}

void markov_tempo(t_markov *x, t_floatarg tempo)
{
    x->tempo = tempo;
    x->event_samples = x->sr * (60.0/x->tempo) / (t_float) x->subdiv;
}

void markov_subdiv(t_markov *x, t_floatarg subdiv)
{
    x->subdiv = (int) subdiv;
    if( subdiv < 1 || subdiv > 128)
        subdiv = 1;
    x->event_samples = x->sr * (60.0/x->tempo) / (t_float) x->subdiv;
}

void *markov_new(t_floatarg event_count)
{
    
    int i;
    
    t_markov *x = (t_markov *)pd_new(markov_class);
    outlet_new(&x->x_obj, gensym("signal"));
    outlet_new(&x->x_obj, gensym("signal"));
    // event_count is MAXIMUM event_count
    if( event_count < 2 || event_count > 256 ) {
        pd_error(0, "maximum event length limited to 256, set to 16 here");
        event_count = 16 ;
    }
    x->maximum_length = event_count;
    
    x->event_count = 4; // default pattern
    x->count = 0;
    
    x->event_weights = (t_float **) getbytes(x->maximum_length * sizeof(t_float *));
    for( i = 0; i < 10; i++ ) {
        x->event_weights[i] = (t_float *) getbytes(x->maximum_length * sizeof(t_float));
    }
    x->values = (t_float *) getbytes(x->maximum_length * sizeof(t_float));
    
    x->current_event = 0;
    x->values[0] = 300;
    x->values[1] = 400;
    x->values[2] = 500;
    x->values[3] = 600;
    // weights
    x->event_weights[0][0] = 0;
    x->event_weights[0][1] = 0.5;
    x->event_weights[0][2] = 0.5;
    x->event_weights[0][3] = 0;
    x->event_weights[1][0] = 0.25;
    x->event_weights[1][1] = 0.5;
    x->event_weights[1][2] = 0.0;
    x->event_weights[1][3] = 0.25;
    x->event_weights[2][0] = 1;
    x->event_weights[2][1] = 0.0;
    x->event_weights[2][2] = 0.0;
    x->event_weights[2][3] = 0.0;
    x->event_weights[3][0] = 0.33;
    x->event_weights[3][1] = 0.33;
    x->event_weights[3][2] = 0.34;
    x->event_weights[3][3] = 0.0;
    x->current_value = x->values[ x->current_event ];
    x->count = 0;
    x->tempo = 60.0;
    x->sr = sys_getsr();
    if( ! x->sr ) {
        pd_error(0, "zero sampling rate - set to 44100");
        x->sr = 44100;
    }
    x->subdiv = 1;
    x->event_samples = x->sr * (60.0/x->tempo) / (t_float) x->subdiv;
    
    x->trigger = 0;
    
    return x;
}


t_int *markov_perform(t_int *w)
{
    
    t_markov *x = (t_markov *) (w[1]);
    t_float *out = (t_float *)(w[3]);
    t_float *sync = (t_float *)(w[4]);
    int n = (int)(w[5]);
    
    int count = x->count;
    int event_samples = x->event_samples;
    int event_count = x->event_count;
    t_float **event_weights = x->event_weights;
    int current_event = x->current_event;
    t_float *values = x->values;
    t_float current_value = x->current_value;
    
    
    if( x->manual_override ) {
        while (n--) {
            if( x->trigger ) {
                current_event = markov_domarkov( current_event, event_weights, event_count );
                current_value = values[ current_event ];
                x->trigger = 0;
            }
            *out++ = current_value;
        }
        x->current_value = current_value;
        x->current_event = current_event;
        return (w+6);
    }
    
    while (n--) {
        
        if( ++count >= event_samples ) {
            current_event = markov_domarkov( current_event, event_weights, event_count );
            current_value = values[ current_event ];
            count = 0;
        }
        
        *sync++ = (t_float) count / (t_float) event_samples;
        
        *out++ = current_value;
        
    }
    
    x->current_value = current_value;
    x->count = count;
    x->current_event = current_event;
    
    return (w+6);
}

int markov_domarkov(int current_event, t_float **event_weights, int event_count)
{
    t_float randval;
    int i;
    
    randval =  rand() % 32768 ;
    randval /= 32768.0;
    
    
    for( i = 0; i < event_count; i++ ) {
        if( randval < event_weights[current_event][i] ) {
            return i;
        }
        randval -= event_weights[current_event][i];
    }
    return 0; // should never happen
}

void markov_dsp(t_markov *x, t_signal **sp)
{
    //  long i;
    if(x->sr!=sp[0]->s_sr) {
        x->sr=sp[0]->s_sr;// BUG!!! in MSP code was !=
        x->event_samples = x->sr * (60.0/x->tempo) / (t_float) x->subdiv;
        x->count = 0;
    }
    dsp_add(markov_perform, 5, x, sp[0]->s_vec , sp[1]->s_vec, sp[2]->s_vec,
            (t_int)sp[0]->s_n);
}
