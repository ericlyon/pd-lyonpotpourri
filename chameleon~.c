#import "MSPd.h"
#include "chameleon_pd.h"
#include "stdlib.h"
/* Pure Data version of chameleon */

/* maximum number of parameters for each pattern */
#define MAX_PARAMETERS 1024
/* maximum number of patterns that a chameleon~ can store */
#define MAX_SLOTS 512
/* number of currently coded DSP processes */
#define PROCESS_COUNT 20
/* turn on for debug print messages*/
#define DEBUG_CHAMELEON (0)

#define OBJECT_NAME "chameleon~"

static t_class *chameleon_class;


static void *chameleon_new(t_symbol *msg, short argc, t_atom *argv);
t_int *chameleon_perform_hosed(t_int *w);
//static void chameleon_dsp(t_chameleon *x, t_signal **sp, short *count);
static void chameleon_dsp(t_chameleon *x, t_signal **sp);
static void chameleon_assist (t_chameleon *x, void *b, long msg, long arg, char *dst);
static void chameleon_dsp_free(t_chameleon *x);
static void chameleon_set_parameters(t_chameleon *x);
static t_int *chameleon_perform(t_int *w);
static void chameleon_copy_to_MSP_buffer(t_chameleon *x, int slot);
static void chameleon_clear_presets(t_chameleon *x);

/*user messages*/
static void chameleon_info(t_chameleon *x);
static void chameleon_dblclick(t_chameleon *x);
static void chameleon_mute(t_chameleon *x, t_floatarg t);
static void chameleon_maximum_process(t_chameleon *x, t_floatarg n);
static void chameleon_minimum_process(t_chameleon *x, t_floatarg n);
static void chameleon_setbuf(t_chameleon *x, t_symbol *wavename);
static void attach_buffer(t_chameleon *x);
static void chameleon_flatodds(t_chameleon *x);
static void chameleon_randodds(t_chameleon *x);
static void chameleon_killproc(t_chameleon *x, t_floatarg p);
static void chameleon_soloproc(t_chameleon *x, t_floatarg p);
static void chameleon_verbose(t_chameleon *x, t_floatarg t);
static void chameleon_block_dsp(t_chameleon *x, t_floatarg t);
static void chameleon_gozero(t_chameleon *x);
static void chameleon_grab(t_chameleon *x);
static void chameleon_setodds(t_chameleon *x,t_symbol *msg, short argc, t_atom *argv);
static void chameleon_tcycle(t_chameleon *x,t_symbol *msg, short argc, t_atom *argv);
static void chameleon_version(t_chameleon *x);
static void chameleon_dsp64(t_chameleon *x, t_object *dsp64, short *count,
                     double samplerate, long maxvectorsize, long flags);
static void chameleon_store(t_chameleon *x, t_floatarg t);
static void chameleon_recall(t_chameleon *x, t_floatarg t);
static void chameleon_recall_parameters_exec(t_chameleon *x);
static void chameleon_loadslot(t_chameleon *x,t_symbol *msg, short argc, t_atom *argv);
static void chameleon_memory(t_chameleon *x);
;
/* CMIX function code */

static void killdc( double *inbuf, int in_frames, int channels, t_chameleon *x);
static void ringmod(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void ringmod4(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void retrograde(t_chameleon *x, int slot, int *pcount);
static void comber(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void flange(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void butterme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void sweepreson(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void truncateme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void slidecomb(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void reverb1(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void ellipseme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void feed1me(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void flam1(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void flam2(t_chameleon *x, int slot, int *pcount);
static void comb4(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void ringfeed(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void resonchameleon(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void stv(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void compdist(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void bendy(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void slideflam(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void bitcrush(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void chameleon_set_parameters_exec(t_chameleon *x);
static void chameleon_print_parameters(t_chameleon *x);
static void chameleon_report(t_chameleon *x);
/* static void chameleon_perform64(t_chameleon *x, t_object *dsp64, double **ins,
                         long numins, double **outs,long numouts, long vectorsize,
                         long flags, void *userparam); */
static t_int *chameleon_perform(t_int *w);
static void chameleon_tweak_parameters(t_chameleon *x, t_floatarg t);
static void clip(double *x, double min, double max);

/* Chameleon DSP utility prototypes - moved to local file here */
/*function prototypes*/

static void putsine (double *arr, long len);
static float boundrand(float min, float max);
static void mycombset(double loopt,double rvt,int init,double *a,double srate);
static double mycomb(double samp,double *a);
static void setweights(float *a, int len);
static void  delset2(double *a,int *l,double xmax, double srate);
static void  delset3(double *a,int *l,double xmax,double srate, double alloc_max);
static void delput2(double x,double *a,int *l);
static double dliget2(double *a,double dwait,int *l,double srate);
static void lobut(double *a, double cutoff,double SR);
static void hibut(double *a, double cutoff, double SR);
static void bpbut(double *a, double formant, double bandwidth, double SR);
static void butter_filter(double *buf, double *a, long frames);
void rsnset2(double cf,double bw,double scl,double xinit,double *a,double srate);
static double reson(double x,double *a);
static void resonadsr(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2);
static void ellipset(double *list, LSTRUCT *eel, int  *nsects, double *xnorm);
static double ellipse(double x, LSTRUCT *eel, int nsects, double xnorm);
static double allpass(double samp,double *a);
static void init_reverb_data(double *a);
static void init_ellipse_data(double **a);
static void killdc( double *inbuf, int in_frames, int channels, t_chameleon *x);
// static void setExpFlamFunc(float *arr, int flen, float v1,float v2,float alpha);
static void setflamfunc1(double *arr, int flen);
static void funcgen1(double *outArray, int outlen, double duration, double outMin, double outMax,
              double speed1, double speed2, double gain1, double gain2, double *phs1, double *phs2,
              double *sine, int sinelen);
static void normtab(double *inarr,double *outarr, double min, double max, int len);
static double mapp(double insamp,double imin,double imax,double omin,double omax);
static double oscil(double amp,double si,double *farray,int len,double *phs);
static void set_dcflt(double *a);
static void set_distortion_table(double *arr, double cut, double max, int len);
static double dlookup(double samp,double *arr,int len);
static void do_compdist(double *in,double *out,int sampFrames,int nchans,int channel,
                 double cutoff,double maxmult,int lookupflag,double *table,int range,double bufMaxamp);
static double getmaxamp(double *arr, int len);
static void buildadsr(CMIXADSR *a);
/*chameleon dsp functions */
static void feed1(double *inbuf, double *outbuf, int in_frames, int out_frames,int channels, double *functab1,
           double *functab2,double *functab3,double *functab4,int funclen,
           double duration, double maxDelay, t_chameleon *x);
static void reverb1me(double *in, double *out, int inFrames, int out_frames, int nchans, int channel, double revtime, double dry, t_chameleon *x);
// static void atom_arg_getfloat(float *c, long idx, long ac, t_atom *av);

/* Main Pure Data code */

void chameleon_tilde_setup(void)
{
    chameleon_class = class_new(gensym("chameleon~"),(t_newmethod)chameleon_new,(t_method)chameleon_dsp_free, sizeof(t_chameleon), 0, A_GIMME,0);
    CLASS_MAINSIGNALIN(chameleon_class,t_chameleon, x_f);
    class_addmethod(chameleon_class,(t_method)chameleon_dsp,gensym("dsp"),A_CANT,0);
    class_addmethod(chameleon_class,(t_method)chameleon_flatodds,gensym("flatodds"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_randodds,gensym("randodds"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_memory,gensym("memory"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_soloproc,gensym("soloproc"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_killproc,gensym("killproc"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_store,gensym("store"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_recall,gensym("recall"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_setodds,gensym("setodds"), A_GIMME, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_loadslot,gensym("loadslot"), A_GIMME, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_set_parameters,gensym("set_parameters"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_tweak_parameters,gensym("tweak_parameters"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_maximum_process,gensym("maximum_process"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_minimum_process,gensym("minimum_process"), A_FLOAT, 0);
    class_addmethod(chameleon_class,(t_method)chameleon_print_parameters,gensym("print_parameters"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_report,gensym("report"), 0);
    class_addmethod(chameleon_class,(t_method)chameleon_clear_presets,gensym("clear_presets"), 0);
    potpourri_announce(OBJECT_NAME);
}


void chameleon_print_parameters(t_chameleon *x){
    int i;
    post("loadslot 9999 %d",x->pcount);
    for(i = 0; i < x->pcount; i++){
        post("%f",x->params[i]);
    }
}

void chameleon_clear_presets(t_chameleon *x){
    int i;
    for(i = 0; i < MAX_SLOTS; i++){
        x->slots[i].pcount = 0;
    }
    x->stored_slot_count = 0;
}

void chameleon_report(t_chameleon *x){
    t_atom *data = x->data;
    t_symbol *loadslot = gensym("loadslot");
    t_symbol *comma = gensym(",");
    t_slot *slots = x->slots;
    
  //  t_atom atom_comma;
    long stored_slot_count = x->stored_slot_count;
    long slots_printed = 0;
    int data_index = 0;
    int i,j;
    

    for(i = 0; i < MAX_SLOTS; i++){
        if( slots[i].pcount > 0 ){
            SETSYMBOL(data+data_index, loadslot); data_index++;
            SETFLOAT(data+data_index, (float)i); data_index++;
            SETFLOAT(data+data_index, (float)slots[i].pcount); data_index++;
            for(j = 0; j < slots[i].pcount; j++){
                SETFLOAT(data+data_index, slots[i].params[j]); data_index++;
            }
            if( slots_printed < (stored_slot_count - 1)){
                (data+data_index)->a_type = A_COMMA;
                (data+data_index)->a_w.w_index = 0;
                data_index++;
            }
            slots_printed++;
        }
    }
    outlet_list(x->listo, 0, data_index, data);
}

void chameleon_maximum_process(t_chameleon *x, t_floatarg n)
{
    if(n < 0){
        error("chameleon~: illegal val to maximum_process");
        return;
    }
    x->max_process_per_note = (int)n;
}

void chameleon_minimum_process(t_chameleon *x, t_floatarg n)
{
    if(n < 0){
        error("chameleon~: illegal val to minimum_process");
        return;
    }
    x->min_process_per_note = (int)n;
}

void chameleon_setodds(t_chameleon *x,t_symbol *msg, short argc, t_atom *argv)
{
    int i;
    
    if(argc > PROCESS_COUNT){
        error("chameleon~: there are only %d processes",PROCESS_COUNT);
        return;
    }
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = 0.0;
    }
    for(i=0;i<argc;i++){
        x->odds[i] = atom_getfloatarg(i,argc,argv);
    }
    setweights(x->odds,PROCESS_COUNT);
}

void chameleon_soloproc(t_chameleon *x, t_floatarg fp)
{
    int i;
    int p = (int) fp;
    if(p < 0 || p >= PROCESS_COUNT){
        error("chameleon~: bad %d",p);
    }
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = 0.0;
    }
    x->odds[p] = 1.0;
    setweights(x->odds,PROCESS_COUNT);
}


void chameleon_store(t_chameleon *x, t_floatarg fp)
{
    long slotnum = (long) fp;
    int i;
    if( (slotnum < 0) || (slotnum >= MAX_SLOTS)){
        pd_error(x, "%ld is not a valid slot number", slotnum);
        return;
    }
    // we're good, store the data
    // post("storing %d pieces of data at slot %d", x->pcount, slotnum);
    x->slots[slotnum].pcount = x->pcount;
    for(i = 0; i < x->pcount; i++){
        x->slots[slotnum].params[i] = x->params[i];
    }
    x->stored_slot_count += 1;
}

void chameleon_loadslot(t_chameleon *x,t_symbol *msg, short argc, t_atom *argv)
{
    int i;
    t_float f_arg;
    int slot, pcount;
    
    atom_arg_getfloat(&f_arg,0,argc, argv);
    slot = (int) f_arg;
    atom_arg_getfloat(&f_arg,1,argc, argv);
    pcount = (int) f_arg;
    
   // post("args are pcount: %d and slot: %d", pcount, slot);
    post("chameleon~: loaded slot %d", slot);
    if(argc < pcount + 2){
        pd_error(x, "wrong number of arguments to loadslot. Should be %d, got %d", pcount, argc);
    }
    for(i = 0; i < pcount; i++){
       // atom_arg_getfloat( &x->slots[slot].params[i], 2 + i, argc, argv);
        x->slots[slot].params[i] = atom_getfloatarg(2 + i, argc, argv);
        // post("data %d: %f", i, x->slots[slot].params[i]);
    }
    x->slots[slot].pcount = pcount;
    x->stored_slot_count += 1;
}


void chameleon_recall(t_chameleon *x, t_floatarg fp)
{
    long slotnum = (long) fp;
    if( (slotnum < 0) || (slotnum >= MAX_SLOTS)){
        pd_error(x, "%ld is not a valid slot number", slotnum);
        return;
    }
    // post("preparing to recall slot %d", slotnum);
    x->recall_slot = slotnum;
    x->recall_parameters_flag = 1;
}

void chameleon_killproc(t_chameleon *x, t_floatarg fp)
{
    int i;
    int p = (int) fp;
    if(p < 0 || p >= PROCESS_COUNT){
        error("chameleon~: bad %d",p);
    }
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = 1.0;
    }
    x->odds[p] = 0.0;
    setweights(x->odds,PROCESS_COUNT);
}

void chameleon_flatodds(t_chameleon *x)
{
    int i;
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = 1.0;
    }
    setweights(x->odds,PROCESS_COUNT);
}

void chameleon_randodds(t_chameleon *x)
{
    int i;
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = boundrand(0.0,1.0);
    }
    setweights(x->odds,PROCESS_COUNT);
}

void *chameleon_new(t_symbol *msg, short argc, t_atom *argv)
{
    t_chameleon *x = (t_chameleon *)pd_new(chameleon_class);
    
    int i,j;
    t_float aLong;
    long max_dsp_units;
    long membytes = 0;
    srand((unsigned int)time(0));
    
    x->sr = sys_getsr();
    x->vs = sys_getblksize();
    // post("chameleon init: sr = %f, vs = %d", x->sr, x->vs);
    if(! x->sr){
        x->sr = 48000;
    }
    // one additional signal inlet
    inlet_new(&x->x_obj, &x->x_obj.ob_pd,gensym("signal"), gensym("signal"));
    // two signal outlets
    outlet_new(&x->x_obj, gensym("signal"));
    outlet_new(&x->x_obj, gensym("signal"));
    // create list outlet
    x->listo = outlet_new(&x->x_obj, gensym("list"));
    
    max_dsp_units = 3; // default value
    if(argc > 0){
        atom_arg_getfloat(&aLong, 0, argc, argv);
        if(aLong > 0){
            max_dsp_units = (long)aLong;
        }
    }
    x->max_dsp_units = max_dsp_units;
    // post("max dsp units: %d", x->max_dsp_units);
    x->sinelen = 65536;
    x->verbose = 0;
    x->stored_slot_count = 0;
    
    
    x->chan1buf = NULL;
    x->chan2buf = NULL;
    
    // OK HERE
    x->maxdelay = 1.0; // in seconds
    x->max_flangedelay = 0.1; // in seconds
    
    // memory allocation
    
    x->data = (t_atom *) getbytes(1024 * sizeof(t_atom)); // need to clean it
    membytes += 1024 * sizeof(t_atom);
    x->sinewave = (double *) getbytes( (x->sinelen + 1) * sizeof(double));
    membytes += (x->sinelen + 1) * sizeof(double);
    x->params = (double *) getbytes(MAX_PARAMETERS * sizeof(double));
    membytes += MAX_PARAMETERS * sizeof(double);
    x->odds = (float *) getbytes(64 * sizeof(float));
    x->distortion_length = 32768;
    x->distortion_function = (double *) getbytes(x->distortion_length * sizeof(double));
    membytes += x->distortion_length * sizeof(double);
    set_distortion_table(x->distortion_function, 0.1, 0.5, x->distortion_length);
    putsine(x->sinewave, x->sinelen);
    x->comb_delay_pool1 = (double **) getbytes(max_dsp_units * sizeof(double *));
    x->comb_delay_pool2 = (double **) getbytes(max_dsp_units * sizeof(double *));
    membytes += 2 * max_dsp_units * sizeof(double *);
    // might want to allow for sample rate change by putting this code into a re-init section
    for(i = 0; i < max_dsp_units; i++){
        x->comb_delay_pool1[i] = (double *) getbytes( ((x->maxdelay * x->sr) + 2) * sizeof(double));
        x->comb_delay_pool2[i] = (double *) getbytes( ((x->maxdelay * x->sr) + 2) * sizeof(double));
        membytes += 2 * ((x->maxdelay * x->sr) + 2) * sizeof(double);
    }
    
    // RINGMOD
    x->ringmod_phases = (double *) getbytes(max_dsp_units * sizeof(double));
    membytes += max_dsp_units * sizeof(double);
    // RINGMOD4
    x->ringmod4_phases = (double *) getbytes(max_dsp_units * sizeof(double));
    membytes += max_dsp_units * sizeof(double);
    
    // FLANGE
    
    x->flange_units = (t_flange_unit *) getbytes(max_dsp_units * sizeof(t_flange_unit));
    for(i = 0; i < max_dsp_units; i++){
        x->flange_units[i].flange_dl1 = (double *) getbytes(((x->maxdelay * x->sr) + 2) * sizeof(double));
        membytes += ((x->maxdelay * x->sr) + 2) * sizeof(double);
        x->flange_units[i].flange_dl2 = (double *) getbytes(((x->maxdelay * x->sr) + 2) * sizeof(double));
        membytes += ((x->maxdelay * x->sr) + 2) * sizeof(double);
        x->flange_units[i].dv1 = (int *) getbytes(2 * sizeof(int));
        x->flange_units[i].dv2 = (int *) getbytes(2 * sizeof(int));
        membytes += 4 * sizeof(int);
        x->flange_units[i].phase = 0.0;
        delset2(x->flange_units[i].flange_dl1, x->flange_units[i].dv1, x->max_flangedelay,x->sr);
        delset2(x->flange_units[i].flange_dl2, x->flange_units[i].dv2, x->max_flangedelay,x->sr);
    }
    
    // BUTTERWORTH
    
    x->butterworth_units = (t_butterworth_unit *) getbytes(max_dsp_units * sizeof(t_butterworth_unit));
    for(i = 0; i < max_dsp_units; i++){
        x->butterworth_units[i].data1 = (double *) getbytes(8 * sizeof(double));
        x->butterworth_units[i].data2 = (double *) getbytes(8 * sizeof(double));
        membytes += 16 * sizeof(double);
    }
    
    // TRUNCATE
    
    x->truncate_units = (t_truncate_unit *) getbytes(max_dsp_units * sizeof(t_truncate_unit));
    membytes += max_dsp_units * sizeof(t_truncate_unit);
    
    // SWEEPRESON
    
    x->sweepreson_units = (t_sweepreson_unit *) getbytes(max_dsp_units * sizeof(t_sweepreson_unit));
    membytes += max_dsp_units * sizeof(t_sweepreson_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->sweepreson_units[i].q1 = (double *) getbytes(5 * sizeof(double));
        x->sweepreson_units[i].q2 = (double *) getbytes(5 * sizeof(double));
        membytes += 10 * sizeof(double);
    }
    
    // SLIDECOMB
    
    x->slidecomb_units = (t_slidecomb_unit *) getbytes(max_dsp_units * sizeof(t_slidecomb_unit));
    membytes += max_dsp_units * sizeof(t_slidecomb_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->slidecomb_units[i].dv1 = (int *) getbytes(2 * sizeof(int));
        x->slidecomb_units[i].dv2 = (int *) getbytes(2 * sizeof(int));
        membytes += 4 * sizeof(int);
        x->slidecomb_units[i].delayline1 = (double *) getbytes( (2 + (MAX_SLIDECOMB_DELAY * x->sr)) * sizeof(double));
        x->slidecomb_units[i].delayline2 = (double *) getbytes( (2 + (MAX_SLIDECOMB_DELAY * x->sr)) * sizeof(double));
        membytes += 2 * (2 + (MAX_SLIDECOMB_DELAY * x->sr)) * sizeof(double);
        delset2(x->slidecomb_units[i].delayline1, x->slidecomb_units[i].dv1, MAX_SLIDECOMB_DELAY, x->sr);
        delset2(x->slidecomb_units[i].delayline2, x->slidecomb_units[i].dv2, MAX_SLIDECOMB_DELAY, x->sr);
    }
    
    // REVERB1
    
    x->reverb1_units = (t_reverb1_unit *) getbytes(max_dsp_units * sizeof(t_reverb1_unit));
    membytes += max_dsp_units * sizeof(t_reverb1_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->reverb1_units[i].dels = (double *) getbytes(4 * sizeof(double));
        x->reverb1_units[i].eel1 = (LSTRUCT *) getbytes(MAXSECTS * sizeof(LSTRUCT));
        x->reverb1_units[i].eel2 = (LSTRUCT *) getbytes(MAXSECTS * sizeof(LSTRUCT));
        x->reverb1_units[i].alpo1 = (double **) getbytes(4 * sizeof(double *));
        membytes += 2 * MAXSECTS * sizeof(LSTRUCT);
        membytes += 8 * sizeof(double *);
        for(j = 0; j < 4 ; j++ ){
            x->reverb1_units[i].alpo1[j] = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double));
            membytes += ((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double);
        }
        x->reverb1_units[i].alpo2 = (double **) getbytes(4 * sizeof(double *));
        for(j = 0; j < 4 ; j++ ){
            x->reverb1_units[i].alpo2[j] = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double));
            membytes += ((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double);
        }
    }
    x->reverb_ellipse_data = (double *) getbytes(16 * sizeof(double));
    membytes += 16 * sizeof(double);
    
    // ELLIPSEME
    
    x->ellipseme_units = (t_ellipseme_unit *) getbytes(max_dsp_units * sizeof(t_ellipseme_unit));
    for(i = 0; i < max_dsp_units; i++){
        x->ellipseme_units[i].eel1 = (LSTRUCT *) getbytes(MAXSECTS * sizeof(LSTRUCT));
        x->ellipseme_units[i].eel2 = (LSTRUCT *) getbytes(MAXSECTS * sizeof(LSTRUCT));
        membytes += 2 * MAXSECTS * sizeof(LSTRUCT);
    }
    
    x->ellipse_data = (double **) getbytes(MAXFILTER * sizeof(double *));
    membytes += MAXFILTER * sizeof(double *);
    for(i=0;i<MAXFILTER;i++){
        x->ellipse_data[i] = (double *) getbytes(MAX_COEF * sizeof(double));
        membytes += MAX_COEF * sizeof(double);
    }
    
    // FEED1
    
    x->feed1_units = (t_feed1_unit *) getbytes(max_dsp_units * sizeof(t_feed1_unit));
    membytes += max_dsp_units * sizeof(t_feed1_unit);
    x->feedfunclen = 8192;
    for(i = 0; i < max_dsp_units; i++){
        x->feed1_units[i].func1 = (double *)getbytes(2 * FEEDFUNCLEN * sizeof(double)); // doubling size, just in case
        x->feed1_units[i].func2 = (double *)getbytes(2 * FEEDFUNCLEN * sizeof(double));
        x->feed1_units[i].func3 = (double *)getbytes(2 * FEEDFUNCLEN * sizeof(double));
        x->feed1_units[i].func4 = (double *)getbytes(2 * FEEDFUNCLEN * sizeof(double));
        membytes += 8 * FEEDFUNCLEN * sizeof(double);
        x->feed1_units[i].delayLine1a = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        x->feed1_units[i].delayLine2a = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        x->feed1_units[i].delayLine1b = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        x->feed1_units[i].delayLine2b = (double *) getbytes(((int)(x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        membytes += 4 * ((int)(x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double);
        x->feed1_units[i].dv1a = (int *) getbytes(2 * sizeof(int));
        x->feed1_units[i].dv2a = (int *) getbytes(2 * sizeof(int));
        x->feed1_units[i].dv1b = (int *) getbytes(2 * sizeof(int));
        x->feed1_units[i].dv2b = (int *) getbytes(2 * sizeof(int));
        membytes += 8 * sizeof(int);
    }
    
    // BITCRUSH
    
    x->bitcrush_factors = (double *) getbytes(max_dsp_units * sizeof(double));
    membytes += max_dsp_units * sizeof(double);
    
    // FLAM1
    
    x->flam1_units = (t_flam1_unit *) getbytes(max_dsp_units * sizeof(t_flam1_unit));
    membytes += max_dsp_units * sizeof(t_flam1_unit);
    
    for(i = 0; i < max_dsp_units; i++){
        x->flam1_units[i].delayline1 = (double *) getbytes(((int)(x->sr * FLAM1_MAX_DELAY) + 5)  * sizeof(double));
        x->flam1_units[i].delayline2 = (double *) getbytes(((int)(x->sr * FLAM1_MAX_DELAY) + 5)  * sizeof(double));
        membytes += 2 * ((int)(x->sr * FLAM1_MAX_DELAY) + 5)  * sizeof(double);
        x->flam1_units[i].dv1 = (int *) getbytes(2 * sizeof(int));
        x->flam1_units[i].dv2 = (int *) getbytes(2 * sizeof(int));
        membytes += 4 * sizeof(int);
    }
    
    // COMB4
    
    x->comb4_units = (t_comb4_unit *) getbytes(max_dsp_units * sizeof(t_comb4_unit));
    membytes += max_dsp_units * sizeof(t_comb4_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->comb4_units[i].combs1 = (double **)getbytes(4 * sizeof(double *));
        x->comb4_units[i].combs2 = (double **)getbytes(4 * sizeof(double *));
        membytes += 8 * sizeof(double *);
        for(j = 0; j < 4; j++){
            x->comb4_units[i].combs1[j] = (double *)getbytes( (5 + TONECOMB_MAX_DELAY * x->sr) * sizeof(double));
            x->comb4_units[i].combs2[j] = (double *)getbytes( (5 + TONECOMB_MAX_DELAY * x->sr) * sizeof(double));
            membytes += 2 * (5 + TONECOMB_MAX_DELAY * x->sr) * sizeof(double);
        }
    }
    
    // RESONFEED
    
    x->resonfeed_units = (t_resonfeed_unit *) getbytes(max_dsp_units * sizeof(t_resonfeed_unit));
    membytes += max_dsp_units * sizeof(t_resonfeed_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->resonfeed_units[i].res1q = (double *)getbytes(5 * sizeof(double));
        x->resonfeed_units[i].res2q = (double *)getbytes(5 * sizeof(double));
        membytes += 10 * sizeof(double);
        x->resonfeed_units[i].comb1arr = (double *)getbytes(TONECOMB_MAX_DELAY * x->sr * sizeof(double));
        x->resonfeed_units[i].comb2arr = (double *)getbytes(TONECOMB_MAX_DELAY * x->sr * sizeof(double));
        membytes += 2 * TONECOMB_MAX_DELAY * x->sr * sizeof(double);
    }
    
    // RESONADSR
    
    x->resonadsr_units = (t_resonadsr_unit *) getbytes(max_dsp_units * sizeof(t_resonadsr_unit));
    membytes += max_dsp_units * sizeof(t_resonadsr_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->resonadsr_units[i].q1 = (double *)getbytes(5 * sizeof(double));
        x->resonadsr_units[i].q2 = (double *)getbytes(5 * sizeof(double));
        membytes += 10 * sizeof(double);
        x->resonadsr_units[i].adsr = (CMIXADSR *)getbytes(sizeof(CMIXADSR));
        membytes += sizeof(CMIXADSR);
        x->resonadsr_units[i].adsr->func = (double *)getbytes(8192 * sizeof(double));
        membytes += 8192 * sizeof(float);
        x->resonadsr_units[i].adsr->len = 8192;
    }
    
    // STV
    
    x->stv_units = (t_stv_unit *) getbytes(max_dsp_units * sizeof(t_stv_unit));
    membytes += max_dsp_units * sizeof(t_stv_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->stv_units[i].delayline1 = (double *) getbytes(((int)(x->sr * STV_MAX_DELAY) + 5)  * sizeof(double));
        x->stv_units[i].delayline2 = (double *) getbytes(((int)(x->sr * STV_MAX_DELAY) + 5)  * sizeof(double));
        membytes += 2 * ((int)(x->sr * STV_MAX_DELAY) + 5)  * sizeof(double);
        x->stv_units[i].dv1 = (int *) getbytes(2 * sizeof(int));
        x->stv_units[i].dv2 = (int *) getbytes(2 * sizeof(int));
        membytes += 4 * sizeof(int);
    }
    
    // BENDY
    
    x->bendy_units = (t_bendy_unit *) getbytes(max_dsp_units * sizeof(t_bendy_unit));
    membytes += max_dsp_units * sizeof(t_bendy_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->bendy_units[i].delayline1 = (double *) getbytes(((int)(x->sr * BENDY_MAXDEL) + 5)  * sizeof(double));
        x->bendy_units[i].delayline2 = (double *) getbytes(((int)(x->sr * BENDY_MAXDEL) + 5)  * sizeof(double));
        membytes += 2 * ((int)(x->sr * BENDY_MAXDEL) + 5)  * sizeof(double);
        x->bendy_units[i].dv1 = (int *) getbytes(8 * sizeof(int));
        x->bendy_units[i].dv2 = (int *) getbytes(8 * sizeof(int));
        membytes += 4 * sizeof(int);
        x->bendy_units[i].val1 = boundrand(0.01, BENDY_MAXDEL);
        x->bendy_units[i].val2 = boundrand(0.01, BENDY_MAXDEL);
        x->bendy_units[i].counter = 0;
    }
    
    // SLIDEFLAM
    
    x->slideflam_units = (t_slideflam_unit *) getbytes(max_dsp_units * sizeof(t_slideflam_unit));
    membytes += max_dsp_units * sizeof(t_slideflam_unit);
    for(i = 0; i < max_dsp_units; i++){
        x->slideflam_units[i].dv1 = (int *) getbytes(2 * sizeof(int));
        x->slideflam_units[i].dv2 = (int *) getbytes(2 * sizeof(int));
        membytes += 4 * sizeof(int);
        x->slideflam_units[i].delayline1 = (double *) getbytes( (10 + (MAX_SLIDEFLAM_DELAY * x->sr)) * sizeof(double));
        x->slideflam_units[i].delayline2 = (double *) getbytes( (10 + (MAX_SLIDEFLAM_DELAY * x->sr)) * sizeof(double));
        membytes += 2 * (10 + (MAX_SLIDEFLAM_DELAY * x->sr)) * sizeof(double);
        delset2(x->slideflam_units[i].delayline1, x->slideflam_units[i].dv1, MAX_SLIDEFLAM_DELAY, x->sr);
        delset2(x->slideflam_units[i].delayline2, x->slideflam_units[i].dv2, MAX_SLIDEFLAM_DELAY, x->sr);
    }
    
    x->tf_len = 1;
    x->tf_len <<= 16;
    
    setflamfunc1(x->flamfunc1,x->flamfunc1len);
    x->max_comb_lpt = 0.15 ;// watch out here
    
    x->adsr = (CMIXADSR *) getbytes(1 * sizeof(CMIXADSR));
    membytes += sizeof(CMIXADSR);
    x->adsr->len = 32768;
    x->adsr->func = (double *) getbytes(x->adsr->len * sizeof(double));
    membytes += x->adsr->len * sizeof(double);
    x->dcflt = (double *) getbytes(16 * sizeof(double));
    membytes += 16 * sizeof(double);
    
    
    // allocate memory to store random patterns
    x->slots = (t_slot *) getbytes( MAX_SLOTS * sizeof(t_slot) );
    membytes += MAX_SLOTS * sizeof(t_slot);
    for(i = 0; i < MAX_SLOTS; i++){
        x->slots[i].params = (double *) getbytes( MAX_PARAMETERS * sizeof(double));
        membytes += MAX_PARAMETERS * sizeof(double);
    }
        
    /* be sure to finish clearing memory */
    set_dcflt(x->dcflt); // WE NEED THIS FILTER!
    init_reverb_data(x->reverb_ellipse_data);
    init_ellipse_data(x->ellipse_data);
    
    for(i=0;i<PROCESS_COUNT;i++){
        x->odds[i] = 1;
    }
    x->min_process_per_note = 1;
    x->max_process_per_note = 1;
    setweights(x->odds,PROCESS_COUNT);
    
    x->ratios = (double *) getbytes(5 * sizeof(double));
    membytes += 5 * sizeof(double);
    x->ratios[0] = 9.0/8.0;
    x->ratios[1] = 4.0/3.0;
    x->ratios[2] = 5.0/4.0;
    x->ratios[3] = 6.0/5.0;
    x->ratios[4] = 3.0/2.0;
    
    x->membytes = membytes;
    
    chameleon_set_parameters(x);
    return x;
}

void chameleon_memory(t_chameleon *x)
{
    post("Memory allocation for this instance of chameleon~ is %.2f MBytes",(float)x->membytes/1000000.0);
}

void chameleon_dsp_free(t_chameleon *x)
{
    int i,j;
    long max_dsp_units = x->max_dsp_units;
    
    freebytes(x->data,1024 * sizeof(t_atom));
    freebytes(x->sinewave,(x->sinelen + 1) * sizeof(double));
    freebytes(x->params,MAX_PARAMETERS * sizeof(double));
    freebytes(x->odds,64 * sizeof(float));
    
    // free combs
    for(i= 0; i < max_dsp_units; i++){
        freebytes( x->comb_delay_pool1[i],((x->maxdelay * x->sr) + 2) * sizeof(double));
        freebytes( x->comb_delay_pool2[i],((x->maxdelay * x->sr) + 2) * sizeof(double));
    }
    freebytes(x->comb_delay_pool1,max_dsp_units * sizeof(double *));
    freebytes(x->comb_delay_pool2,max_dsp_units * sizeof(double *));
    freebytes(x->eel, MAXSECTS * sizeof(LSTRUCT));
    
    // free ringmod
    freebytes(x->ringmod_phases, max_dsp_units * sizeof(double));
    
    // free ringmod4
    freebytes(x->ringmod4_phases, max_dsp_units * sizeof(double));

    // free flange
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->flange_units[i].flange_dl1,((x->maxdelay * x->sr) + 2) * sizeof(double));
        freebytes(x->flange_units[i].flange_dl2,((x->maxdelay * x->sr) + 2) * sizeof(double));
        freebytes(x->flange_units[i].dv1, 2 * sizeof(int));
        freebytes(x->flange_units[i].dv2, 2 * sizeof(int));
    }
    freebytes(x->flange_units, max_dsp_units * sizeof(t_flange_unit));
    
    // free butterworth
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->butterworth_units[i].data1, 8 * sizeof(double));
        freebytes(x->butterworth_units[i].data2, 8 * sizeof(double));
    }
    freebytes(x->butterworth_units,max_dsp_units * sizeof(t_butterworth_unit));
    
    // free truncate
    freebytes(x->truncate_units,max_dsp_units * sizeof(t_truncate_unit));
    
    // free sweepreson
    for(i = 0; i < max_dsp_units; i++){
        freebytes( x->sweepreson_units[i].q1,5 * sizeof(double));
        freebytes( x->sweepreson_units[i].q2,5 * sizeof(double));
    }
    freebytes(x->sweepreson_units,max_dsp_units * sizeof(t_sweepreson_unit));
    
    // free ellipse
    freebytes(x->reverb_ellipse_data,16 * sizeof(double));
    for(i=0;i<MAXFILTER;i++){
        freebytes(x->ellipse_data[i],MAX_COEF * sizeof(double));
    }
    freebytes(x->ellipse_data, MAXFILTER * sizeof(double *));
    // free bendy
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->bendy_units[i].delayline1, ((int)(x->sr * BENDY_MAXDEL) + 5)  * sizeof(double));
        freebytes(x->bendy_units[i].delayline2, ((int)(x->sr * BENDY_MAXDEL) + 5)  * sizeof(double));
        freebytes(x->bendy_units[i].dv1,8 * sizeof(int));
        freebytes(x->bendy_units[i].dv2,8 * sizeof(int));
    }
    freebytes(x->bendy_units,max_dsp_units * sizeof(t_bendy_unit));
    
    // free slidecomb
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->slidecomb_units[i].delayline1,(2 + (MAX_SLIDECOMB_DELAY * x->sr)) * sizeof(double));
        freebytes(x->slidecomb_units[i].delayline2,(2 + (MAX_SLIDECOMB_DELAY * x->sr)) * sizeof(double));
        freebytes(x->slidecomb_units[i].dv1,2 * sizeof(int));
        freebytes(x->slidecomb_units[i].dv2,2 * sizeof(int));
    }
    freebytes(x->slidecomb_units, max_dsp_units * sizeof(t_slidecomb_unit));
    
    // free reverb
    
    for(i = 0; i < max_dsp_units; i++){
        for(j = 0; j < 4 ; j++ ){
            freebytes(x->reverb1_units[i].alpo1[j], ((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double));
        }
        freebytes(x->reverb1_units[i].alpo1,4 * sizeof(double *));
        for(j = 0; j < 4 ; j++ ){
            freebytes(x->reverb1_units[i].alpo2[j],((int)(x->sr * MAX_MINI_DELAY) + 1)  * sizeof(double));
        }
        freebytes(x->reverb1_units[i].alpo2,4 * sizeof(double *));
        freebytes(x->reverb1_units[i].dels,4 * sizeof(double));
        freebytes(x->reverb1_units[i].eel1,MAXSECTS * sizeof(LSTRUCT));
        freebytes(x->reverb1_units[i].eel2,MAXSECTS * sizeof(LSTRUCT));
    }

/*
 ANOMALY - this free statement is a crasher:
    freebytes(x->reverb_ellipse_data,16 * sizeof(double));
     */
    freebytes(x->reverb1_units,max_dsp_units * sizeof(t_reverb1_unit));


    // free feed1
    for(i = 0; i < max_dsp_units; i++){

        freebytes(x->feed1_units[i].func1,2 * FEEDFUNCLEN * sizeof(double));
        freebytes(x->feed1_units[i].func2,2 * FEEDFUNCLEN * sizeof(double));
        freebytes(x->feed1_units[i].func3,2 * FEEDFUNCLEN * sizeof(double));
        freebytes(x->feed1_units[i].func4,2 * FEEDFUNCLEN * sizeof(double));

        freebytes(x->feed1_units[i].delayLine1a,(int)((x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        freebytes(x->feed1_units[i].delayLine2a,(int)((x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        freebytes(x->feed1_units[i].delayLine1b,(int)((x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        freebytes(x->feed1_units[i].delayLine2b,(int)((x->sr * MAX_MINI_DELAY) + 5)  * sizeof(double));
        freebytes(x->feed1_units[i].dv1a,2 * sizeof(int));
        freebytes(x->feed1_units[i].dv2a,2 * sizeof(int));
        freebytes(x->feed1_units[i].dv1b,2 * sizeof(int));
        freebytes(x->feed1_units[i].dv2b,2 * sizeof(int));

    }

    freebytes(x->feed1_units,max_dsp_units * sizeof(t_feed1_unit));

    // free bitcrush
    freebytes(x->bitcrush_factors,max_dsp_units * sizeof(double));
 
    // free flam1
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->flam1_units[i].delayline1,((int)(x->sr * FLAM1_MAX_DELAY) + 5)  * sizeof(double));
        freebytes(x->flam1_units[i].delayline2,((int)(x->sr * FLAM1_MAX_DELAY) + 5)  * sizeof(double));
        freebytes(x->flam1_units[i].dv1, 2 * sizeof(int));
        freebytes(x->flam1_units[i].dv2, 2 * sizeof(int));
    }
    freebytes(x->flam1_units,max_dsp_units * sizeof(t_flam1_unit));
    
    // free comb4
    
    for(i = 0; i < max_dsp_units; i++){

        for(j = 0; j < 4; j++){
            freebytes(x->comb4_units[i].combs1[j],(5 + TONECOMB_MAX_DELAY * x->sr) * sizeof(double) );
            freebytes(x->comb4_units[i].combs2[j],(5 + TONECOMB_MAX_DELAY * x->sr) * sizeof(double));
        }
     freebytes(x->comb4_units[i].combs1,4 * sizeof(double *));
     freebytes(x->comb4_units[i].combs2,4 * sizeof(double *));
    }
    freebytes(x->comb4_units, max_dsp_units * sizeof(t_comb4_unit));
    
    // free resonfeed
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->resonfeed_units[i].res1q,5 * sizeof(double));
        freebytes(x->resonfeed_units[i].res2q,5 * sizeof(double));
        freebytes(x->resonfeed_units[i].comb1arr,TONECOMB_MAX_DELAY * x->sr * sizeof(double));
        freebytes(x->resonfeed_units[i].comb2arr,TONECOMB_MAX_DELAY * x->sr * sizeof(double));
    }
    freebytes(x->resonfeed_units,max_dsp_units * sizeof(t_resonfeed_unit));

    // free resonadsr
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->resonadsr_units[i].q1,5 * sizeof(double));
        freebytes(x->resonadsr_units[i].q2,5 * sizeof(double));
        freebytes(x->resonadsr_units[i].adsr,sizeof(CMIXADSR));
        freebytes(x->resonadsr_units[i].adsr->func,8192 * sizeof(double));
    }
    freebytes(x->resonadsr_units,max_dsp_units * sizeof(t_resonadsr_unit));
    
    // free stv
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->stv_units[i].delayline1,((int)(x->sr * STV_MAX_DELAY) + 5)  * sizeof(double));
        freebytes(x->stv_units[i].delayline2,((int)(x->sr * STV_MAX_DELAY) + 5)  * sizeof(double));
        freebytes(x->stv_units[i].dv1,2 * sizeof(int));
        freebytes(x->stv_units[i].dv2,2 * sizeof(int));
    }
    freebytes(x->stv_units,max_dsp_units * sizeof(t_stv_unit));
    

    // free slideflam
    
    for(i = 0; i < max_dsp_units; i++){
        freebytes(x->slideflam_units[i].delayline1,(10 + (MAX_SLIDEFLAM_DELAY * x->sr)) * sizeof(double));
        freebytes(x->slideflam_units[i].delayline2,(10 + (MAX_SLIDEFLAM_DELAY * x->sr)) * sizeof(double));
        freebytes(x->slideflam_units[i].dv1,2 * sizeof(int));
        freebytes(x->slideflam_units[i].dv2,2 * sizeof(int));
    }
    freebytes(x->slideflam_units,max_dsp_units * sizeof(t_slideflam_unit));

    // free a few other misc items
    freebytes( x->adsr,sizeof(CMIXADSR));
    freebytes(x->adsr->func,x->adsr->len * sizeof(double));
    freebytes(x->dcflt,16 * sizeof(double));
 
    // free slot memory
    for(i = 0; i < MAX_SLOTS; i++){
        freebytes(x->slots[i].params, MAX_PARAMETERS * sizeof(double));
    }
    freebytes(x->slots, MAX_SLOTS * sizeof(t_slot));
}


/*
void chameleon_perform64(t_chameleon *x, t_object *dsp64, double **ins,
                         long numins, double **outs,long numouts, long n,
                         long flags, void *userparam)
*/
t_int *chameleon_perform(t_int *w)
{
    t_chameleon *x = (t_chameleon *) (w[1]);
    t_float *OGchan1 = (t_float *) (w[2]);
    t_float *OGchan2 = (t_float *) (w[3]);
    t_float *outchanL = (t_float *) (w[4]);
    t_float *outchanR = (t_float *) (w[5]);
    int n = (int) w[6];
    long i;
    long pcount;
    double *params;
    long curarg;
    t_double *chan1, *chan2;
    if( x->set_parameters_flag == 1){
        chameleon_set_parameters_exec(x);
    }
    else if(x->recall_parameters_flag == 1)  {
        chameleon_recall_parameters_exec(x);
    }
    
    /* NOTE: operating directly on the input vectors resulted in those vectors (dry signal) being replaced by
     processed signal. We need to copy to local arrays to avoid signal contamination. */
    
    if( x->chan1buf == NULL){
        // post("allocated vector inside perform routine");
        x->chan1buf = (t_double *) getbytes(sizeof(t_double) * n);
        x->chan2buf = (t_double *) getbytes(sizeof(t_double) * n);
    }
    chan1 = x->chan1buf;
    chan2 = x->chan2buf;
    
    
    
    /*
     Copy input vectors to work vectors (chan1 & chan2). Also, promote single-precision floats from Pure Data to
     double precision for Chameleon processing.
     */
    for(i = 0; i < n; i++){
        chan1[i] = (double)OGchan1[i];
        chan2[i] = (double)OGchan2[i];
    }
    
    pcount = x->pcount;
    params = x->params;
    curarg = 0;
    while(curarg < pcount){
        if(params[curarg] == BENDY){
            bendy(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == RINGMOD){
            ringmod(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == BITCRUSH){
            bitcrush(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == COMB){
            comber(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == FLANGE){
            flange(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == BUTTER){
            butterme(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == TRUNCATE){
            truncateme(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == SWEEPRESON){
            sweepreson(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == SLIDECOMB){
            slidecomb(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == REVERB1){
            reverb1(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == ELLIPSE){
            ellipseme(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == FEED1){
            feed1me(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == FLAM1){
            flam1(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == SLIDEFLAM){
            slideflam(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == RINGMOD4){
            ringmod4(x,&curarg,chan1,chan2);
        }
        else if(params[curarg] == COMB4){
            comb4(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == COMPDIST){
            compdist(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == RINGFEED){
            ringfeed(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == RESONADSR){
            resonadsr(x, &curarg, chan1, chan2);
        }
        else if(params[curarg] == STV){
            stv(x, &curarg, chan1, chan2);
        }
        else {
            error("chameleon~: programming error - deploy missing branch for %ld", curarg);
            goto panic;
        }
    }
    // HERE GOES THE OUTPUT. Downgrade Chameleon doubles to Pd single-precision floats
    for(i = 0; i < n; i++){
        outchanL[i] = (float) chan1[i];
        outchanR[i] = (float) chan2[i];
    }
panic:
    ;
    return w + 7;
}

void chameleon_set_parameters(t_chameleon *x){
    x->set_parameters_flag = 1;
}

void clip(double *x, double min, double max){
    if(*x < min){
        *x = min;
    } else if (*x > max){
        *x = max;
    }
}

void chameleon_tweak_parameters(t_chameleon *x, t_floatarg tdev){
    int i, j;
    int ftype;
    double cf, bw;//, bw;
    double delay, revtime;
    double *params = x->params;
    long pcount = x->pcount;
    int comb_dl_count = 0; // where in the comb pool to grab memory
    int flange_count = 0;
    int truncate_count = 0;
    int butterworth_count = 0;
    int sweepreson_count = 0;
    int slidecomb_count = 0;
    int reverb1_count = 0;
    int ellipseme_count = 0;
    int feed1_count = 0;
    int flam1_count = 0;
    int comb4_count = 0;
    int ringfeed_count = 0;
    int bendy_count = 0;
    int ringmod_count = 0;
    int ringmod4_count = 0;
    int slideflam_count = 0;
    int bitcrush_count = 0;
    int resonadsr_count = 0;
    int stv_count = 0;
    double raw_wet;
    double sr = x->sr;
    double *dels;
    double **alpo1, **alpo2;
    double speed1, speed2, mindelay, maxdelay, duration;
    double basefreq;
    double rvt = boundrand(0.1,0.98);
    double lpt;
    double notedur, sust;
    double phz1, phz2, minfeedback = 0.1, maxfeedback = 0.7;
    long slotnum = x->recall_slot;
    double *rparams; // parameter set to recall
    long rpcount; // number of parameters to read
    long slot_pcount;
    double multmin, multmax, tmp;
    
    tmp = tdev;
    clip(&tmp, 0.001, 0.5);
    tdev = tmp;
    multmin = 1.0 - tdev;
    multmax = 1.0 + tdev;
    /* Now iterate through the recall, but with no setting of memory */
    slot_pcount = x->slots[slotnum].pcount;
    
    if( slot_pcount <= 0){
        //    post("Aborting reload of slot %d", slotnum);
        return;
    }
    rparams = x->slots[slotnum].params;
    pcount = 0;
    rpcount = 0;

    while( rpcount < slot_pcount ){
        j = rparams[rpcount++];
        if(j == COMB){
            params[pcount++] = COMB;
            delay = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&delay, 0.01, 0.35);
            params[pcount++] = delay;
            revtime = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&revtime,0.5,0.98);
            params[pcount++] = revtime;
            params[pcount++] = comb_dl_count = rparams[rpcount++]; // Possible Bug Here???
            mycombset(delay,revtime,0,x->comb_delay_pool1[comb_dl_count],x->sr);
            mycombset(delay,revtime,0,x->comb_delay_pool2[comb_dl_count],x->sr);
        }
        else if(j == RINGMOD) {
            // post("Added RINGMOD unit");
            params[pcount++] = RINGMOD;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 40.0, 3000.0 );
            params[pcount++] = tmp;
            params[pcount++] = ringmod_count = rparams[rpcount++];
        }
        else if(j == RINGMOD4) {
            // post("Added RINGMOD4 unit");
            params[pcount++] = RINGMOD4;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 40.0, 3000.0 );
            params[pcount++] = tmp;
            params[pcount++] = ringmod4_count = rparams[rpcount++];
        }
        else if(j == BENDY){
            params[pcount++] = BENDY;
            params[pcount++] = bendy_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, BENDY_MAXDEL);
            params[pcount++] = x->bendy_units[bendy_count].val1 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, BENDY_MAXDEL);
            params[pcount++] = x->bendy_units[bendy_count].val2 = tmp;
            // danger here
            x->bendy_units[bendy_count].counter = 0;
            delset2(x->bendy_units[bendy_count].delayline1, x->bendy_units[bendy_count].dv1, x->bendy_units[bendy_count].val1,x->sr);
            delset2(x->bendy_units[bendy_count].delayline2, x->bendy_units[bendy_count].dv2, x->bendy_units[bendy_count].val2,x->sr);
        }
        else if(j == FLANGE){
            params[pcount++] = FLANGE;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 80.0, 400.0);
            params[pcount++] = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 600.0, 4000.0);
            params[pcount++] = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.1, 2.0);
            params[pcount++] = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.1, 0.95);
            params[pcount++] = tmp;
            params[pcount++] = flange_count = rparams[rpcount++];
        }
        else if(j == BUTTER){
            /*
             params[pcount++] = cf = boundrand(70.0,3000.0);
             params[pcount++] = bw = cf * boundrand(0.05,0.6);
             */
            params[pcount++] = BUTTER;
            params[pcount++] = ftype = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 70.0, 3000.0);
            params[pcount++] = cf = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.05, 0.6);
            params[pcount++] = bw = tmp;
            params[pcount++] = butterworth_count = rparams[rpcount++];
            if( ftype == LOPASS) {
                lobut(x->butterworth_units[butterworth_count].data1, cf, sr);
                lobut(x->butterworth_units[butterworth_count].data2, cf, sr);
            } else if (ftype == HIPASS){
                hibut(x->butterworth_units[butterworth_count].data1, cf, sr);
                hibut(x->butterworth_units[butterworth_count].data2, cf, sr);
            }
            else if(ftype == BANDPASS){
                bpbut(x->butterworth_units[butterworth_count].data1, cf, bw, sr);
                bpbut(x->butterworth_units[butterworth_count].data2, cf, bw, sr);
                x->butterworth_units[butterworth_count].bw = bw;
            }
            x->butterworth_units[butterworth_count].cf = cf;
            x->butterworth_units[butterworth_count].ftype = ftype;
        }
        else if(j == TRUNCATE){
            // nothing to change here
            params[pcount++] = TRUNCATE;
            params[pcount++] = truncate_count = rparams[rpcount++];
            /*
            x->truncate_units[truncate_count].counter = 0;
            x->truncate_units[truncate_count].state = 0;
            x->truncate_units[truncate_count].segsamples = 1;
            */
        }
        else if(j == SWEEPRESON){
            params[pcount++] = SWEEPRESON;
            params[pcount++] = sweepreson_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 100.0, 300.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].minfreq = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 600.0, 6000.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].maxfreq = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, 0.2);
            params[pcount++] = x->sweepreson_units[sweepreson_count].bwfac = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.05, 2.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].speed = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.0, 0.5);
            params[pcount++] = x->sweepreson_units[sweepreson_count].phase = tmp;
            x->sweepreson_units[sweepreson_count].q1[3] = 0;
            x->sweepreson_units[sweepreson_count].q1[4] = 0;
            x->sweepreson_units[sweepreson_count].q2[3] = 0;
            x->sweepreson_units[sweepreson_count].q2[4] = 0;
        }
        else if(j == SLIDECOMB){
            params[pcount++] = SLIDECOMB;
            params[pcount++] = slidecomb_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.001,MAX_SLIDECOMB_DELAY * 0.95);
            params[pcount++] = x->slidecomb_units[slidecomb_count].start_delay = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.001,MAX_SLIDECOMB_DELAY * 0.95);
            params[pcount++] = x->slidecomb_units[slidecomb_count].end_delay = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.7,0.99);
            params[pcount++] = x->slidecomb_units[slidecomb_count].feedback = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.1,2.0);
            params[pcount++] = x->slidecomb_units[slidecomb_count].sample_length = tmp;
            // scale length to SR
            x->slidecomb_units[slidecomb_count].sample_length *= x->sr;
            x->slidecomb_units[slidecomb_count].counter = 0;
            delset2(x->slidecomb_units[slidecomb_count].delayline1,x->slidecomb_units[slidecomb_count].dv1,MAX_SLIDECOMB_DELAY, x->sr);
            delset2(x->slidecomb_units[slidecomb_count].delayline2,x->slidecomb_units[slidecomb_count].dv2,MAX_SLIDECOMB_DELAY, x->sr);
        }
        else if(j == REVERB1){
            params[pcount++] = REVERB1;
            params[pcount++] = reverb1_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.25,0.99);
            params[pcount++] = revtime = x->reverb1_units[reverb1_count].revtime = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.2,0.8);
            params[pcount++] = raw_wet = tmp;
            x->reverb1_units[reverb1_count].wet = sin(1.570796 * raw_wet);
            x->reverb1_units[reverb1_count].dry = cos(1.570796 * raw_wet);
            dels = x->reverb1_units[reverb1_count].dels;
            for(i = 0; i < 4; i++){
                tmp = rparams[rpcount++] * boundrand(multmin, multmax);
                clip(&tmp, 0.005, 0.1);
                params[pcount++] = dels[i] = tmp;
            }
            alpo1 = x->reverb1_units[reverb1_count].alpo1;
            alpo2 = x->reverb1_units[reverb1_count].alpo2;
            for( i = 0; i < 4; i++ ){
                if(dels[i] < .005 || dels[i] > 0.1) {
                    pd_error((t_object *)x,"reverb1: bad random delay time: %f",dels[i]);
                    dels[i] = .05;
                }
                // could be dangerous!
                mycombset(dels[i], revtime, 0, alpo1[i], x->sr);
                mycombset(dels[i], revtime, 0, alpo2[i], x->sr);
            }
        }
        else if(j == ELLIPSE){
           // nothing to change here
            params[pcount++] = ELLIPSE;
            params[pcount++] = ellipseme_count = rparams[rpcount++];
            params[pcount++] = x->ellipseme_units[ellipseme_count].filtercode = rparams[rpcount++];
            
            if( x->ellipseme_units[ellipseme_count].filtercode >= ELLIPSE_FILTER_COUNT ){
                error("chameleon~: there is no %d ellipse data",x->ellipseme_units[ellipseme_count].filtercode);
                return;
            };
        }
        else if(j == FEED1){
            params[pcount++] = FEED1;
            params[pcount++] = feed1_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.001, 0.1);
            params[pcount++] = mindelay = x->feed1_units[feed1_count].mindelay = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, mindelay, 0.1);
            params[pcount++] = maxdelay = x->feed1_units[feed1_count].maxdelay = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, 0.5);
            params[pcount++] = speed1 = x->feed1_units[feed1_count].speed1 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, speed1, 0.5);
            params[pcount++] = speed2 = x->feed1_units[feed1_count].speed2 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.05, 1.0);
            params[pcount++] = duration = x->feed1_units[feed1_count].duration = tmp;
            
            funcgen1(x->feed1_units[feed1_count].func1, FEEDFUNCLEN,duration,
                     mindelay,maxdelay, speed1, speed2, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func2, FEEDFUNCLEN,duration,
                     mindelay * 0.5,maxdelay * 2.0, speed1 * 1.25, speed2 * 0.75, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func3, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.35, speed2*1.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func4, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.55, speed2*2.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
        }
        else if(j == BITCRUSH){
            params[pcount++] = BITCRUSH;
            params[pcount++] = bitcrush_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 2.0, 8.0);
            params[pcount++] = x->bitcrush_factors[bitcrush_count] = tmp;
        }
        else if(j == FLAM1){
            params[pcount++] = FLAM1;
            params[pcount++] = flam1_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax); // not yet scaled to SR
            clip(&tmp, 1.0, 4.0);
            params[pcount++] = x->flam1_units[flam1_count].sample_length = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.08, 0.2);
            params[pcount++] = x->flam1_units[flam1_count].dt = tmp;
            
            x->flam1_units[flam1_count].counter = 0;
            // scale to SR
            x->flam1_units[flam1_count].sample_length *= x->sr;
            delset2(x->flam1_units[flam1_count].delayline1, x->flam1_units[flam1_count].dv1, MAX_MINI_DELAY, x->sr);
            delset2(x->flam1_units[flam1_count].delayline2, x->flam1_units[flam1_count].dv2, MAX_MINI_DELAY, x->sr);
        }
        else if(j == SLIDEFLAM){
            double sf_del1, sf_del2;
            params[pcount++] = SLIDEFLAM;
            params[pcount++] = slideflam_count = rparams[rpcount++];
            sf_del1 = rparams[rpcount++] * boundrand(multmin, multmax);
            sf_del2 = rparams[rpcount++] * boundrand(multmin, multmax);
            if( sf_del2 > sf_del1 ){
                clip(&sf_del1, 0.01, 0.05);
                clip(&sf_del2, 0.1,MAX_SLIDEFLAM_DELAY * 0.95);
            } else {
                clip(&sf_del2, 0.01, 0.05);
                clip(&sf_del1, 0.1,MAX_SLIDEFLAM_DELAY * 0.95);
            }
            params[pcount++] = x->slideflam_units[slideflam_count].dt1 = sf_del1;
            params[pcount++] = x->slideflam_units[slideflam_count].dt2 = sf_del2;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.7, 0.99);
            params[pcount++] = x->slideflam_units[slideflam_count].feedback = tmp;
            params[pcount++] = x->slideflam_units[slideflam_count].sample_length = rparams[rpcount++] * boundrand(multmin, multmax); // not yet scaled to SR
            
            x->slideflam_units[slideflam_count].counter = 0;
            // scale to SR
            x->slideflam_units[slideflam_count].sample_length *= x->sr;
            delset2(x->slideflam_units[slideflam_count].delayline1,x->slideflam_units[slideflam_count].dv1,MAX_SLIDEFLAM_DELAY, x->sr);
            delset2(x->slideflam_units[slideflam_count].delayline2,x->slideflam_units[slideflam_count].dv2,MAX_SLIDEFLAM_DELAY, x->sr);
        }
        else if(j == COMB4){
            params[pcount++] = COMB4;
            params[pcount++] = comb4_count = rparams[rpcount++];
            rvt = 0.99;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 100.0,400.0);
            params[pcount++] = basefreq = tmp;
            // post("comb4: count: %d, basefreq: %f",comb4_count, basefreq);
            for( i = 0; i < 4; i++ ){
                tmp = rparams[rpcount++] * boundrand(multmin, multmax);
                clip(&tmp, 0.0001,0.05);
                params[pcount++] = lpt = tmp;
                // dangerous ?
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs1[i], x->sr);
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs2[i], x->sr);
            }
        }
        else if(j == COMPDIST){
            /* nothing to change */
            params[pcount++] = COMPDIST;
        }
        else if(j == RINGFEED){
            params[pcount++] = RINGFEED;
            params[pcount++] = ringfeed_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 90.0, 1500.0);
            params[pcount++] = cf = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.02 * cf, 0.4 * cf);
            params[pcount++] = bw = tmp;
            
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res1q, x->sr);
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res2q, x->sr);
            x->resonfeed_units[ringfeed_count].osc1phs = 0.0;
            x->resonfeed_units[ringfeed_count].osc2phs = 0.0;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 90.0, 1500.0);
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc1si = tmp; // not yet scaled to SR
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 90.0, 1500.0);
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc2si = tmp; // not yet scaled to SR
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 1.0/1500.0, 1.0/90.0);
            params[pcount++] = lpt = tmp;
            tmp = rparams[rpcount++];
            clip(&tmp, 0.01, 0.8);
            params[pcount++] = rvt = tmp;
            // scale to SR
            x->resonfeed_units[ringfeed_count].osc1si *= ((double)x->sinelen / x->sr);
            x->resonfeed_units[ringfeed_count].osc2si *= ((double)x->sinelen / x->sr);
            // dangerous ?
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb1arr, x->sr);
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb2arr, x->sr);

        }
        else if(j == RESONADSR){
            params[pcount++] = RESONADSR;
            params[pcount++] = resonadsr_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, 0.1);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->a = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.01, 0.05);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->d = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.05, 0.5);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->r = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 150.0, 4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v1 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 150.0, 4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v2 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 150.0, 4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v3 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 150.0, 4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v4 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.03,0.7);
            params[pcount++] = x->resonadsr_units[resonadsr_count].bwfac = tmp;
            tmp =  rparams[rpcount++] * boundrand(multmin, multmax);
            params[pcount++] =  notedur = tmp;
            clip(&tmp, 0.7,1.2);
            x->resonadsr_units[resonadsr_count].phs = 0.0;
            sust = notedur - (x->resonadsr_units[resonadsr_count].adsr->a + x->resonadsr_units[resonadsr_count].adsr->d +  x->resonadsr_units[resonadsr_count].adsr->r);
            x->resonadsr_units[resonadsr_count].adsr->s = sust;
            buildadsr(x->resonadsr_units[resonadsr_count].adsr);
            x->resonadsr_units[resonadsr_count].si = ((double)x->resonadsr_units[resonadsr_count].adsr->len / x->sr) / notedur;
        }
        else if(j == STV){

            params[pcount++] = STV;
            params[pcount++] = stv_count = rparams[rpcount++];
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.025,0.5);
            params[pcount++] = speed1 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.025,0.5);
            params[pcount++] = speed2 = tmp;
            tmp = rparams[rpcount++] * boundrand(multmin, multmax);
            clip(&tmp, 0.001,0.01);
            params[pcount++] = maxdelay = tmp;
            
            x->stv_units[stv_count].osc1phs = 0.0;
            x->stv_units[stv_count].osc2phs = 0.0;
            x->stv_units[stv_count].fac2 = 0.5 * (maxdelay - 0.001);
            x->stv_units[stv_count].fac1 = 0.001 + x->stv_units[stv_count].fac2;
            x->stv_units[stv_count].osc1si = ((double)x->sinelen / x->sr) * speed1;
            x->stv_units[stv_count].osc2si = ((double)x->sinelen / x->sr) * speed2;
            delset2(x->stv_units[stv_count].delayline1, x->stv_units[stv_count].dv1, maxdelay, x->sr);
            delset2(x->stv_units[stv_count].delayline2, x->stv_units[stv_count].dv2, maxdelay, x->sr);
        }
        else {
            error("el.chameleon~: could not find a process for %d",j);
        }
    }
    x->pcount = pcount;
}

void chameleon_recall_parameters_exec(t_chameleon *x)
{
    int i, j;
    int ftype;
    double cf, bw;//, bw;
    double delay, revtime;
    double *params = x->params;
    long pcount = x->pcount;
    int comb_dl_count = 0; // where in the comb pool to grab memory
    int flange_count = 0;
    int truncate_count = 0;
    int butterworth_count = 0;
    int sweepreson_count = 0;
    int slidecomb_count = 0;
    int reverb1_count = 0;
    int ellipseme_count = 0;
    int feed1_count = 0;
    int flam1_count = 0;
    int comb4_count = 0;
    int ringfeed_count = 0;
    int bendy_count = 0;
    int ringmod_count = 0;
    int ringmod4_count = 0;
    int slideflam_count = 0;
    int bitcrush_count = 0;
    int resonadsr_count = 0;
    int stv_count = 0;
    double raw_wet;
    double sr = x->sr;
    double *dels;
    double **alpo1, **alpo2;
    double *fltdata;
    double xnorm;
    int nsects;
    LSTRUCT *eel1, *eel2;
    double speed1, speed2, mindelay, maxdelay, duration;
    double basefreq;
    double rvt = boundrand(0.1,0.98);
    double lpt;
    double notedur, sust;
    double phz1, phz2, minfeedback = 0.1, maxfeedback = 0.7;
    long slotnum = x->recall_slot;
    double *rparams; // parameter set to recall
    long rpcount; // number of parameters to read
    long slot_pcount;
        
    if(x->recall_parameters_flag == 1 ){
        x->recall_parameters_flag = 0;
    } else {
        return;
    }

    slot_pcount = x->slots[slotnum].pcount;
    
    if( slot_pcount <= 0){
    //    post("Aborting reload of slot %d", slotnum);
        return;
    }
    rparams = x->slots[slotnum].params;
    pcount = 0;
    rpcount = 0;
    // read the pattern here
   // post("loading %d parameters", slot_pcount);
    while( rpcount < slot_pcount ){
        j = rparams[rpcount++];
        if(j == COMB){
            // post("Added COMB unit");
            params[pcount++] = COMB;
            params[pcount++] = delay = rparams[rpcount++];
            params[pcount++] = revtime = rparams[rpcount++];
            params[pcount++] = comb_dl_count = rparams[rpcount++]; // Possible Bug Here???
            mycombset(delay,revtime,0,x->comb_delay_pool1[comb_dl_count],x->sr);
            mycombset(delay,revtime,0,x->comb_delay_pool2[comb_dl_count],x->sr);
        }
        else if(j == RINGMOD) {
            // post("Added RINGMOD unit");
            params[pcount++] = RINGMOD;
            params[pcount++] = rparams[rpcount++]; //need a log version
            params[pcount++] = ringmod_count = rparams[rpcount++];
            x->ringmod_phases[ringmod_count] = 0.0;
        }
        else if(j == RINGMOD4) {
            // post("Added RINGMOD4 unit");
            params[pcount++] = RINGMOD4;
            params[pcount++] = rparams[rpcount++];
            params[pcount++] = ringmod4_count = rparams[rpcount++];
            x->ringmod4_phases[ringmod4_count] = 0.0;
        }
        else if(j == BENDY){
            // post("Added BENDY unit");
            params[pcount++] = BENDY;
            params[pcount++] = bendy_count = rparams[rpcount++];
            params[pcount++] = x->bendy_units[bendy_count].val1 = rparams[rpcount++];
            params[pcount++] = x->bendy_units[bendy_count].val2 = rparams[rpcount++];
            x->bendy_units[bendy_count].counter = 0;
            delset2(x->bendy_units[bendy_count].delayline1, x->bendy_units[bendy_count].dv1, x->bendy_units[bendy_count].val1,x->sr);
            delset2(x->bendy_units[bendy_count].delayline2, x->bendy_units[bendy_count].dv2, x->bendy_units[bendy_count].val2,x->sr);
        }
        else if(j == FLANGE){
            // post("Added FLANGE unit");
            params[pcount++] = FLANGE;
            params[pcount++] = rparams[rpcount++];
            params[pcount++] = rparams[rpcount++];
            params[pcount++] = rparams[rpcount++];
            params[pcount++] = rparams[rpcount++];
            params[pcount++] = flange_count = rparams[rpcount++];
            x->flange_units[flange_count].phase = boundrand(0.0,0.5) * (double)x->sinelen / x->sr; // maybe should have stored this too
        }
        else if(j == BUTTER){
            // post("Added BUTTER unit");
            params[pcount++] = BUTTER;
            params[pcount++] = ftype = rparams[rpcount++];
            params[pcount++] = cf = rparams[rpcount++];
            params[pcount++] = bw = rparams[rpcount++];
            params[pcount++] = butterworth_count = rparams[rpcount++];
            if( ftype == LOPASS) {
                lobut(x->butterworth_units[butterworth_count].data1, cf, sr);
                lobut(x->butterworth_units[butterworth_count].data2, cf, sr);
            } else if (ftype == HIPASS){
                hibut(x->butterworth_units[butterworth_count].data1, cf, sr);
                hibut(x->butterworth_units[butterworth_count].data2, cf, sr);
            }
            else if(ftype == BANDPASS){
                bpbut(x->butterworth_units[butterworth_count].data1, cf, bw, sr);
                bpbut(x->butterworth_units[butterworth_count].data2, cf, bw, sr);
                x->butterworth_units[butterworth_count].bw = bw;
            }
            x->butterworth_units[butterworth_count].cf = cf;
            x->butterworth_units[butterworth_count].ftype = ftype;
        }
        else if(j == TRUNCATE){
            //post("Added TRUNCATE unit");
            params[pcount++] = TRUNCATE;
            params[pcount++] = truncate_count = rparams[rpcount++];
            x->truncate_units[truncate_count].counter = 0;
            x->truncate_units[truncate_count].state = 0;
            x->truncate_units[truncate_count].segsamples = 1;
        }
        else if(j == SWEEPRESON){
            // post("Added SWEEPRESON unit");
            params[pcount++] = SWEEPRESON;
            params[pcount++] = sweepreson_count = rparams[rpcount++];
            params[pcount++] = x->sweepreson_units[sweepreson_count].minfreq = rparams[rpcount++];
            params[pcount++] = x->sweepreson_units[sweepreson_count].maxfreq = rparams[rpcount++];
            params[pcount++] = x->sweepreson_units[sweepreson_count].bwfac = rparams[rpcount++];
            params[pcount++] = x->sweepreson_units[sweepreson_count].speed = rparams[rpcount++];
            params[pcount++] = x->sweepreson_units[sweepreson_count].phase = rparams[rpcount++];
            x->sweepreson_units[sweepreson_count].q1[3] = 0;
            x->sweepreson_units[sweepreson_count].q1[4] = 0;
            x->sweepreson_units[sweepreson_count].q2[3] = 0;
            x->sweepreson_units[sweepreson_count].q2[4] = 0;
        }
        else if(j == SLIDECOMB){
            // post("Added SLIDECOMB unit");
            params[pcount++] = SLIDECOMB;
            params[pcount++] = slidecomb_count = rparams[rpcount++];
            params[pcount++] = x->slidecomb_units[slidecomb_count].start_delay = rparams[rpcount++];
            params[pcount++] = x->slidecomb_units[slidecomb_count].end_delay = rparams[rpcount++];
            params[pcount++] = x->slidecomb_units[slidecomb_count].feedback = rparams[rpcount++];
            params[pcount++] = x->slidecomb_units[slidecomb_count].sample_length = rparams[rpcount++];
            // scale length to SR
            x->slidecomb_units[slidecomb_count].sample_length *= x->sr;
            x->slidecomb_units[slidecomb_count].counter = 0;
            delset2(x->slidecomb_units[slidecomb_count].delayline1,x->slidecomb_units[slidecomb_count].dv1,MAX_SLIDECOMB_DELAY, x->sr);
            delset2(x->slidecomb_units[slidecomb_count].delayline2,x->slidecomb_units[slidecomb_count].dv2,MAX_SLIDECOMB_DELAY, x->sr);
        }
        else if(j == REVERB1){
            // post("Added REVERB1 unit");
            params[pcount++] = REVERB1;
            params[pcount++] = reverb1_count = rparams[rpcount++];
            params[pcount++] = revtime = x->reverb1_units[reverb1_count].revtime = rparams[rpcount++];
            params[pcount++] = raw_wet = rparams[rpcount++];
            x->reverb1_units[reverb1_count].wet = sin(1.570796 * raw_wet);
            x->reverb1_units[reverb1_count].dry = cos(1.570796 * raw_wet);
            dels = x->reverb1_units[reverb1_count].dels;
            for(i = 0; i < 4; i++){
                 params[pcount++] = dels[i] = rparams[rpcount++];
            }
            alpo1 = x->reverb1_units[reverb1_count].alpo1;
            alpo2 = x->reverb1_units[reverb1_count].alpo2;
            for( i = 0; i < 4; i++ ){
                if(dels[i] < .005 || dels[i] > 0.1) {
                    post("reverb1: bad random delay time: %f",dels[i]);
                    dels[i] = .05;
                }
                mycombset(dels[i], revtime, 0, alpo1[i], x->sr);
                mycombset(dels[i], revtime, 0, alpo2[i], x->sr);
            }
            ellipset(x->reverb_ellipse_data,x->reverb1_units[reverb1_count].eel1,&x->reverb1_units[reverb1_count].nsects,&x->reverb1_units[reverb1_count].xnorm);
            ellipset(x->reverb_ellipse_data,x->reverb1_units[reverb1_count].eel2,&x->reverb1_units[reverb1_count].nsects,&x->reverb1_units[reverb1_count].xnorm);
           //  reverb1_count = (reverb1_count + 1) % max_dsp_units;
        }
        else if(j == ELLIPSE){
           // post("Added ELLIPSE unit");
            params[pcount++] = ELLIPSE;
            params[pcount++] = ellipseme_count = rparams[rpcount++];
            params[pcount++] = x->ellipseme_units[ellipseme_count].filtercode = rparams[rpcount++];
            
            if( x->ellipseme_units[ellipseme_count].filtercode >= ELLIPSE_FILTER_COUNT ){
                error("there is no %d ellipse data",x->ellipseme_units[ellipseme_count].filtercode);
                return;
            };
            fltdata = x->ellipse_data [x->ellipseme_units[ellipseme_count].filtercode];
            eel1 = x->ellipseme_units[ellipseme_count].eel1;
            eel2 = x->ellipseme_units[ellipseme_count].eel2;
            ellipset(fltdata,eel1,&nsects,&xnorm);
            ellipset(fltdata,eel2,&nsects,&xnorm);
            x->ellipseme_units[ellipseme_count].nsects = nsects;
            x->ellipseme_units[ellipseme_count].xnorm = xnorm;
            // ellipseme_count = (ellipseme_count + 1) % max_dsp_units;
        }
        else if(j == FEED1){
          //  post("Added FEED1 unit");
            params[pcount++] = FEED1;
            params[pcount++] = feed1_count = rparams[rpcount++];
            params[pcount++] = mindelay = x->feed1_units[feed1_count].mindelay = rparams[rpcount++];
            params[pcount++] = maxdelay = x->feed1_units[feed1_count].maxdelay = rparams[rpcount++];
            params[pcount++] = speed1 = x->feed1_units[feed1_count].speed1 = rparams[rpcount++];
            params[pcount++] = speed2 = x->feed1_units[feed1_count].speed2 = rparams[rpcount++];
            params[pcount++] = duration = x->feed1_units[feed1_count].duration = rparams[rpcount++];
            
            funcgen1(x->feed1_units[feed1_count].func1, FEEDFUNCLEN,duration,
                     mindelay,maxdelay, speed1, speed2, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func2, FEEDFUNCLEN,duration,
                     mindelay * 0.5,maxdelay * 2.0, speed1 * 1.25, speed2 * 0.75, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func3, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.35, speed2*1.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func4, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.55, speed2*2.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            delset2(x->feed1_units[feed1_count].delayLine1a, x->feed1_units[feed1_count].dv1a, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine2a, x->feed1_units[feed1_count].dv2a, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine1b, x->feed1_units[feed1_count].dv1b, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine2b, x->feed1_units[feed1_count].dv2b, MAX_MINI_DELAY, x->sr);
            // feed1_count = (feed1_count + 1) % max_dsp_units;
            
        }
        else if(j == BITCRUSH){
         //   post("Added BITCRUSH unit");
            params[pcount++] = BITCRUSH;
            params[pcount++] = bitcrush_count = rparams[rpcount++];
            params[pcount++] = x->bitcrush_factors[bitcrush_count] = rparams[rpcount++];
            // bitcrush_count = (bitcrush_count + 1) % max_dsp_units;
        }
        else if(j == FLAM1){
            params[pcount++] = FLAM1;
            params[pcount++] = flam1_count = rparams[rpcount++];
            params[pcount++] = x->flam1_units[flam1_count].sample_length = rparams[rpcount++]; // not yet scaled to SR
            params[pcount++] = x->flam1_units[flam1_count].dt = rparams[rpcount++];
            
            x->flam1_units[flam1_count].counter = 0;
            // scale to SR
            x->flam1_units[flam1_count].sample_length *= x->sr;
            
            delset2(x->flam1_units[flam1_count].delayline1, x->flam1_units[flam1_count].dv1, MAX_MINI_DELAY, x->sr);
            delset2(x->flam1_units[flam1_count].delayline2, x->flam1_units[flam1_count].dv2, MAX_MINI_DELAY, x->sr);
           //  flam1_count = (flam1_count + 1) % max_dsp_units;
        }
        else if(j == SLIDEFLAM){
         //   post("Added SLIDEFLAM unit");
            params[pcount++] = SLIDEFLAM;
            params[pcount++] = slideflam_count = rparams[rpcount++];
            params[pcount++] = x->slideflam_units[slideflam_count].dt1 = rparams[rpcount++];
            params[pcount++] = x->slideflam_units[slideflam_count].dt2 = rparams[rpcount++];
            params[pcount++] = x->slideflam_units[slideflam_count].feedback = rparams[rpcount++];
            params[pcount++] = x->slideflam_units[slideflam_count].sample_length = rparams[rpcount++]; // not yet scaled to SR
            
            x->slideflam_units[slideflam_count].counter = 0;
            // scale to SR
            x->slideflam_units[slideflam_count].sample_length *= x->sr;
            
            delset2(x->slideflam_units[slideflam_count].delayline1,x->slideflam_units[slideflam_count].dv1,MAX_SLIDEFLAM_DELAY, x->sr);
            delset2(x->slideflam_units[slideflam_count].delayline2,x->slideflam_units[slideflam_count].dv2,MAX_SLIDEFLAM_DELAY, x->sr);
            // slideflam_count = (slideflam_count + 1) % max_dsp_units;
        }
        else if(j == COMB4){
        //    post("Added COMB4 unit");
            params[pcount++] = COMB4;
            params[pcount++] = comb4_count = rparams[rpcount++];
            rvt = 0.99;
            params[pcount++] = basefreq = rparams[rpcount++];
            // post("comb4: count: %d, basefreq: %f",comb4_count, basefreq);
            for( i = 0; i < 4; i++ ){
                params[pcount++] = lpt = rparams[rpcount++];
                // post("comb4: lpt: %f",lpt);
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs1[i], x->sr);
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs2[i], x->sr);
            }
        }
        else if(j == COMPDIST){
            // post("Added COMPDIST unit");
            params[pcount++] = COMPDIST;
        }
        else if(j == RINGFEED){
            // post("Added RINGFEED unit");
            params[pcount++] = RINGFEED;
            params[pcount++] = ringfeed_count = rparams[rpcount++];
            params[pcount++] = cf = rparams[rpcount++];
            params[pcount++] = bw = rparams[rpcount++];
            
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res1q, x->sr);
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res2q, x->sr);
            x->resonfeed_units[ringfeed_count].osc1phs = 0.0;
            x->resonfeed_units[ringfeed_count].osc2phs = 0.0;
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc1si = rparams[rpcount++]; // not yet scaled to SR
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc2si = rparams[rpcount++]; // not yet scaled to SR
            params[pcount++] = lpt = rparams[rpcount++];
            params[pcount++] = rvt = rparams[rpcount++];
            
            // scale to SR
            x->resonfeed_units[ringfeed_count].osc1si *= ((double)x->sinelen / x->sr);
            x->resonfeed_units[ringfeed_count].osc2si *= ((double)x->sinelen / x->sr);
            
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb1arr, x->sr);
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb2arr, x->sr);
            // is this line superfluous?
            // ringfeed_count = (ringfeed_count + 1) % max_dsp_units;
        }
        else if(j == RESONADSR){
         //   post("Added RESONADSR unit");
            params[pcount++] = RESONADSR;
            params[pcount++] = resonadsr_count = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->a = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->d = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->r = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v1 = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v2 = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v3 = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v4 = rparams[rpcount++];
            params[pcount++] = x->resonadsr_units[resonadsr_count].bwfac = rparams[rpcount++];
            params[pcount++] =  notedur = rparams[rpcount++];
            x->resonadsr_units[resonadsr_count].phs = 0.0;
            sust = notedur - (x->resonadsr_units[resonadsr_count].adsr->a + x->resonadsr_units[resonadsr_count].adsr->d +  x->resonadsr_units[resonadsr_count].adsr->r);
            x->resonadsr_units[resonadsr_count].adsr->s = sust;
            buildadsr(x->resonadsr_units[resonadsr_count].adsr);
            x->resonadsr_units[resonadsr_count].si = ((double)x->resonadsr_units[resonadsr_count].adsr->len / x->sr) / notedur;
            // resonadsr_count = (resonadsr_count + 1) % max_dsp_units;
        }
        else if(j == STV){
         //   post("Added STV unit");
            params[pcount++] = STV;
            params[pcount++] = stv_count = rparams[rpcount++];
            params[pcount++] = speed1 = rparams[rpcount++];
            params[pcount++] = speed2 = rparams[rpcount++];
            params[pcount++] = maxdelay = rparams[rpcount++];

            x->stv_units[stv_count].osc1phs = 0.0;
            x->stv_units[stv_count].osc2phs = 0.0;
            x->stv_units[stv_count].fac2 = 0.5 * (maxdelay - 0.001);
            x->stv_units[stv_count].fac1 = 0.001 + x->stv_units[stv_count].fac2;
            x->stv_units[stv_count].osc1si = ((double)x->sinelen / x->sr) * speed1;
            x->stv_units[stv_count].osc2si = ((double)x->sinelen / x->sr) * speed2;
            delset2(x->stv_units[stv_count].delayline1, x->stv_units[stv_count].dv1, maxdelay, x->sr);
            delset2(x->stv_units[stv_count].delayline2, x->stv_units[stv_count].dv2, maxdelay, x->sr);
            // stv_count = (stv_count + 1) % max_dsp_units;
        }
        else {
            error("el.chameleon~: could not find a process for %d",j);
        }
    }
    
//    events = floor( boundrand( (float)minproc, (float) maxproc) );

    x->pcount = pcount;
}

void chameleon_set_parameters_exec(t_chameleon *x)
{
    long max_dsp_units = x->max_dsp_units;
    float rval;
    int events;
    int i, j;
    int ftype;
    double cf, bw;//, bw;
    float *odds = x->odds;
    int maxproc  = x->max_process_per_note;
    int minproc = x->min_process_per_note;
    double delay, revtime;
    double *params = x->params;
    long pcount = x->pcount;
    int comb_dl_count = 0; // where in the comb pool to grab memory
    int flange_count = 0;
    int truncate_count = 0;
    int butterworth_count = 0;
    int sweepreson_count = 0;
    int slidecomb_count = 0;
    int reverb1_count = 0;
    int ellipseme_count = 0;
    int feed1_count = 0;
    int flam1_count = 0;
    int comb4_count = 0;
    int ringfeed_count = 0;
    int bendy_count = 0;
    int ringmod_count = 0;
    int ringmod4_count = 0;
    int slideflam_count = 0;
    int bitcrush_count = 0;
    int resonadsr_count = 0;
    int stv_count = 0;
    double raw_wet;
    double sr = x->sr;
    double *dels;
    double **alpo1, **alpo2;
    double *fltdata;
    double xnorm;
    int nsects;
    LSTRUCT *eel1, *eel2;
    double speed1, speed2, mindelay, maxdelay, duration;
    double basefreq;
    double *ratios = x->ratios;
    double rvt = boundrand(0.1,0.98);
    double lpt;
    int dex;
    double notedur, sust;
    double phz1, phz2, minfeedback = 0.1, maxfeedback = 0.7;
    
    if(x->set_parameters_flag == 0){
        return;
    } else {
        x->set_parameters_flag = 0;
    }
    if(maxproc <= 0){
        return;
    }
    
    events = floor( boundrand( (float)minproc, (float) maxproc) );

    pcount = 0;
    if( DEBUG_CHAMELEON ){
        post("*** EVENT LIST ***");
    }
    for(i = 0; i < events; i++){
        rval = boundrand(0.0,1.0);
        j = 0;
        while(rval > odds[j]){
            j++;
        }
        if( DEBUG_CHAMELEON ){
            post("event: %d", j);
        }
        if(j == COMB){
            params[pcount++] = COMB;
            params[pcount++] = delay = boundrand(0.01,0.35);
            params[pcount++] = revtime = boundrand(0.5,0.98);
            params[pcount++] = comb_dl_count;
            mycombset(delay,revtime,0,x->comb_delay_pool1[comb_dl_count],x->sr);
            mycombset(delay,revtime,0,x->comb_delay_pool2[comb_dl_count],x->sr);
            ++comb_dl_count;
            if( comb_dl_count >= max_dsp_units){
                comb_dl_count = 0;
            }
        }
        else if(j == RINGMOD) {
            params[pcount++] = RINGMOD;
            params[pcount++] = boundrand(100.0,2000.0); //need a log version
            params[pcount++] = ringmod_count;
            x->ringmod_phases[ringmod_count] = 0.0;
            ringmod4_count = (ringmod_count + 1) % max_dsp_units;
        }
        else if(j == RINGMOD4) {
            params[pcount++] = RINGMOD4;
            params[pcount++] = boundrand(100.0,2000.0); //need a log version
            params[pcount++] = ringmod4_count;
            x->ringmod4_phases[ringmod4_count] = 0.0;
            ringmod4_count = (ringmod4_count + 1) % max_dsp_units;
        }
        else if(j == BENDY){
            params[pcount++] = BENDY;
            params[pcount++] = bendy_count;
            params[pcount++] = x->bendy_units[bendy_count].val1 = boundrand(0.01, BENDY_MAXDEL);
            params[pcount++] = x->bendy_units[bendy_count].val2 = boundrand(0.01, BENDY_MAXDEL);
            x->bendy_units[bendy_count].counter = 0;
            delset3(x->bendy_units[bendy_count].delayline1, x->bendy_units[bendy_count].dv1, x->bendy_units[bendy_count].val1,x->sr, BENDY_MAXDEL);
            delset3(x->bendy_units[bendy_count].delayline2, x->bendy_units[bendy_count].dv2, x->bendy_units[bendy_count].val2,x->sr, BENDY_MAXDEL);
            bendy_count = (bendy_count + 1) % max_dsp_units;
        }
        else if(j == FLANGE){
            params[pcount++] = FLANGE;
            params[pcount++] = boundrand(80.0,400.0);
            params[pcount++] = boundrand(600.0,4000.0);
            params[pcount++] = boundrand(0.1,2.0);
            params[pcount++] = boundrand(0.1,0.95);
            params[pcount++] = flange_count;
            x->flange_units[flange_count].phase = boundrand(0.0,0.5) * (double)x->sinelen / x->sr;
            flange_count = (flange_count + 1) % max_dsp_units;
            /* ++flange_count;
            if( flange_count >= max_dsp_units){
                flange_count = 0;
            }*/
        }
        else if(j == BUTTER){
            params[pcount++] = BUTTER;
            params[pcount++] = ftype = rand() % 3;
            params[pcount++] = cf = boundrand(70.0,3000.0);
            params[pcount++] = bw = cf * boundrand(0.05,0.6);
            params[pcount++] = butterworth_count;
            if( ftype == LOPASS) {
                lobut(x->butterworth_units[butterworth_count].data1, cf, sr);
                lobut(x->butterworth_units[butterworth_count].data2, cf, sr);
            } else if (ftype == HIPASS){
                hibut(x->butterworth_units[butterworth_count].data1, cf, sr);
                hibut(x->butterworth_units[butterworth_count].data2, cf, sr);
            }
            else if(ftype == BANDPASS){
                bpbut(x->butterworth_units[butterworth_count].data1, cf, bw, sr);
                bpbut(x->butterworth_units[butterworth_count].data2, cf, bw, sr);
                x->butterworth_units[butterworth_count].bw = bw;
            }
            x->butterworth_units[butterworth_count].cf = cf;
            x->butterworth_units[butterworth_count].ftype = ftype;
            
            butterworth_count = (butterworth_count + 1) % max_dsp_units;
        }
        else if(j == TRUNCATE){
            params[pcount++] = TRUNCATE;
            params[pcount++] = truncate_count;
            x->truncate_units[truncate_count].counter = 0;
            x->truncate_units[truncate_count].state = 0;
            x->truncate_units[truncate_count].segsamples = 1;
            truncate_count = (truncate_count + 1) % max_dsp_units;
        }
        else if(j == SWEEPRESON){
            params[pcount++] = SWEEPRESON;
            params[pcount++] = sweepreson_count;
            params[pcount++] = x->sweepreson_units[sweepreson_count].minfreq = boundrand(100.0,300.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].maxfreq = boundrand(600.0,6000.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].bwfac = boundrand(0.01,0.2);
            params[pcount++] = x->sweepreson_units[sweepreson_count].speed = boundrand(0.05,2.0);
            params[pcount++] = x->sweepreson_units[sweepreson_count].phase = boundrand(0.0,0.5);
            // scale phase to SR
            x->sweepreson_units[sweepreson_count].phase *= ((double)x->sinelen / x->sr);
            sweepreson_count = (sweepreson_count + 1) % max_dsp_units;
        }
        else if(j == SLIDECOMB){
            params[pcount++] = SLIDECOMB;
            params[pcount++] = slidecomb_count;
            params[pcount++] = x->slidecomb_units[slidecomb_count].start_delay = boundrand(0.001,MAX_SLIDECOMB_DELAY * 0.95);
            params[pcount++] = x->slidecomb_units[slidecomb_count].end_delay = boundrand(0.001,MAX_SLIDECOMB_DELAY * 0.95);
            params[pcount++] = x->slidecomb_units[slidecomb_count].feedback = boundrand(0.7,0.99);
            params[pcount++] = x->slidecomb_units[slidecomb_count].sample_length = boundrand(0.1,2.0);
            // scale length to SR
            x->slidecomb_units[slidecomb_count].sample_length *= x->sr;
            x->slidecomb_units[slidecomb_count].counter = 0;
            delset2(x->slidecomb_units[slidecomb_count].delayline1,x->slidecomb_units[slidecomb_count].dv1,MAX_SLIDECOMB_DELAY, x->sr);
            delset2(x->slidecomb_units[slidecomb_count].delayline2,x->slidecomb_units[slidecomb_count].dv2,MAX_SLIDECOMB_DELAY, x->sr);
            slidecomb_count = (slidecomb_count + 1) % max_dsp_units;
        }
        else if(j == REVERB1){
            params[pcount++] = REVERB1;
            params[pcount++] = reverb1_count;
            params[pcount++] = revtime = x->reverb1_units[reverb1_count].revtime = boundrand(0.25,0.99);
            params[pcount++] = raw_wet = boundrand(0.2,0.8);
            
            x->reverb1_units[reverb1_count].wet = sin(1.570796 * raw_wet);
            x->reverb1_units[reverb1_count].dry = cos(1.570796 * raw_wet);
            dels = x->reverb1_units[reverb1_count].dels;
            
            for(i = 0; i < 4; i++){
                 params[pcount++] = dels[i] = boundrand(.005, .1 );
            }
            
            alpo1 = x->reverb1_units[reverb1_count].alpo1;
            alpo2 = x->reverb1_units[reverb1_count].alpo2;
            for( i = 0; i < 4; i++ ){
                // dels[i] = boundrand(.005, .1 );
                if(dels[i] < .005 || dels[i] > 0.1) {
                    post("reverb1: bad random delay time: %f",dels[i]);
                    dels[i] = .05;
                }
                mycombset(dels[i], revtime, 0, alpo1[i], x->sr);
                mycombset(dels[i], revtime, 0, alpo2[i], x->sr);
            }
            ellipset(x->reverb_ellipse_data,x->reverb1_units[reverb1_count].eel1,&x->reverb1_units[reverb1_count].nsects,&x->reverb1_units[reverb1_count].xnorm);
            ellipset(x->reverb_ellipse_data,x->reverb1_units[reverb1_count].eel2,&x->reverb1_units[reverb1_count].nsects,&x->reverb1_units[reverb1_count].xnorm);
            reverb1_count = (reverb1_count + 1) % max_dsp_units;
        }
        else if(j == ELLIPSE){
            params[pcount++] = ELLIPSE;
            params[pcount++] = ellipseme_count;
            params[pcount++] = x->ellipseme_units[ellipseme_count].filtercode = rand() % ELLIPSE_FILTER_COUNT;
            
            if( x->ellipseme_units[ellipseme_count].filtercode >= ELLIPSE_FILTER_COUNT ){
                error("there is no %d ellipse data",x->ellipseme_units[ellipseme_count].filtercode);
                return;
            };
            fltdata = x->ellipse_data [x->ellipseme_units[ellipseme_count].filtercode];
            eel1 = x->ellipseme_units[ellipseme_count].eel1;
            eel2 = x->ellipseme_units[ellipseme_count].eel2;
            ellipset(fltdata,eel1,&nsects,&xnorm);
            ellipset(fltdata,eel2,&nsects,&xnorm);
            x->ellipseme_units[ellipseme_count].nsects = nsects;
            x->ellipseme_units[ellipseme_count].xnorm = xnorm;
            ellipseme_count = (ellipseme_count + 1) % max_dsp_units;
        }
        else if(j == FEED1){
            params[pcount++] = FEED1;
            params[pcount++] = feed1_count;
            params[pcount++] = mindelay = x->feed1_units[feed1_count].mindelay = boundrand(.001,0.1);
            params[pcount++] = maxdelay = x->feed1_units[feed1_count].maxdelay = boundrand(mindelay,0.1);
            params[pcount++] = speed1 = x->feed1_units[feed1_count].speed1 = boundrand(.01,0.5);
            params[pcount++] = speed2 = x->feed1_units[feed1_count].speed2 = boundrand(speed1,0.5);
            params[pcount++] = duration = x->feed1_units[feed1_count].duration = boundrand(.05,1.0);
            
            funcgen1(x->feed1_units[feed1_count].func1, FEEDFUNCLEN,duration,
                     mindelay,maxdelay, speed1, speed2, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func2, FEEDFUNCLEN,duration,
                     mindelay * 0.5,maxdelay * 2.0, speed1 * 1.25, speed2 * 0.75, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func3, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.35, speed2*1.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            phz1 /= (double) FEEDFUNCLEN; phz2 /= (double) FEEDFUNCLEN;
            funcgen1(x->feed1_units[feed1_count].func4, FEEDFUNCLEN,duration,
                     minfeedback, maxfeedback, speed1*.55, speed2*2.25, 1.0, 1.0,&phz1, &phz2, x->sinewave, x->sinelen);
            delset2(x->feed1_units[feed1_count].delayLine1a, x->feed1_units[feed1_count].dv1a, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine2a, x->feed1_units[feed1_count].dv2a, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine1b, x->feed1_units[feed1_count].dv1b, MAX_MINI_DELAY, x->sr);
            delset2(x->feed1_units[feed1_count].delayLine2b, x->feed1_units[feed1_count].dv2b, MAX_MINI_DELAY, x->sr);
            feed1_count = (feed1_count + 1) % max_dsp_units;
            
        }
        else if(j == BITCRUSH){
            params[pcount++] = BITCRUSH;
            params[pcount++] = bitcrush_count;
            params[pcount++] = x->bitcrush_factors[bitcrush_count] = boundrand(2.0,8.0);
            bitcrush_count = (bitcrush_count + 1) % max_dsp_units;
        }
        else if(j == FLAM1){
            params[pcount++] = FLAM1;
            params[pcount++] = flam1_count;
            params[pcount++] = x->flam1_units[flam1_count].sample_length = boundrand(1.0,4.0); // not yet scaled to SR
            params[pcount++] = x->flam1_units[flam1_count].dt = boundrand(0.08, 0.2);
            
            x->flam1_units[flam1_count].counter = 0;
            // scale to SR
            x->flam1_units[flam1_count].sample_length *= x->sr;
            delset2(x->flam1_units[flam1_count].delayline1, x->flam1_units[flam1_count].dv1, MAX_MINI_DELAY, x->sr);
            delset2(x->flam1_units[flam1_count].delayline2, x->flam1_units[flam1_count].dv2, MAX_MINI_DELAY, x->sr);
            flam1_count = (flam1_count + 1) % max_dsp_units;
        }
        else if(j == SLIDEFLAM){
            params[pcount++] = SLIDEFLAM;
            params[pcount++] = slideflam_count;
            if( boundrand(0.0,1.0) > 0.1 ){
                params[pcount++] = x->slideflam_units[slideflam_count].dt1 = boundrand(0.01,0.05);
                params[pcount++] = x->slideflam_units[slideflam_count].dt2 = boundrand(0.1,MAX_SLIDEFLAM_DELAY * 0.95);
            } else {
                params[pcount++] = x->slideflam_units[slideflam_count].dt1 = boundrand(0.1,MAX_SLIDEFLAM_DELAY * 0.95);
                params[pcount++] = x->slideflam_units[slideflam_count].dt2 = boundrand(0.01,0.05);
            }
            params[pcount++] = x->slideflam_units[slideflam_count].feedback = boundrand(0.7,0.99);
            params[pcount++] = x->slideflam_units[slideflam_count].sample_length = boundrand(0.25,4.0); // not yet scaled to SR
            
            x->slideflam_units[slideflam_count].counter = 0;
            // scale to SR
            x->slideflam_units[slideflam_count].sample_length *= x->sr;
            
            delset2(x->slideflam_units[slideflam_count].delayline1,x->slideflam_units[slideflam_count].dv1,MAX_SLIDEFLAM_DELAY, x->sr);
            delset2(x->slideflam_units[slideflam_count].delayline2,x->slideflam_units[slideflam_count].dv2,MAX_SLIDEFLAM_DELAY, x->sr);
            slideflam_count = (slideflam_count + 1) % max_dsp_units;
        }
        else if(j == COMB4){
            params[pcount++] = COMB4;
            params[pcount++] = comb4_count;
            rvt = 0.99;
            params[pcount++] = basefreq = boundrand(100.0,400.0);
            for( i = 0; i < 4; i++ ){
                dex = (int)floor( boundrand(0.0, 4.0));
                lpt = 1. / basefreq;
                if( dex > 4) { dex = 4; };
                basefreq *= ratios[dex];
                params[pcount++] = lpt;
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs1[i], x->sr);
                mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs2[i], x->sr);
            }
            comb4_count = (comb4_count + 1) % max_dsp_units;
        }
        else if(j == COMPDIST){
            params[pcount++] = COMPDIST;
        }
        else if(j == RINGFEED){
            params[pcount++] = RINGFEED;
            params[pcount++] = ringfeed_count;
            params[pcount++] = cf = boundrand(90.0,1500.0);
            params[pcount++] = bw = boundrand(0.02,0.4) * cf;
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res1q, x->sr);
            rsnset2(cf, bw, RESON_NO_SCL, 0., x->resonfeed_units[ringfeed_count].res2q, x->sr);
            x->resonfeed_units[ringfeed_count].osc1phs = 0.0;
            x->resonfeed_units[ringfeed_count].osc2phs = 0.0;
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc1si = boundrand(90.0,1500.0); // not yet scaled to SR
            params[pcount++] = x->resonfeed_units[ringfeed_count].osc2si = boundrand(90.0,1500.0); // not yet scaled to SR
            params[pcount++] = lpt = 1.0 / boundrand(90.0,1500.0);
            params[pcount++] = rvt = boundrand(.01,.8);
            
            // scale to SR
            x->resonfeed_units[ringfeed_count].osc1si *= ((double)x->sinelen / x->sr);
            x->resonfeed_units[ringfeed_count].osc2si *= ((double)x->sinelen / x->sr);
            
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb1arr, x->sr);
            mycombset(lpt, rvt, 0, x->resonfeed_units[ringfeed_count].comb2arr, x->sr);
            ringfeed_count = (ringfeed_count + 1) % max_dsp_units;
        }
        else if(j == RESONADSR){
            params[pcount++] = RESONADSR;
            params[pcount++] = resonadsr_count;
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->a = boundrand(0.01,0.1);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->d = boundrand(0.01,0.05);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->r = boundrand(0.05,0.5);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v1 = boundrand(150.0,4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v2 = boundrand(150.0,4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v3 = boundrand(150.0,4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v4 = boundrand(150.0,4000.0);
            params[pcount++] = x->resonadsr_units[resonadsr_count].bwfac = boundrand(0.03,0.7);
            params[pcount++] =  notedur = boundrand(0.7, 1.2);
            x->resonadsr_units[resonadsr_count].phs = 0.0;
            sust = notedur - (x->resonadsr_units[resonadsr_count].adsr->a + x->resonadsr_units[resonadsr_count].adsr->d +  x->resonadsr_units[resonadsr_count].adsr->r);
            x->resonadsr_units[resonadsr_count].adsr->s = sust;
            buildadsr(x->resonadsr_units[resonadsr_count].adsr);
            x->resonadsr_units[resonadsr_count].si = ((double)x->resonadsr_units[resonadsr_count].adsr->len / x->sr) / notedur;
            resonadsr_count = (resonadsr_count + 1) % max_dsp_units;
        }
        else if(j == STV){
            params[pcount++] = STV;
            params[pcount++] = stv_count;
            params[pcount++] = speed1 = boundrand(0.025,0.5);
            params[pcount++] = speed2 = boundrand(0.025,0.5);
            params[pcount++] = maxdelay = boundrand(0.001,0.01);
            x->stv_units[stv_count].osc1phs = 0.0;
            x->stv_units[stv_count].osc2phs = 0.0;
            
            x->stv_units[stv_count].fac2 = 0.5 * (maxdelay - 0.001);
            x->stv_units[stv_count].fac1 = 0.001 + x->stv_units[stv_count].fac2;
            x->stv_units[stv_count].osc1si = ((double)x->sinelen / x->sr) * speed1;
            x->stv_units[stv_count].osc2si = ((double)x->sinelen / x->sr) * speed2;
            delset2(x->stv_units[stv_count].delayline1, x->stv_units[stv_count].dv1, maxdelay, x->sr);
            delset2(x->stv_units[stv_count].delayline2, x->stv_units[stv_count].dv2, maxdelay, x->sr);
            stv_count = (stv_count + 1) % max_dsp_units;
        }
        else {
            error("el.chameleon~: could not find a process for %d",j);
        }
    }
    x->pcount = pcount;
}

//void chameleon_dsp64(t_chameleon *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)

void chameleon_dsp(t_chameleon *x, t_signal **sp)
{
    float samplerate;
    x->vs = sp[0]->s_n;
    samplerate = sys_getsr();
    // post("chameleon - dsp blocksize %d, and sample rate %f", x->vs, samplerate);
    // need more comprehensive reinitialization if sr changes
    if(x->sr != samplerate){
        x->sr = samplerate;
        x->delayline1 = (double *)getbytes(((x->maxdelay * x->sr) + 2) * sizeof(double));
        x->delayline2 = (double *)getbytes(((x->maxdelay * x->sr) + 2) * sizeof(double));
    }
    
    //object_method(dsp64, gensym("dsp_add64"),x,chameleon_perform64,0,NULL);
    dsp_add(chameleon_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, (t_int)sp[0]->s_n);
}

// MaxMSP code, not needed
void chameleon_assist (t_chameleon *x, void *b, long msg, long arg, char *dst)
{
    if (msg==1) {
        switch (arg) {
            case 0: sprintf(dst,"(signal) Channel 1 Input"); break;
            case 1: sprintf(dst,"(signal) Channel 2 Input"); break;
        }
    }
    else if (msg==2) {
        switch(arg){
            case 0: sprintf(dst,"(signal) Channel 1 Output"); break;
            case 1: sprintf(dst,"(signal) Channel 2 Output"); break;
            case 2: sprintf(dst,"(list) Data Bank of Stored Settings"); break;
        }
    }
}

/* Code from ellipse.c */

double ellipse(double x, LSTRUCT *eel, int nsects, double xnorm)
{
  register int m;
  double op;
  
  for(m=0;m<nsects;m++) {
    op = x + eel[m].c0 * eel[m].ps0 + eel[m].c2 * eel[m].ps1
      - eel[m].c1 * eel[m].ps2 - eel[m].c3 * eel[m].ps3;
    eel[m].ps1 = eel[m].ps0;
    eel[m].ps0 = x;
    eel[m].ps3 = eel[m].ps2;
    eel[m].ps2 = op;
    x = op;
  }
  return(x*xnorm);
}

void ellipset(double *list, LSTRUCT *eel, int  *nsects, double *xnorm)
{
/* the first argument in the list is the number of sections */
  int m,i;
  *nsects = (int)list[0];
  if(*nsects > MAXSECTS) {
    error("sorry, only configured for %d sections",MAXSECTS);
    return;
  }
  i=1;
  for(m=0;m<*nsects;m++) {
    eel[m].c0 = list[i++];
    eel[m].c1 = list[i++];
    eel[m].c2 = list[i++];
    eel[m].c3 = list[i++];
    eel[m].ps0 = eel[m].ps1 = eel[m].ps2 = eel[m].ps3 = 0;
  }
  *xnorm = list[i];
}
/*set biquad coefficients one time*/
void init_ellipse_data(double **a)
{
    /* 0: hipass at 200 */
        a[0][0] = 4;
        a[0][1] = 1.5156562;
        a[0][2] = -1.9958239;
        a[0][3] = 1;
        a[0][4] = 0.9965234;
        a[0][5] = -1.9996996;
        a[0][6] = 0.97229244;
        a[0][7] = 1;
        a[0][8] = 0.89313463;
        a[0][9] = 1.8828678;
        a[0][10] = -1.9561138;
        a[0][11] = 1;
        a[0][12] = 0.95845959;
        a[0][13] = -1.9999342;
        a[0][14] = -0.065893592;
        a[0][15] = 1;
        a[0][16] = 0.39565826;
        a[0][17] = 0.26225143;
    /* 1: hipass at 500 */
        a[1][0] = 3;
        a[1][1] = -1.9950633;
        a[1][2] = -1.9910746;
        a[1][3] = 1;
        a[1][4] = 0.99729187;
        a[1][5] = -1.9964182;
        a[1][6] = -1.9742933;
        a[1][7] = 1;
        a[1][8] = 0.98250468;
        a[1][9] = -1.999304;
        a[1][10] = -1.8149534;
        a[1][11] = 1;
        a[1][12] = 0.84353385;
        a[1][13] = 0.90419364;
    /* 2: bandpass 280 - 700 */
        a[2][0] = 4;
        a[2][1] = -1.9989934;
        a[2][2] = -1.9946771;
        a[2][3] = 1;
        a[2][4] = 0.99626146;
        a[2][5] = -1.9843098;
        a[2][6] = -1.9807532;
        a[2][7] = 1;
        a[2][8] = 0.99066977;
        a[2][9] = -1.9996779;
        a[2][10] = -1.9795816;
        a[2][11] = 1;
        a[2][12] = 0.9820447;
        a[2][13] = -1.9513627;
        a[2][14] = -1.965153;
        a[2][15] = 1;
        a[2][16] = 0.97142923;
        a[2][17] = 0.013949928;
    /* 3: lopass at 500 */
        a[3][0] = 3;
        a[3][1] = -1.9922858;
        a[3][2] = -1.9903447;
        a[3][3] = 1;
        a[3][4] = 0.99525722;
        a[3][5] = -1.9849712;
        a[3][6] = -1.9765264;
        a[3][7] = 1;
        a[3][8] = 0.97923558;
        a[3][9] = 1;
        a[3][10] = -0.98180316;
        a[3][11] = 0;
        a[3][12] = 0;
        a[3][13] = 0.0014021298;
    /* 4: lopass at 2K */
        a[4][0] = 3;
        a[4][1] = -1.9170388;
        a[4][2] = -1.9264647;
        a[4][3] = 1;
        a[4][4] = 0.99064223;
        a[4][5] = -1.8850187;
        a[4][6] = -1.9092573;
        a[4][7] = 1;
        a[4][8] = 0.95627234;
        a[4][9] = -1.4613313;
        a[4][10] = -1.8821271;
        a[4][11] = 0.99999996;
        a[4][12] = 0.89422169;
        a[4][13] = 0.0071020706;
    /* 5:
       f1,f2,f3= 400.0     900.0     2500.     ripple= 1.000     db= 40.00    */
        a[5][0] = 4;
        a[5][1] = -1.9594111;
        a[5][2] = -1.9900446;
        a[5][3] = 1;
        a[5][4] = 0.9932218;
        a[5][5] = -1.9986938;
        a[5][6] = -1.967975;
        a[5][7] = 1;
        a[5][8] = 0.98456911;
        a[5][9] = -1.8436241;
        a[5][10] = -1.9696535;
        a[5][11] = 1;
        a[5][12] = 0.9745592;
        a[5][13] = -1.9996708;
        a[5][14] = -1.9523041;
        a[5][15] = 1;
        a[5][16] = 0.96284515;
        a[5][17] = 0.0027629927;
    /* 6:  hipass at 500 */
        a[6][0] = 3;
        a[6][1] = -1.9950633;
        a[6][2] = -1.9910746;
        a[6][3] = 1;
        a[6][4] = 0.99729187;
        a[6][5] = -1.9964182;
        a[6][6] = -1.9742933;
        a[6][7] = 1;
        a[6][8] = 0.98250468;
        a[6][9] = -1.999304;
        a[6][10] = -1.8149534;
        a[6][11] = 1;
        a[6][12] = 0.84353385;
        a[6][13] = 0.90419364;
    /* 7: bandpass 1000-4000-6000 */
        a[7][0] = 4;
        a[7][1] = -1.9924766;
        a[7][2] = -1.9573893;
        a[7][3] = 1;
        a[7][4] = 0.97714717;
        a[7][5] = -1.2471186;
        a[7][6] = -1.6082873;
        a[7][7] = 1;
        a[7][8] = 0.91482231;
        a[7][9] = -1.9983103;
        a[7][10] = -1.8512918;
        a[7][11] = 1;
        a[7][12] = 0.88910013;
        a[7][13] = 0.033335148;
        a[7][14] = -1.6413378;
        a[7][15] = 0.99999998;
        a[7][16] = 0.78964041;
        a[7][17] = 0.0087452226;
    /* 8: hipass 4000-10000 */
        a[8][0] = 1;
        a[8][1] = -1.9896868;
        a[8][2] = -1.3953066;
        a[8][3] = 1;
        a[8][4] = 0.58943112;
        a[8][5] = 0.74811328;
     /* 9: bandpass 500-2500-3500 */
        a[9][0] = 6;
        a[9][1] = -1.9975736;
        a[9][2] = -1.9902167;
        a[9][3] = 1;
        a[9][4] = 0.99529287;
        a[9][5] = -1.7460823;
        a[9][6] = -1.853476;
        a[9][7] = 1;
        a[9][8] = 0.97721553;
        a[9][9] = -1.9984481;
        a[9][10] = -1.9737545;
        a[9][11] = 1;
        a[9][12] = 0.98056598;
        a[9][13] = -1.6166383;
        a[9][14] = -1.8408836;
        a[9][15] = 0.99999999;
        a[9][16] = 0.93097271;
        a[9][17] = -1.9997426;
        a[9][18] = -1.9320458;
        a[9][19] = 1;
        a[9][20] = 0.94629262;
        a[9][21] = -0.44018748;
        a[9][22] = -1.8664352;
        a[9][23] = 0.99999993;
        a[9][24] = 0.90871633;
        a[9][25] = 0.00044746789;
    /* 10: bp-300-400-2500-70dB */
        a[10][0] = 8;
        a[10][1] = -1.7823256;
        a[10][2] = -1.9938863;
        a[10][3] = 1;
        a[10][4] = 0.99712611;
        a[10][5] = -1.9981713;
        a[10][6] = -1.8579881;
        a[10][7] = 1;
        a[10][8] = 0.9825214;
        a[10][9] = -1.7151492;
        a[10][10] = -1.9844167;
        a[10][11] = 1;
        a[10][12] = 0.9884184;
        a[10][13] = -1.9986272;
        a[10][14] = -1.8447412;
        a[10][15] = 1;
        a[10][16] = 0.94374559;
        a[10][17] = -1.3382862;
        a[10][18] = -1.9602273;
        a[10][19] = 1;
        a[10][20] = 0.96717992;
        a[10][21] = -1.9994689;
        a[10][22] = -1.8529558;
        a[10][23] = 1;
        a[10][24] = 0.90889168;
        a[10][25] = 1;
        a[10][26] = -1.903171;
        a[10][27] = 0;
        a[10][28] = 0.92280038;
        a[10][29] = -1;
        a[10][30] = 0;
        a[10][31] = 0;
        a[10][32] = 0;
        a[10][33] = 0.00022546378;
}

/* Code from chameleon_helper.c */


void putsine (double *arr, long len)
{
    long i;
    double twopi;
    twopi = 8.0 * atan2(1.,1.);
    for ( i = 0; i < len ; i++) {
        *(arr + i) = sin(twopi * i / len);
    }
}


float boundrand(float min, float max)
{
    return min + (max-min) * ((float)rand()/MY_MAX);
}


void mycombset(double loopt,double rvt,int init,double *a,double srate)
{
    int j;
    
    a[0] =  (3.0 + (loopt * srate + .5));
    a[1] = rvt;
    if(!init) {
        for(j=3; j<(int)*a; j++) {
            a[j] = 0;
        }
        a[2] = 3;
    }
}

double mycomb(double samp,double *a)
{
    double temp,*aptr;
    if ( a[2] >= (int) a[0])
        a[2] = 3;
    aptr = a + (int)a[2];
    a[2]++;
    temp = *aptr;
    *aptr = *aptr * a[1] + samp;
    return(temp);
}

void setweights(float *a, int len)
{
    float sum = 0.0;
    int i;
    for(i=0;i<len;i++)
        sum += a[i];
    if(sum == 0.0){
        error("zero odds sum");
    }
    for(i=0;i<len;i++)
        a[i] /= sum;
    for(i=1;i<len;i++)
        a[i] += a[i-1];
}

// delset2(x->slidecomb_units[i].delayline1, x->slidecomb_units[i].dv1, MAX_SLIDECOMB_DELAY, x->sr);

void  delset2(double *a,int *l,double xmax, double srate)
{
    /* delay initialization.  a is address of float array, l is size-2 int
     * array for bookkeeping variables, xmax, is maximum expected delay */
    
    int i;

    *l = 0;
    // *(l+1) = (int)(xmax * srate + .5);
    // attempted protection:
    *(l+1) = (int)(xmax * srate + .5) - 1;
    for(i = 0; i < *(l+1); i++) *(a+i) = 0;
}

void  delset3(double *a,int *l,double xmax,double srate, double alloc_max)
{
    /* delay initialization.  a is address of float array, l is size-2 int
     * array for bookkeeping variables, xmax, is maximum expected delay */
    
    int i;
    
    if( xmax > alloc_max ) {
        post("chameleon~ delset3 error: %d is too big, compared to %d. Resetting delay length.", xmax, alloc_max );
        xmax = alloc_max;
    }
    l[0] = 0;
    l[1] = (int)(xmax * srate + .5) - 1;
    for(i = 0; i < l[1]; i++) *(a+i) = 0;
}


void delput2(double x,double *a,int *l)
{
    
    /* put value in delay line. See delset. x is float */
    
    if( ((*l) >= 0)  && ((*l) < *(l+1)) ){
        *(a + (*l)++) = x;
    }
    while(*(l) >= *(l+1)){
        *l -= *(l+1);
    }
}

double dliget2(double *a,double wait,int *l,double srate)
{
    /* get interpolated value from delay line, wait seconds old */
    int im1;
    double x = wait * srate;
    int i = x;
    double frac = x - i; // this could be the bug
    i = *l - i;
    im1 = i - 1;
    if(i <= 0) {
        if(i < 0) i += *(l+1);
        if(i < 0) return(0.);
        if(im1 < 0) im1 += *(l+1);
    }
    return(*(a+i) + frac * (*(a+im1) - *(a+i)));
}

void butset(float *a)
{
    a[6] = a[7] = 0.0;
}

void lobut(double *a, double cutoff, double SR)
{
    register float     c;
    
    c = 1.0 / tan( PI * cutoff / SR);
    a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
    a[2] = a[1] + a[1];
    a[3] = a[1];
    a[4] = 2.0 * ( 1.0 - c*c) * a[1];
    a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    a[6] = a[7] = 0.0;
}

void hibut(double *a, double cutoff, double SR)
{
    register float    c;
    
    c = tan( PI * cutoff / SR);
    a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
    a[2] = -2.0 * a[1];
    a[3] = a[1];
    a[4] = 2.0 * ( c*c - 1.0) * a[1];
    a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    a[6] = a[7] = 0.0;
    
}

void bpbut(double *a, double formant, double bandwidth, double SR)
{
    register float  c, d;
    
    c = 1.0 / tan( PI * bandwidth / SR);
    d = 2.0 * cos( 2.0 * PI * formant / SR);
    a[1] = 1.0 / ( 1.0 + c);
    a[2] = 0.0;
    a[3] = -a[1];
    a[4] = - c * d * a[1];
    a[5] = ( c - 1.0) * a[1];
    a[6] = a[7] = 0.0;
}
/* in array can == out array */

void butter_filter(double *buf, double *a, long frames)
{
    
    int i;
    double t,y;
    for(i = 0 ; i < frames; i++) {
        t = *(buf + i) - a[4] * a[6] - a[5] * a[7];
        y = t * a[1] + a[2] * a[6] + a[3] * a[7];
        a[7] = a[6];
        a[6] = t;
        *(buf + i) = y;
    }
}

void rsnset2(double cf,double bw,double scl,double xinit,double *a,double srate)
{
    double c,temp;
    if(!xinit) {
        a[4] = 0;
        a[3] = 0;
    }
    a[2] = exp(-PI2 * bw/srate);
    temp = 1. - a[2];
    c = a[2] + 1;
    a[1] = 4. * a[2]/c * cos(PI2 * cf/srate);
    if(scl < 0) a[0] = 1;
    if(scl) a[0] = sqrt(temp/c*(c*c-a[1]*a[1]));
    if(!scl) a[0] = temp*sqrt(1.-a[1]*a[1]/(4.*a[2]));
}

double reson(double x,double *a)
{
    double temp;
    temp = *a * x + *(a+1) * *(a+3) - *(a+2) * *(a+4);
    *(a+4) = *(a+3);
    *(a+3) = temp;
    return(temp);
}

double allpass(double samp,double *a)
{
    double temp,*aptr;
    if ( a[STARTM1] >= (int) a[0]) a[STARTM1] = START;
    aptr = a + (int)a[STARTM1];
    a[STARTM1] ++;
    temp = *aptr;
    *aptr = *aptr * a[1] + samp;
    return(temp - a[1] * *aptr);
}

void init_reverb_data(double *a)
{
    a[0] = 2;
    a[1] = -0.61043329;
    a[2] = -1.4582246;
    a[3] = 1;
    a[4] = 0.75887003;
    a[5] = 1;
    a[6] = -0.6922953;
    a[7] = 0;
    a[8] = 0;
    a[9] = 0.035888535;
}

void setflamfunc1(double *arr, int flen)
{
    int i;
    double x;
    for ( i = 0; i < flen; i++){
        x = (double)i / (double) flen ;
        *(arr + i) = ((x - 1) / (x + 1)) * -1.  ;
        
    }
}


void setExpFlamFunc(float *arr, int flen, float v1,float v2,float alpha)
{
    int i;
//    double exp();
    
    if( alpha == 0 )
        alpha = .00000001 ;
    
    for ( i = 0; i < flen; i++){
        *(arr + i) = v1 + (v2-v1) * ((1-exp((float)i*alpha/((float)flen-1.)))/(1-exp(alpha)));
    }
}

void funcgen1(double *outArray, int outlen, double duration, double outMin, double outMax,
              double speed1, double speed2, double gain1, double gain2, double *phs1, double *phs2,
              double *sine, int sinelen)
{
    double si1, si2;
    double localSR;
    int i;
    
    // distrust all of these parameters!
    localSR = duration * (double) outlen ;
    *phs1 *= (double) sinelen;
    *phs2 *= (double) sinelen;
    si1 = ((double)sinelen/localSR)  * speed1;
    si2 = ((double)sinelen/localSR)  * speed2;
    
    // crash situation  - maybe outArray is too short...
    for( i = 0; i < outlen; i++ ){
        *(outArray + i) = oscil(gain1, si1, sine, sinelen, phs1) ; // this is a crasher
        *(outArray + i) += oscil(gain2, si2, sine, sinelen, phs2) ;
    }
    normtab( outArray, outArray, outMin, outMax, outlen);
}

void normtab(double *inarr,double *outarr, double min, double max, int len)
{
    int i;
    
    float imin=9999999999., imax=-9999999999.;
    
    for(i = 0; i < len ; i++){
        if( imin > inarr[i] )
            imin = inarr[i];
        if( imax < inarr[i] )
            imax = inarr[i];
    }
    for(i = 0; i < len; i++ )
        outarr[i] = mapp(inarr[i], imin, imax, min, max);
    
}

double mapp(double in,double imin,double imax,double omin,double omax)
{
    if( imax == 0.0 )
    {
        return 0.0 ;
    } else if(imax == imin){
        return imin;
    }
    return( omin+((omax-omin)*((in-imin)/(imax-imin))) );
}

double oscil(double amp,double si,double *farray,int len,double *phs)
{
    register int i =  *phs;
    *phs += si;
    while(*phs >= len){
        *phs -= len;
    }
    while(*phs < 0.){
        *phs += len;
    }
    if(i < 0){
        i = 0;
    } else if(i >= len){
        i = len - 1;
    }
    return(*(farray+i) * amp);
}

void set_dcflt(double *a)
{
    a[0] = 3;
    a[1] = -1.9999924;
    a[2] = -1.9992482;
    a[3] = 1;
    a[4] = 0.99928019;
    a[5] = -1.9999956;
    a[6] = -1.996408;
    a[7] = 1;
    a[8] = 0.99645999;
    a[9] = -1.9999994;
    a[10] = -1.9805074;
    a[11] = 1;
    a[12] = 0.98069401;
    a[13] = 0.98817413;
}

void set_distortion_table(double *arr, double cut, double max, int len)
{
    int i, j, len2;
    double samp;
    
    len2 = len / 2; // len>>1 ;
    // clean array;
    
    for(i = 0; i < len; i++){
        arr[i] = 0.0;
    }
    // first half
    for( i = len2; i < len; i++ ){
        samp = (double)(i - len2) / (double) len2 ;
        if( samp > cut ){
            samp = mapp( samp, cut, 1.0,  cut, max );
        }
        arr[i] = samp;
    }
    for( i = len2 - 1, j = len2; i > 0; i--, j++ )
        arr[i] = arr[j] * -1.0;
}

double dlookup(double samp,double *arr,int len)
{
    double raw = (samp + 1.0) / 2.0;
    double scaled = raw * (double) len;
    int index = floor(scaled);
    if( index > len - 1) {
        index = len - 1;
    } else if( index < 0 ){
        index = 0;
    }
    return arr[index];
}

void do_compdist(double *in,double *out,int sampFrames,int nchans,int channel,
                 double cutoff,double maxmult,int lookupflag,double *table,int range,double bufMaxamp)
{
    
    int i;
    
    float rectsamp;
    
    for( i = channel ; i < sampFrames * nchans; i+= nchans )
    {
        
        if( lookupflag){
            *(out + i) = dlookup( *(in + i)/bufMaxamp, table, range );
        } else {
            rectsamp = fabs( *(in + i) ) / bufMaxamp;
            if( rectsamp > cutoff ){
                *(in + i) = *(out + i) *
                mapp( rectsamp, cutoff, 1.0, cutoff, maxmult);
            }
        }
    }
}

double getmaxamp(double *arr, int len)
{
    int i;
    double max = 0;
    
    for(i = 0; i < len; i++ ){
        if( fabs(arr[i]) > max )
            max = fabs(arr[i]);
    }
    return max;
}

void buildadsr(CMIXADSR *a)
{
    double A = a->a;
    double D = a->d;
    double S = a->s;
    double R = a->r;
    double f1 = a->v1;
    double f2 = a->v2;
    double f3 = a->v3;
    double f4 = a->v4;
    
    int funclen = a->len;
    double *func = a->func;
    double total;
    int ipoint = 0;
    int i;
    int segs[4];
    double m1,m2;
    total = A + D + S + R ;
    
    segs[0] = (A/total) * funclen;
    segs[1] = (D/total) * funclen;
    segs[2] = (S/total) * funclen;
    segs[3] = funclen - (segs[0]+segs[1]+segs[2]);
    
    if( f1 > 20000. || f1 < -20000. ){
        f1 = 250.0;
    }
    if( f2 > 20000. || f2 < -20000. ){
        f2 = 1250.0;
    }
    if( f3 > 20000. || f3 < -20000. ){
        f3 = 950.0;
    }
    if( f4 > 20000. || f4 < -20000. ){
        f4 = f1;
    }
    
    if( segs[0] <= 0 || segs[1] <= 0 || segs[2] <= 0 || segs[3] <= 0 ){
        
        for( i = 0; i < 4; i++ ){
            segs[i] = funclen / 4;
        }
    }
    
    for( i = 0 ; i < segs[0]; i++ ){
        m1 = 1.-(float)i/(float)(segs[0]);
        m2 = 1. - m1;
        *(func +i ) = f1 * m1 + f2 * m2;
    }
    ipoint = i;
    
    for( i = 0 ; i < segs[1]; i++ ){
        m1 = 1.-(float)i/(float)(segs[1]);
        m2 = 1. - m1;
        *(func + i + ipoint) = f2 * m1 + f3 * m2;
    }
    ipoint += i;
    
    for( i = 0 ; i < segs[2]; i++ ){
        m1 = 1.-(float)i/(float)(segs[2]);
        m2 = 1. - m1;
        *(func + i + ipoint) = f3;
    }
    ipoint += i;
    
    for( i = 0 ; i < segs[3]; i++ ){
        m1 = 1.-(float)i/(float)(segs[3]);
        m2 = 1. - m1;
        *(func + ipoint + i) = f3 * m1 + f4 * m2;
    }
    ipoint += i;
    
}
/* Code from chameleon_dsp.c */

void ringmod(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *sinewave = x->sinewave;
    int sinelen = x->sinelen;
    double *params = x->params;
    float srate = x->sr;
    long vs = x->vs;
    int i;
    double phase = 0.0;
    long iphase;
    double si;
    double rmodFreq;
    int ringmod_count;
    
    ++(*pcount);
    rmodFreq = params[(*pcount)++];
    ringmod_count = params[(*pcount)++];
    phase = x->ringmod_phases[ringmod_count];
    
    si = ((double) sinelen / srate) * rmodFreq;
    
    for(i = 0; i < vs; i++ ){
        iphase = (long) phase;
        buf1[i] = buf1[i] * sinewave[iphase];
        buf2[i] = buf2[i] * sinewave[iphase];
        phase += si;
        while( phase > sinelen ){
            phase -= sinelen;
        }
    }
    x->ringmod_phases[ringmod_count] = phase;
}

void ringmod4(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *sinewave = x->sinewave;
    int sinelen = x->sinelen;
    double *params = x->params;
    float srate = x->sr;
    long vs = x->vs;
    int i;
    double phase = 0.0;
    long iphase;
    double si;
    double rmodFreq;
    int ringmod_count;
    
    ++(*pcount);
    rmodFreq = params[(*pcount)++];
    ringmod_count =params[(*pcount)++];
    phase = x->ringmod4_phases[ringmod_count];
    
    si = ((double) sinelen / srate) * rmodFreq;
    
    // ((a*a *b) - (a*b*b))
    for(i = 0; i < vs; i++ ){
        iphase = (long) phase;
        buf1[i] = (buf1[i] * buf1[i] * sinewave[iphase]) - (buf1[i] * sinewave[iphase] * sinewave[iphase]);
        buf2[i] = (buf2[i] * buf2[i] * sinewave[iphase]) - (buf2[i] * sinewave[iphase] * sinewave[iphase]);
        phase += si;
        while( phase > sinelen ){
            phase -= sinelen;
        }
    }
    x->ringmod_phases[ringmod_count] = phase;
}

void comber(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *delayline1;
    double *delayline2;
    double delay, revtime;
    long vs = x->vs;
    int i;
    int comb_dl_count;
    double *params = x->params;
    /******************************/
    ++(*pcount);
    delay = params[(*pcount)++];
    revtime = params[(*pcount)++];
    comb_dl_count = params[(*pcount)++];
    
    delayline1 = x->comb_delay_pool1[comb_dl_count];
    delayline2 = x->comb_delay_pool2[comb_dl_count];
    // ADD IN ORIGINAL SIGNAL
    for( i = 0; i < vs; i++){
        buf1[i] += mycomb(buf1[i], delayline1);
        buf2[i] += mycomb(buf2[i], delayline2);
    }
}

void flange(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    double si;
    double mindel, maxdel;
    double fac1, fac2;
    double delsamp1, delsamp2;
    double delay_time;
    double speed, feedback, minres, maxres;
    double *params = x->params;
    double srate = x->sr;
    long vs = x->vs;
    double *flange_dl1;
    double *flange_dl2;
    double max_delay = x->max_flangedelay;
    double *sinewave = x->sinewave;
    int sinelen = x->sinelen;
    int *dv1 = x->flange_dv1;
    int *dv2 = x->flange_dv2;
    double phase = x->flange_phase;
    int flange_count;
    
    ++(*pcount);
    minres = params[(*pcount)++];
    maxres = params[(*pcount)++];
    speed = params[(*pcount)++];
    feedback = params[(*pcount)++];
    flange_count = params[(*pcount)++];
    
    flange_dl1 = x->flange_units[flange_count].flange_dl1;
    flange_dl2 = x->flange_units[flange_count].flange_dl2;
    dv1 = x->flange_units[flange_count].dv1;
    dv2 = x->flange_units[flange_count].dv2;
    phase = x->flange_units[flange_count].phase;
    
    if( minres <= 0. || maxres <= 0. ){
        error("flange: got zero frequency resonances as input");
        return;
    }
    mindel = 1.0/maxres;
    maxdel = 1.0/minres;
    // added safety
    if( maxdel > max_delay * 0.99 ){
        maxdel = max_delay * 0.99;
        error("flange: excessive delay time shortened");
    }
    
    si = ((double)sinelen/srate) * speed;
    delsamp1 = delsamp2 = 0;
    fac2 = .5 * (maxdel - mindel);
    fac1 = mindel + fac2;
    
    for(i = 0; i < vs; i++ ){
        /* homemade oscillator */
        delay_time = fac1 + fac2 *  sinewave[(int) phase];
        if( delay_time < .00001 ){
            delay_time = .00001;
        } else if(delay_time >= maxdel) {
            delay_time = maxdel;
        }
        phase += si;
        while( phase > sinelen ){
            phase -= sinelen;
        }
        delput2(buf1[i] + delsamp1 * feedback, flange_dl1, dv1);
        delsamp1 = dliget2(flange_dl1, delay_time, dv1, srate);
        buf1[i] += delsamp1;
 
        delput2(buf2[i] + delsamp2 * feedback, flange_dl2, dv2);
        delsamp2 = dliget2(flange_dl2, delay_time, dv2, srate);
        buf2[i] += delsamp2;
    }
    x->flange_units[flange_count].phase = phase;
}

void butterme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *params = x->params;
    int butterworth_count;
    double *data1, *data2;
    long vs = x->vs;
    int ftype;
    double bw, cf;
    
    ++(*pcount);
    ftype = params[(*pcount)++];
    cf = params[(*pcount)++];
    bw = params[(*pcount)++];
    butterworth_count = params[(*pcount)++];
    
    data1 = x->butterworth_units[butterworth_count].data1;
    data2 = x->butterworth_units[butterworth_count].data2;

    butter_filter(buf1, data1, vs);
    butter_filter(buf2, data2, vs);
}



void bitcrush(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *params = x->params;
    int bitcrush_count;
    double bitcrush_factor;
    long vs = x->vs;
    int i;
    ++(*pcount);
    bitcrush_count = params[(*pcount)++];
    bitcrush_factor = params[(*pcount)++];

    bitcrush_factor = x->bitcrush_factors[bitcrush_count];
    for(i = 0; i < vs; i++){
        buf1[i] = floor(buf1[i] * bitcrush_factor) / bitcrush_factor;
        buf2[i] = floor(buf2[i] * bitcrush_factor) / bitcrush_factor;
    }
}


void truncateme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double *params = x->params;
    double env_gain;
    int truncate_count;
    long vs = x->vs;
    int i;
    long counter;
    long state;
    long segsamples;
    double samp1, samp2;
    ++(*pcount);
    truncate_count = params[(*pcount)++];

    counter = x->truncate_units[truncate_count].counter;
    state = x->truncate_units[truncate_count].state;
    segsamples = x->truncate_units[truncate_count].segsamples;

    for(i = 0; i < vs; i++){
        samp1 = buf1[i];
        samp2 = buf2[i];
        counter++;
        if( (state == 0) || (state == 1) || (state == 2) ){
            buf1[i] = 0.05; // changed so not totally silent
            buf2[i] = 0.05;
        }
        else if( state == 3 ){
            env_gain = (double)counter / (double)segsamples;
            buf1[i] = samp1 * env_gain;
            buf2[i] = samp2 * env_gain;
        }
        /*
        else if( (state == 4) || (state == 5) || (state == 6) ) {
         just output unprocessed sound
        }
        */
        else if( state == 7 ){
            env_gain = 1.0 - ((double)counter / (double)segsamples);
            buf1[i] = samp1 * env_gain;
            buf2[i] = samp2 * env_gain;
        }
        
        if( counter >= segsamples ){
            if( state == 0){
                segsamples = x->sr * boundrand(0.1, 2.0); // silence
            }
            else if( (state == 2) || (state == 6) ){
                segsamples = x->sr * 0.05; // transition
            }
            else if(state == 4){
                segsamples = x->sr * boundrand(1.0, 7.0); // sustain
            }
            else if( (state == 1 || (state == 3) || (state == 5) || (state == 7)) ){
                segsamples = 1; // single sample transition to next state
            }
            counter = 0;
            state = (state + 1) % 8;
        }
    }
    
    x->truncate_units[truncate_count].state = state;
    x->truncate_units[truncate_count].segsamples = segsamples;
    x->truncate_units[truncate_count].counter = counter; // update count
}

void sweepreson(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    double *params = x->params;
    double bwfac;
    double minfreq, maxfreq, speed, phase;
    double *q1, *q2;
    double cf, bw;
    double si;
    double fac1, fac2;
    double srate = x->sr;
    long vs = x->vs;
    double *sinewave = x->sinewave;
    int sinelen = x->sinelen ;
    long sweepreson_count;

    ++(*pcount);
    sweepreson_count = params[(*pcount)++]; // note swapped location
    minfreq = params[(*pcount)++];
    maxfreq = params[(*pcount)++];
    bwfac = params[(*pcount)++];
    speed = params[(*pcount)++];
    phase = params[(*pcount)++]; // cannot use this, phase is reentrant

    minfreq = x->sweepreson_units[sweepreson_count].minfreq;
    maxfreq = x->sweepreson_units[sweepreson_count].maxfreq;
    bwfac = x->sweepreson_units[sweepreson_count].bwfac;
    speed = x->sweepreson_units[sweepreson_count].speed;
    phase = x->sweepreson_units[sweepreson_count].phase;
    q1 = x->sweepreson_units[sweepreson_count].q1;
    q2 = x->sweepreson_units[sweepreson_count].q2;
    
    si = ((double) sinelen / srate) * speed;
    fac2 = .5 * (maxfreq - minfreq) ;
    fac1 = minfreq + fac2;
    cf = fac1 + fac2 * sinewave[(int) phase];
    bw = bwfac * cf;
    
    for(i = 0; i < vs; i++ ){
        phase += si;
        while( phase >= sinelen ){
            phase -= sinelen;
        }
        fac2 = .5 * (maxfreq - minfreq) ;
        fac1 = minfreq + fac2;
        if( phase < 0.0 || phase > sinelen - 1){
            phase = 0;
            // post("sweepreson - sine phase out of bounds");
        }
        cf = fac1 + fac2 * sinewave[(int) phase];
        bw = bwfac * cf;
        if(cf < 10 || cf > 8000 || bw < 1 || srate < 100){
            post("sweepreson - danger values, cf %f bw %f sr %f",cf, bw, srate);
        }
        
        // void rsnset2(double cf,double bw,double scl,double xinit,double *a,double srate)
        rsnset2( cf, bw, 2.0, 1.0, q1, srate );
        rsnset2( cf, bw, 2.0, 1.0, q2, srate );
        buf1[i] = reson(buf1[i], q1);
        buf2[i] = reson(buf2[i], q2);
    }
    x->sweepreson_units[sweepreson_count].phase = phase;
}

void slidecomb(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double feedback, delay1, delay2;
    int i;
    int *dv1, *dv2;        /* cmix bookkeeping */
    double delsamp1 = 0, delsamp2 = 0;
    double m1, m2;
    double delay_time;
    double *params = x->params;
    double srate = x->sr;
    double *delayline1;
    double *delayline2;
    int slidecomb_count;
    long vs = x->vs;
    long sample_length;
    long counter;
    
    ++(*pcount);
    slidecomb_count = params[(*pcount)++];
    delay1 = params[(*pcount)++];
    delay2 = params[(*pcount)++];
    feedback = params[(*pcount)++];
    sample_length = params[(*pcount)++];

    delay1 = x->slidecomb_units[slidecomb_count].start_delay;
    delay2 = x->slidecomb_units[slidecomb_count].end_delay;
    sample_length = x->slidecomb_units[slidecomb_count].sample_length;
    counter = x->slidecomb_units[slidecomb_count].counter;
    feedback = x->slidecomb_units[slidecomb_count].feedback;
    dv1 = x->slidecomb_units[slidecomb_count].dv1;
    dv2 = x->slidecomb_units[slidecomb_count].dv2;
    delayline1 = x->slidecomb_units[slidecomb_count].delayline1;
    delayline2 = x->slidecomb_units[slidecomb_count].delayline2;

    for( i = 0; i < vs; i++){
        if(counter < sample_length){
            m2 = (double)counter / (double)sample_length;
            m1 = 1. - m2;
            delay_time = delay1 * m1 + delay2 * m2;
            counter++;
        } else {
            // implement switchback here
            double tmp;
            // post("swapping out slidecomb");
            counter = 0;
            tmp = x->slidecomb_units[slidecomb_count].start_delay;
            delay1 = x->slidecomb_units[slidecomb_count].start_delay = x->slidecomb_units[slidecomb_count].end_delay;
            delay2 = x->slidecomb_units[slidecomb_count].end_delay = tmp;
            // delay_time = delay2;
        }

        if(delay_time > MAX_SLIDECOMB_DELAY){
            delay_time = MAX_SLIDECOMB_DELAY;
        } else if(delay_time < 0.0){
            delay_time = 0.0;
        }
        
        delput2(buf1[i] + (delsamp1 * feedback), delayline1, dv1);
        delsamp1 = dliget2(delayline1, delay_time, dv1, srate);
        buf1[i] = delsamp1;
        
        delput2(buf2[i] + delsamp2*feedback, delayline2, dv2);
        delsamp2 = dliget2(delayline2, delay_time, dv2, srate);
        buf2[i] = delsamp2;
    }
    x->slidecomb_units[slidecomb_count].dv1 = dv1;
    x->slidecomb_units[slidecomb_count].dv2 = dv2;
    x->slidecomb_units[slidecomb_count].counter = counter;
}

void slideflam(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    double feedback, delay1, delay2;
    int i;
    int *dv1, *dv2;        /* cmix bookkeeping */
    double delsamp1 = 0, delsamp2 = 0;
    double m1, m2;
    double delay_time;
    double *params = x->params;
    double srate = x->sr;
    double *delayline1;
    double *delayline2;
    int slideflam_count;
    long vs = x->vs;
    long sample_length;
    long counter;
    
    ++(*pcount);
    slideflam_count = params[(*pcount)++];
    delay1 = params[(*pcount)++];
    delay2 = params[(*pcount)++];
    feedback = params[(*pcount)++];
    sample_length = params[(*pcount)++];
 /*
  params[pcount++] = SLIDEFLAM;
  params[pcount++] = slideflam_count;
  if( boundrand(0.0,1.0) > 0.1 ){
      params[pcount++] = x->slideflam_units[slideflam_count].dt1 = boundrand(0.01,0.05);
      params[pcount++] = x->slideflam_units[slideflam_count].dt2 = boundrand(0.1,MAX_SLIDEFLAM_DELAY * 0.95);
  } else {
      params[pcount++] = x->slideflam_units[slideflam_count].dt1 = boundrand(0.1,MAX_SLIDEFLAM_DELAY * 0.95);
      params[pcount++] = x->slideflam_units[slideflam_count].dt2 = boundrand(0.01,0.05);
  }
  params[pcount++] = x->slideflam_units[slideflam_count].feedback = boundrand(0.7,0.99);
  params[pcount++] = x->slideflam_units[slideflam_count].sample_length = boundrand(0.25,4.0) * x->sr;
  */
    delay1 = x->slideflam_units[slideflam_count].dt1;
    delay2 = x->slideflam_units[slideflam_count].dt2;
    sample_length = x->slideflam_units[slideflam_count].sample_length;
    counter = x->slideflam_units[slideflam_count].counter;
    feedback = x->slideflam_units[slideflam_count].feedback;
    dv1 = x->slideflam_units[slideflam_count].dv1;
    dv2 = x->slideflam_units[slideflam_count].dv2;
    delayline1 = x->slideflam_units[slideflam_count].delayline1;
    delayline2 = x->slideflam_units[slideflam_count].delayline2;
    for( i = 0; i < vs; i++){
        if(counter < sample_length){
            m2 = (double)counter / (double)sample_length;
            m1 = 1. - m2;
            delay_time = delay1 * m1 + delay2 * m2;
            counter++;
        } else {
            double tmp;
            counter = 0;
            // delay_time = delay2;
            counter = 0;
            tmp = x->slideflam_units[slideflam_count].dt1;
            delay1 = x->slideflam_units[slideflam_count].dt1 = x->slideflam_units[slideflam_count].dt2;
            delay2 = x->slideflam_units[slideflam_count].dt2 = tmp;
        }
        if( (delay_time < 0.0) || (delay_time >= MAX_SLIDEFLAM_DELAY) ){
            if(DEBUG_CHAMELEON){
                post("chameleon~: bad delay time: %f\n", delay_time);
            }
            delay_time = MAX_SLIDEFLAM_DELAY;
        }
        delput2(buf1[i] + (delsamp1 * feedback), delayline1, dv1);
        delsamp1 = dliget2(delayline1, delay_time, dv1, srate);
        buf1[i] = delsamp1;
        
        delput2(buf2[i] + delsamp2*feedback, delayline2, dv2);
        delsamp2 = dliget2(delayline2, delay_time, dv2, srate);
        buf2[i] = delsamp2;
    }
    x->slideflam_units[slideflam_count].dv1 = dv1;
    x->slideflam_units[slideflam_count].dv2 = dv2;
    x->slideflam_units[slideflam_count].counter = counter;
}


void bendy(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    // double feedback, delay1, delay2;
    int i;
    int *dv1, *dv2;        /* cmix bookkeeping */
    double delsamp1 = 0, delsamp2 = 0;
    double delay_time;
    double *params = x->params;
    double srate = x->sr;
    double *delayline1;
    double *delayline2;
    int bendy_count;
    long vs = x->vs;
    long counter;
    long segment_samples;
    double frak;
    double val1, val2;
    double segdur;
    ++(*pcount);
    bendy_count = params[(*pcount)++];
    val1 = params[(*pcount)++];
    val2 = params[(*pcount)++];
    
    counter = x->bendy_units[bendy_count].counter;
    val1 = x->bendy_units[bendy_count].val1;
    val2 = x->bendy_units[bendy_count].val2;
    dv1 = x->bendy_units[bendy_count].dv1;
    dv2 = x->bendy_units[bendy_count].dv2;
    delayline1 = x->bendy_units[bendy_count].delayline1;
    delayline2 = x->bendy_units[bendy_count].delayline2;
    segment_samples = x->bendy_units[bendy_count].segment_samples;
    
    for( i = 0; i < vs; i++){
        if(counter <= 0){
            val1 = val2;
            val2 = boundrand(0.01,BENDY_MAXDEL);
            segdur = boundrand(BENDY_MINSEG, BENDY_MAXSEG);
            counter = (long) (segdur * srate);
            segment_samples = counter;
            // post("v1 %f v2 %f segsamps %d\n", val1,val2, segment_samples);
        }
        frak = 1.0 - ((double)counter / (double) segment_samples);
        delay_time = val2 + (frak * (val1 - val2));
        if(delay_time > BENDY_MAXDEL){
            delay_time = BENDY_MAXDEL;
        } else if(delay_time < 0.0){
            delay_time = 0.0;
        }
        delput2(buf1[i], delayline1, dv1); // this is a crasher
        delsamp1 = dliget2(delayline1, delay_time, dv1, srate);
        buf1[i] = delsamp1;
        
        delput2(buf2[i], delayline2, dv2);
        delsamp2 = dliget2(delayline2, delay_time, dv2, srate);
        buf2[i] = delsamp2;
        counter--;
    }
    x->bendy_units[bendy_count].counter = counter;
    x->bendy_units[bendy_count].val1 = val1;
    x->bendy_units[bendy_count].val2 = val2;
    x->bendy_units[bendy_count].dv1 = dv1;
    x->bendy_units[bendy_count].dv2 = dv2;
    x->bendy_units[bendy_count].segment_samples = segment_samples;
}


void reverb1(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    float revtime;
    double *params = x->params;
    long vs = x->vs;
    int reverb1_count;
    double wet, dry, raw_wet;
    double a1,a2,a3,a4;
    double **alpo1, **alpo2;
    int nsects;
    double xnorm;
    LSTRUCT *eel1, *eel2;
    double *dels;
    int i;
    
    ++(*pcount);
    reverb1_count = params[(*pcount)++];
    revtime = params[(*pcount)++];
    raw_wet = params[(*pcount)++];
    wet = sin(1.570796 * raw_wet);
    dry = cos(1.570796 * raw_wet);
    
    dels = x->reverb1_units[reverb1_count].dels;
    for(i = 0; i < 4; i++){
         dels[i] = params[(*pcount)++];
    }
    
    wet = x->reverb1_units[reverb1_count].wet;
    dry = x->reverb1_units[reverb1_count].dry;
    revtime = x->reverb1_units[reverb1_count].revtime;
    alpo1 = x->reverb1_units[reverb1_count].alpo1;
    alpo2 = x->reverb1_units[reverb1_count].alpo2;
    nsects = x->reverb1_units[reverb1_count].nsects;
    xnorm = x->reverb1_units[reverb1_count].xnorm;
    eel1 = x->reverb1_units[reverb1_count].eel1;
    eel2 = x->reverb1_units[reverb1_count].eel2;
    
    if( revtime >= 1. ){
        error("reverb1 does not like feedback values over 1.");
        revtime = .99 ;
    }

    for( i = 0 ; i < vs; i++){
        a1 = allpass(buf1[i], alpo1[0]);
        a2 = allpass(buf1[i], alpo1[1]);
        a3 = allpass(buf1[i], alpo1[2]);
        a4 = allpass(buf1[i], alpo1[3]);
        buf1[i] = (buf1[i] * dry) + (ellipse((a1+a2+a3+a4), eel1, nsects,xnorm) * wet);
        
        
        a1 = allpass(buf2[i], alpo2[0]);
        a2 = allpass(buf2[i], alpo2[1]);
        a3 = allpass(buf2[i], alpo2[2]);
        a4 = allpass(buf2[i], alpo2[3]);
        buf2[i] = (buf2[i] * dry) + (ellipse((a1+a2+a3+a4), eel2, nsects,xnorm) * wet);
    }
}

void ellipseme(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    int nsects;
    double xnorm;
    int ellipseme_count;
    double *params = x->params;
    long vs = x->vs;
    int filtercode;
    LSTRUCT *eel1, *eel2;
    
    ++(*pcount);
    ellipseme_count = params[(*pcount)++];
    filtercode = params[(*pcount)++];

    eel1 = x->ellipseme_units[ellipseme_count].eel1;
    eel2 = x->ellipseme_units[ellipseme_count].eel2;
    nsects = x->ellipseme_units[ellipseme_count].nsects;
    xnorm = x->ellipseme_units[ellipseme_count].xnorm;

    for( i = 0; i < vs; i++){
        buf1[i] = ellipse(buf1[i], eel1, nsects, xnorm);
        buf2[i] = ellipse(buf2[i], eel2, nsects, xnorm);
    }
}

void feed1me(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    int feed1_count;
    double *params = x->params;
    long vs = x->vs;
    double mindelay, maxdelay, speed1, speed2;
    
    double srate = x->sr;
    double *func1;
    double *func2;
    double *func3;
    double *func4;
    double *delayLine1a;
    double *delayLine2a;
    double *delayLine1b;
    double *delayLine2b;
    int *dv1a, *dv2a;        /* cmix bookkeeping */
    int *dv1b, *dv2b;        /* cmix bookkeeping */
    double delsamp1a=0, delsamp2a=0 ;
    double delsamp1b=0, delsamp2b=0 ;
    double delay1, delay2, feedback1, feedback2;
    double funcSi, funcPhs;
    double putsamp;
    double duration;
    long ddl_func_len = (long) (x->sr * MAX_MINI_DELAY);
    
    ++(*pcount);
    feed1_count = params[(*pcount)++];
    mindelay = params[(*pcount)++];
    maxdelay = params[(*pcount)++];
    speed1 = params[(*pcount)++];
    speed2 = params[(*pcount)++];
    duration = params[(*pcount)++];

    func1 = x->feed1_units[feed1_count].func1;
    func2 = x->feed1_units[feed1_count].func2;
    func3 = x->feed1_units[feed1_count].func3;
    func4 = x->feed1_units[feed1_count].func4;
    delayLine1a = x->feed1_units[feed1_count].delayLine1a;
    delayLine2a = x->feed1_units[feed1_count].delayLine2a;
    delayLine1b = x->feed1_units[feed1_count].delayLine1b;
    delayLine2b = x->feed1_units[feed1_count].delayLine2b;
    dv1a = x->feed1_units[feed1_count].dv1a;
    dv2a = x->feed1_units[feed1_count].dv2a;
    dv1b = x->feed1_units[feed1_count].dv1b;
    dv2b = x->feed1_units[feed1_count].dv2b;
    
    mindelay = x->feed1_units[feed1_count].mindelay;
    maxdelay = x->feed1_units[feed1_count].maxdelay;
    speed1 = x->feed1_units[feed1_count].speed1;
    speed2 = x->feed1_units[feed1_count].speed2;
    duration = x->feed1_units[feed1_count].duration;
    
    if( maxdelay > MAX_MINI_DELAY ){
        error("feed1: too high max delay, adjusted");
        maxdelay = MAX_MINI_DELAY;
    }
    
    funcSi = ((double) FEEDFUNCLEN / srate) / duration;
    funcPhs = x->feed1_units[feed1_count].funcPhs;
    
    for(i = 0; i < vs; i++){
        delay1 = func1[ (int) funcPhs ];
        delay2 = func2[ (int) funcPhs ];
        if(delay1 < 0.0){
            delay1 = 0.0;
        } else if (delay1 > ddl_func_len){
            delay1 = ddl_func_len - 1;
        }
        if(delay2 < 0.0){
            delay2 = 0.0;
        } else if (delay2 > ddl_func_len){
            delay2 = ddl_func_len - 1;
        }
        feedback1 = func3[ (int) funcPhs ];
        feedback2 = func4[ (int) funcPhs ];
        funcPhs += funcSi;
        if( funcPhs >= (double) FEEDFUNCLEN ){
            funcPhs = 0.0;
        }
        buf1[i] = buf1[i] + delsamp1a * feedback1;
        delput2( putsamp, delayLine1a, dv1a);
        delsamp1a = dliget2(delayLine1a, delay1, dv1a, srate);
        putsamp = delsamp1a + (delsamp2a * feedback2);
        delput2( putsamp, delayLine2a, dv2a);
        delsamp2a = dliget2(delayLine2a, delay2, dv2a, srate);
        buf1[i] += delsamp2a;
        putsamp = buf2[i] + delsamp1b * feedback1;
        buf2[i] = putsamp; // zero instead ??
        delput2( putsamp, delayLine1b, dv1b);
        delsamp1b = dliget2(delayLine1b, delay1, dv1b,srate);
        putsamp = delsamp1b + (delsamp2b * feedback2);
        delput2( putsamp, delayLine2b, dv2b);
        delsamp2b = dliget2(delayLine2b, delay2, dv2b, srate);
        buf2[i] += delsamp2b;
    }
    x->feed1_units[feed1_count].funcPhs = funcPhs;
}

void flam1(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    int flam1_count;
    int *dv1, *dv2;
    double *delayline1, *delayline2, delaytime;
    long counter, sample_length;
    double *params = x->params;
    double srate = x->sr;
    long vs = x->vs;
    double feedback;
    double delsamp1 = 0, delsamp2 = 0;
    
    ++(*pcount);
    flam1_count = params[(*pcount)++];
    sample_length = params[(*pcount)++];
    delaytime = params[(*pcount)++];

    delaytime = x->flam1_units[flam1_count].dt;
    sample_length = x->flam1_units[flam1_count].sample_length;
    counter = x->flam1_units[flam1_count].counter;
    dv1 = x->flam1_units[flam1_count].dv1;
    dv2 = x->flam1_units[flam1_count].dv2;
    delayline1 = x->flam1_units[flam1_count].delayline1;
    delayline2 = x->flam1_units[flam1_count].delayline2;

    if(delaytime >= FLAM1_MAX_DELAY){
        delaytime = FLAM1_MAX_DELAY * 0.99;
    } else if(delaytime < 0.0){
        delaytime = 0.0;
    }
    for( i = 0; i < vs; i++){

        if( counter < sample_length ){
            feedback = 1.0 - ((double)counter / (double)sample_length);
            delput2(buf1[i] + (delsamp1 * feedback), delayline1, dv1);
            delsamp1 = dliget2(delayline1, delaytime, dv1, srate);
            buf1[i] += delsamp1;
            delput2(buf2[i] + delsamp2*feedback, delayline2, dv2);
            delsamp2 = dliget2(delayline2, delaytime, dv2, srate);
            buf2[i] += delsamp2;
        }
        counter++;
        if( counter >= sample_length){
            counter = 0;
        }
    }
    x->flam1_units[flam1_count].counter = counter;
}

void comb4(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    
    int i, j;
    double insamp;
    double *params = x->params;
    long vs = x->vs;
    int comb4_count;
    double basefreq;
    int dex; // throwaway
    
    ++(*pcount);
    comb4_count = params[(*pcount)++];
    basefreq = params[(*pcount)++];
    for(i = 0; i < 4; i++){
        dex = params[(*pcount)++];
    }
/*
 
 params[pcount++] = COMB4;
 params[pcount++] = comb4_count;
 rvt = 0.99;
 params[pcount++] = basefreq = boundrand(100.0,400.0);
 for( i = 0; i < 4; i++ ){
     dex = (int)floor( boundrand(0.0, 4.0));
     lpt = 1. / basefreq;
     if( dex > 4) { dex = 4; };
     params[pcount++] = dex;
     basefreq *= ratios[dex];
     mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs1[i], x->sr);
     mycombset(lpt, rvt, 0, x->comb4_units[comb4_count].combs2[i], x->sr);
 }
 */
    for(i = 0; i < vs; i++){
        insamp = buf1[i];
        buf1[i] = 0.0;
        for(j = 0; j < 4; j++){
            buf1[i] += mycomb(insamp, x->comb4_units[comb4_count].combs1[j]) * 0.25;
        }
        insamp = buf2[i];
        buf2[i] = 0.0;
        for(j = 0; j < 4; j++){
            buf2[i] += mycomb(insamp, x->comb4_units[comb4_count].combs2[j]) * 0.25;
        }
    }
}

void compdist(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    
    int i;
    // double *params = x->params;
    long vs = x->vs;
    double *distortion_function = x->distortion_function;
    int distortion_length = x->distortion_length;
    // double insamp;
    
    ++(*pcount);
    
    for(i = 0; i < vs; i++){
        buf1[i] = dlookup(buf1[i], distortion_function, distortion_length);
        buf2[i] = dlookup(buf2[i], distortion_function, distortion_length);
    }
    
}

void ringfeed(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i; //, k;
    /* main variables */
    double *params = x->params;
    long vs = x->vs;
    /*function specific*/
    int ringfeed_count;
    double si1, si2, osc1phs, osc2phs;
    double rescale = 1.5;
    double cf, bw, lpt, rvt;
    
    ++(*pcount);
    ringfeed_count = params[(*pcount)++];
    cf = params[(*pcount)++];
    bw = params[(*pcount)++];
    si1 = params[(*pcount)++];
    si2 = params[(*pcount)++];
    lpt = params[(*pcount)++];
    rvt = params[(*pcount)++];
    /*
     params[pcount++] = RINGFEED;
     params[pcount++] = ringfeed_count;
     params[pcount++] = cf = boundrand(90.0,1500.0);
     params[pcount++] = bw = boundrand(0.02,0.4) * cf;

     params[pcount++] = x->resonfeed_units[ringfeed_count].osc1si = boundrand(90.0,1500.0) * ((double)x->sinelen / x->sr);
     params[pcount++] = x->resonfeed_units[ringfeed_count].osc2si = boundrand(90.0,1500.0) * ((double)x->sinelen / x->sr);
     params[pcount++] = lpt = 1.0 / boundrand(90.0,1500.0);
     params[pcount++] = rvt = boundrand(.01,.8);
     */
    
    si1 = x->resonfeed_units[ringfeed_count].osc1si;
    si2 = x->resonfeed_units[ringfeed_count].osc2si;
    osc1phs = x->resonfeed_units[ringfeed_count].osc1phs;
    osc2phs = x->resonfeed_units[ringfeed_count].osc2phs;

    for(i = 0; i < vs; i++){
        buf1[i] *= oscil(1.0, si1, x->sinewave, x->sinelen, &osc1phs);
        buf1[i] += mycomb(buf1[i], x->resonfeed_units[ringfeed_count].comb1arr);
        buf1[i] = reson(buf1[i], x->resonfeed_units[ringfeed_count].res1q) * rescale;
        buf2[i] *= oscil(1.0, si2,x->sinewave, x->sinelen, &osc2phs);
        buf2[i] += mycomb(buf2[i], x->resonfeed_units[ringfeed_count].comb2arr);
        buf2[i] = reson(buf2[i], x->resonfeed_units[ringfeed_count].res2q) * rescale;
    }
    x->resonfeed_units[ringfeed_count].osc1phs = osc1phs;
    x->resonfeed_units[ringfeed_count].osc2phs = osc2phs;

}

void resonadsr(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    float bwfac;
    double *q1, *q2;
    float cf, bw;
    float si;
    double phase;
    double *params = x->params;
    float srate = x->sr;
    long vs = x->vs;
    int resonadsr_count;
    double a,d,r,v1,v2,v3,v4,notedur;
    
    /*function specific*/

    int funclen;
    double *adsrfunc;
    
    ++(*pcount);
    resonadsr_count = params[(*pcount)++];
    a = params[(*pcount)++];
    d = params[(*pcount)++];
    r = params[(*pcount)++];
    v1 = params[(*pcount)++];
    v2 = params[(*pcount)++];
    v3 = params[(*pcount)++];
    v4 = params[(*pcount)++];
    bwfac = params[(*pcount)++];
    notedur = params[(*pcount)++];
 
    /*
     params[pcount++] = RESONADSR;
     params[pcount++] = resonadsr_count;
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->a = boundrand(.01,.1);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->d = boundrand(.01,.05);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->r = boundrand(.05,.5);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v1 = boundrand(150.0,4000.0);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v2 = boundrand(150.0,4000.0);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v3 = boundrand(150.0,4000.0);
     params[pcount++] = x->resonadsr_units[resonadsr_count].adsr->v4 = boundrand(150.0,4000.0);
     params[pcount++] = x->resonadsr_units[resonadsr_count].bwfac = boundrand(.03,.7);
     params[pcount++] =  notedur = boundrand(0.7, 1.2);
     
     x->resonadsr_units[resonadsr_count].phs = 0.0;
     */
    q1 = x->resonadsr_units[resonadsr_count].q1;
    q2 = x->resonadsr_units[resonadsr_count].q2;
    phase = x->resonadsr_units[resonadsr_count].phs;
    adsrfunc = x->resonadsr_units[resonadsr_count].adsr->func;
    funclen = x->resonadsr_units[resonadsr_count].adsr->len;
    si = x->resonadsr_units[resonadsr_count].si;
    bwfac = x->resonadsr_units[resonadsr_count].bwfac;

    for(i = 0; i < vs; i++){
        cf = adsrfunc[ (int) phase ];
        bw = bwfac * cf;
        phase += si;
        // this could be DANGEROUS:
        if( phase > funclen - 1){
            phase = funclen - 1;
        }
        rsnset2( cf, bw, 2.0, 1.0, q1, srate );
        rsnset2( cf, bw, 2.0, 1.0, q2, srate );
        buf1[i] = reson(buf1[i], q1);
        buf2[i] = reson(buf2[i], q2);
    }
}

void stv(t_chameleon *x, long *pcount, t_double *buf1, t_double *buf2)
{
    int i;
    /* main variables */
    double *params = x->params;
    float srate = x->sr;
    long vs = x->vs;
    int stv_count;

    /*function specific*/
    double *sinewave = x->sinewave;
    int sinelen = x->sinelen ;
    double *delayline1;
    double *delayline2;
    double fac1, fac2;
    int *dv1, *dv2; /* cmix bookkeeping */
    double delay_time;
    double si1, si2, speed1, speed2, maxdelay;
    double osc1phs, osc2phs;
    
    
    ++(*pcount);
    stv_count = params[(*pcount)++];
    speed1 = params[(*pcount)++];
    speed2 = params[(*pcount)++];
    maxdelay = params[(*pcount)++];
    dv1 = x->stv_units[stv_count].dv1;
    dv2 = x->stv_units[stv_count].dv2;
    delayline1 = x->stv_units[stv_count].delayline1;
    delayline2 = x->stv_units[stv_count].delayline2;
    si1 = x->stv_units[stv_count].osc1si;
    si2 = x->stv_units[stv_count].osc2si;
    osc1phs = x->stv_units[stv_count].osc1phs;
    osc2phs = x->stv_units[stv_count].osc2phs;
    fac1 = x->stv_units[stv_count].fac1;
    fac2 = x->stv_units[stv_count].fac2;

    for(i = 0; i < vs; i++){
        delay_time = fac1 + oscil(fac2, si1, sinewave, sinelen, &osc1phs);
        if(delay_time < 0.0){
            delay_time = 0.0;
        } else if( delay_time > STV_MAX_DELAY){
            delay_time = STV_MAX_DELAY;
        }
        delput2( buf1[i], delayline1, dv1);
        buf1[i] = dliget2(delayline1, delay_time, dv1,srate);

        delay_time = fac1 + oscil(fac2, si2, sinewave, sinelen, &osc2phs);
        if(delay_time < 0.0){
            delay_time = 0.0;
        } else if( delay_time > STV_MAX_DELAY){
            delay_time = STV_MAX_DELAY;
        }
        delput2( buf2[i], delayline2, dv2);
        buf2[i] = dliget2(delayline2, delay_time, dv2,srate);
    }

    x->stv_units[stv_count].osc1phs = osc1phs;
    x->stv_units[stv_count].osc2phs = osc2phs;
    x->stv_units[stv_count].dv1 = dv1;
    x->stv_units[stv_count].dv2 = dv2;
}

