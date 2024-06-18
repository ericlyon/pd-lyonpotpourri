#include "MSPd.h"
#include "fftease.h"

// adapted for Pd

#define CUSHION_FRAMES (128) // pad for grabbing
#define MAX_CHANNELS (2)

#define MAX_RMS_BUFFER (0.250)
#define MIN_RMS_BUFFER (.001)
#define MAX_EVENTS (1024)
#define MAX_RMS_FRAMES (32768)

#define OBJECT_NAME "buffet~"


static t_class *buffet_class;


typedef struct {
    t_word *b_samples;
    int b_valid;
    int b_nchans;
    int b_frames;
} t_guffer; // stuff we care about from garrays and buffers

typedef struct _buffet
{
    t_object x_obj;
    t_float x_f;
    t_symbol *wavename; // name of waveform buffer
    t_guffer *wavebuf; // holds waveform samples
    t_guffer *destbuf; // for copying to another buffer
    
    t_float sr; // sampling rate
    short hosed; // buffers are bad
    t_float minframes; // minimum replacement block in sample frames
    t_float maxframes; // maximum replacement block in sample frames
    long storage_maxframes; // maxframe limit that current memory can handle
    t_float *storage; //temporary memory to store replacement block (set to maxframes * channels)
    long storage_bytes; // amount of currently allocated memory
    t_float fade; // fadein/fadeout time in sample frames
    t_float sync; // input from groove sync signal
    long swapframes; // number of frames in swap block
    long r1startframe; //start frame for block 1
    long r2startframe; // start frame for block 2
    t_float dc_coef; // filter coefficient
    t_float dc_gain; // normalization factor
    short initialized; // first time or not
    t_float *rmsbuf;// for onset analysis
    t_float rmschunk; // store lowest rms value in buffer
    void *list; // for start/end list
    void *bang; // completion bang
    void *floater; // outputs noise floor
    t_atom *listdata;// to report est. start/stop times of events in buffer
    t_float *analbuf; // contain overall envelope
    t_float *onset; // contain attack times for percussive evaluations
    short autoredraw; // to kill redraw if it impacts performance
} t_buffet;

static void buffet_setbuf(t_buffet *x, t_symbol *wavename);
static void *buffet_new(t_symbol *msg, int argc, t_atom *argv);
static t_int *buffet_perform(t_int *w);
static void buffet_dsp(t_buffet *x, t_signal **sp);
static t_float buffet_boundrand(t_float min, t_float max);
static void buffet_assist (t_buffet *x, void *b, long msg, long arg, char *dst);
static void buffet_dsp_free(t_buffet *x);
static void buffet_swap(t_buffet *x);
static void buffet_specswap(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_retroblock(t_buffet *x);
static void buffet_nakedswap(t_buffet *x);
static void buffet_overlap(t_buffet *x, t_floatarg f);
static void buffet_minswap(t_buffet *x, t_floatarg f);
static void buffet_maxswap(t_buffet *x, t_floatarg f);
static void buffet_nosync_setswap(t_buffet *x);
static void buffet_info(t_buffet *x);
static void buffet_killdc(t_buffet *x);
static void buffet_rmschunk(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_fadein(t_buffet *x, t_floatarg f);
static void buffet_fadeout(t_buffet *x, t_floatarg f);
static void buffet_internal_fadein(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_internal_fadeout(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_dc_gain(t_buffet *x, t_floatarg f);
static void buffet_dc_coef(t_buffet *x, t_floatarg f);
static void buffet_normalize(t_buffet *x, t_floatarg f);
static void buffet_rotatetozero(t_buffet *x, t_floatarg f);
static void buffet_autoredraw(t_buffet *x, t_floatarg f);
static void buffet_erase(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_events(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_pevents(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_detect_onsets(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_detect_subband_onsets(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static void buffet_copy_to_buffer(t_buffet *x, t_symbol *msg, int argc, t_atom *argv);
static int buffet_setdestbuf(t_buffet *x, t_symbol *wavename);
static void buffet_init(t_buffet *x, short initialized);
static void buffet_reverse(t_buffet *x);
static void buffet_redraw(t_buffet *x);
static void buffet_redraw_named(t_buffet *x, t_symbol *arrayname);
static void buffet_update(t_buffet *x);


void buffet_tilde_setup(void)
{
    
    buffet_class = class_new(gensym("buffet~"),(t_newmethod)buffet_new,(t_method)buffet_dsp_free, sizeof(t_buffet),
                             0, A_GIMME,0);
    CLASS_MAINSIGNALIN(buffet_class,t_buffet, x_f );
    
    
    class_addmethod(buffet_class,(t_method)buffet_dsp,gensym("dsp"),A_CANT,0);
    class_addmethod(buffet_class,(t_method)buffet_swap,gensym("swap"),0);
    class_addmethod(buffet_class,(t_method)buffet_specswap,gensym("specswap"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_events,gensym("events"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_pevents,gensym("pevents"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_detect_onsets,gensym("detect_onsets"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_detect_subband_onsets,gensym("detect_subband_onsets"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_retroblock,gensym("retroblock"),0);
    class_addmethod(buffet_class,(t_method)buffet_minswap,gensym("minswap"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_maxswap,gensym("maxswap"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_fadein,gensym("fadein"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_fadeout,gensym("fadeout"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_overlap,gensym("overlap"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_dc_gain,gensym("dc_gain"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_dc_coef,gensym("dc_coef"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_normalize,gensym("normalize"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_rotatetozero,gensym("rotatetozero"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_autoredraw,gensym("autoredraw"),A_FLOAT,0);
    class_addmethod(buffet_class,(t_method)buffet_killdc,gensym("killdc"),0);
    class_addmethod(buffet_class,(t_method)buffet_nakedswap,gensym("nakedswap"),0);
    class_addmethod(buffet_class,(t_method)buffet_rmschunk,gensym("rmschunk"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_erase,gensym("erase"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_copy_to_buffer,gensym("copy_to_buffer"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_internal_fadein,gensym("internal_fadein"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_internal_fadeout,gensym("internal_fadeout"),A_GIMME,0);
    class_addmethod(buffet_class,(t_method)buffet_setbuf,gensym("setbuf"),A_SYMBOL,0);
    class_addmethod(buffet_class,(t_method)buffet_reverse,gensym("reverse"),0);
    class_addmethod(buffet_class,(t_method)buffet_info,gensym("info"),0);
    potpourri_announce(OBJECT_NAME);
}



void buffet_info(t_buffet *x)
{
    long totalframes;
    buffet_setbuf(x, x->wavename);
    
    if( x->hosed ) {
        pd_error(0, "buffet~ needs a valid buffer");
        return;
    }
    
    if( ! x->sr) {
        pd_error(0, "zero sample rate!");
        return;
    }
    totalframes = x->wavebuf->b_frames;
    post("minswap: %f, maxswap: %f", 1000. * x->minframes / x->sr, 1000. * x->maxframes / x->sr);
    post("buffer size: %f", 1000. * totalframes / x->sr);
}

void buffet_overlap(t_buffet *x, t_floatarg f)
{
    if( f < 0.01 ) {
        pd_error(0, "minimum fade time is 0.01 milliseconds");
        return;
    }
    x->fade = f * .001 * x->sr;
}

void buffet_events(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_word *b_samples;
    int b_nchans;
    int b_frames;
    t_atom *listdata = x->listdata;
    
    t_float bufsize;
    t_float onthresh;
    t_float offthresh;
    long bufsamps;
    long aframes; // frames to analyze
    t_float tadv;
    
    long i,j;
    t_float meansq;
    t_float rmsval;
    long bindex;
    t_float ipos;
    short activated = 0;
    t_float realtime = 0.0;
    int event_count = 0;
    t_float buffer_duration;
    
    buffet_setbuf(x, x->wavename);
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    // duration in ms.
    buffer_duration = 1000.0 * (t_float)b_frames / x->sr;
    
    bufsize = .001 * atom_getfloatarg(0,argc,argv);
    if(bufsize > MAX_RMS_BUFFER) {
        bufsize = MAX_RMS_BUFFER;
        post("%s: setting analysis buffer to maximum: %f",OBJECT_NAME, MAX_RMS_BUFFER * 1000.0);
    } else if(bufsize < MIN_RMS_BUFFER) {
        bufsize = MIN_RMS_BUFFER;
        post("%s: setting analysis buffer to maximum: %f",OBJECT_NAME, MIN_RMS_BUFFER * 1000.0);
    }
    
    onthresh = atom_getfloatarg(1,argc,argv);
    offthresh = atom_getfloatarg(2,argc,argv);
    bufsamps = x->sr * bufsize;
    
    bufsamps = bufsize * x->sr;
    tadv = (t_float)bufsamps / x->sr;
    //  post("actual window size: %f",tadv);
    aframes = (long) ( (t_float) b_frames / (t_float)bufsamps ) - 1;
    if(aframes < 2) {
        pd_error(0, "%s: this buffer is too short to analyze",OBJECT_NAME);
        return;
    }
    //  post("analyzing %d frames",aframes);
    for(i = 0; i < aframes; i++) {
        meansq = 0.0;
        ipos = b_nchans * i * bufsamps;
        /* only analyze first channel */
        for(j = 0; j < bufsamps; j+= b_nchans) {
            bindex = ipos + j;
            meansq += b_samples[bindex].w_float * b_samples[bindex].w_float;
        }
        meansq /= (t_float) bufsamps;
        rmsval = sqrt(meansq);
        realtime += tadv;
        if(rmsval > onthresh && ! activated) {
            activated = 1;
            
            //      post("event %d starts at %f",event_count+1, realtime);
            
            if(event_count >= MAX_EVENTS) {
                pd_error(0, "%s: exceeded maximum of %d events",OBJECT_NAME, MAX_EVENTS);
                break;
            }
            SETFLOAT(listdata+(event_count*2), realtime * 1000.0);
            
        }
        else if( rmsval < offthresh && activated ) {
            activated = 0;
            //      post("event %d ends at %f",event_count, realtime);
            SETFLOAT(listdata+((event_count*2) + 1), realtime * 1000.0);
            ++event_count;
        }
    }
    if(activated) {
        post("%s: missed the end of the last event; setting to end of buffer",OBJECT_NAME);
        SETFLOAT(listdata+((event_count*2) + 1), buffer_duration);
        ++event_count;
    }
    outlet_list(x->list, 0, event_count * 2, listdata);
}


void buffet_pevents(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_word *b_samples;
    int b_nchans;
    int b_frames;
    t_atom *listdata = x->listdata;
    
    t_float bufsize;
    //t_float onthresh;
    //t_float offthresh;
    t_float diffthresh;
    long bufsamps;
    long aframes; // frames to analyze
    t_float tadv;
    //long current_frame = 0;
    long i,j;
    t_float meansq;
    t_float rmsval;
    long bindex;
    t_float ipos;
    //short activated = 0;
    t_float realtime = 0.0;
    int event_count = 0;
//    t_float buffer_duration;
    t_float mindiff, maxdiff, absdiff,rmsdiff;
    t_float *analbuf = x->analbuf;
    t_float *onset = x->onset;
    //t_float endtime;
    
    buffet_setbuf(x, x->wavename);
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    // duration in ms.
//    buffer_duration = 1000.0 * (t_float)b_frames / x->sr;
    
    bufsize = .001 * atom_getfloatarg(0,argc,argv);
    if(bufsize > MAX_RMS_BUFFER) {
        bufsize = MAX_RMS_BUFFER;
        post("%s: setting analysis buffer to maximum: %f",OBJECT_NAME, MAX_RMS_BUFFER * 1000.0);
    } else if(bufsize < MIN_RMS_BUFFER) {
        bufsize = MIN_RMS_BUFFER;
        post("%s: setting analysis buffer to maximum: %f",OBJECT_NAME, MIN_RMS_BUFFER * 1000.0);
    }
    
    diffthresh = atom_getfloatarg(1,argc,argv);
    bufsamps = x->sr * bufsize;
    
    bufsamps = bufsize * x->sr;
    tadv = (t_float)bufsamps / x->sr;
    
    aframes = (long) ( (t_float) b_frames / (t_float)bufsamps );
    if(aframes < 2) {
        pd_error(0, "%s: this buffer is too short to analyze",OBJECT_NAME);
        return;
    }
    if(aframes > MAX_RMS_FRAMES) {
        post("too many frames: try a larger buffer size");
        return;
    }
    analbuf[0] = 0; // for first comparison
    for(i = 1; i < aframes; i++) {
        meansq = 0.0;
        ipos = b_nchans * i * bufsamps;
        /* only analyze first channel */
        for(j = 0; j < bufsamps; j+= b_nchans) {
            bindex = ipos + j;
            meansq += b_samples[bindex].w_float * b_samples[bindex].w_float;
        }
        meansq /= (t_float) bufsamps;
        analbuf[i] = rmsval = sqrt(meansq);
        realtime += tadv;
        
    }
    
    realtime = 0;
    mindiff = 9999.;
    maxdiff = 0.0;
    /* LPF - use elsewhere?
     for(i = 1; i < aframes; i++ ) {
     analbuf[i] = (analbuf[i] + analbuf[i-1]) * 0.5;
     }
     */
    // look for big changes in direction
    for(i = 1; i < aframes; i++ ) {
        rmsdiff = analbuf[i] - analbuf[i-1];
        absdiff = fabs(rmsdiff);
        if(absdiff > maxdiff)
            maxdiff = absdiff;
        if(absdiff < mindiff)
            mindiff = absdiff;
        
        if( rmsdiff > diffthresh ) {
            // new
            
            onset[event_count] = (realtime + bufsize) * 1000.0;
            if(onset[event_count] < 0)
                onset[event_count] = 0;
            ++event_count;
            //      post("rt %f diff %f",realtime * 1000.0,rmsdiff);
        }
        realtime += tadv;
    }
    //  post("mindiff %f maxdiff %f",mindiff,maxdiff);
    if(event_count == 0) {
        post("%s: no events found",OBJECT_NAME);
    }
    
    
    for(i = 0; i < event_count; i++) {
        SETFLOAT(listdata + i, onset[i]);
    }
    outlet_list(x->list, 0, event_count, listdata);
}

void buffet_internal_fadeout(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    long fadeframes;
    long totalframes;
    int i,j,k;
    t_float env;
    t_word *b_samples;
    long b_nchans;
    long startframe;
    long endframe;
    
    if( ! x->sr) {
        pd_error(0, "zero sample rate!");
        return;
    }
    
    
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(argc < 2) {
        post("%s: internal_fadeout requires start and end times",OBJECT_NAME);
        return;
    }
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    
    if(startframe < 0 || endframe > totalframes || endframe <= startframe) {
        pd_error(0, "%s: bad frame boundaries to internal_fadeout: %ld and %ld",OBJECT_NAME,startframe, endframe);
        return;
    }
    fadeframes = endframe - startframe;
    
    for(i = (endframe-1) * b_nchans , k = 0; k < fadeframes ; i -= b_nchans, k++ ) {
        env = (t_float) k / (t_float) fadeframes;
        for(j = 0; j < b_nchans; j++) {
            b_samples[i + j].w_float *= env;
        }
    }
    buffet_update(x);
}

void buffet_internal_fadein(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    long fadeframes;
    long totalframes;
    int i,j,k;
    t_float env;
    t_word *b_samples;
    long b_nchans;
    long startframe;
    long endframe;
    
    if( ! x->sr) {
        pd_error(0, "zero sample rate!");
        return;
    }
    
    
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(argc < 2) {
        post("%s: internal_fadeout requires start and end times",OBJECT_NAME);
        return;
    }
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    
    if(startframe < 0 || endframe > totalframes || endframe <= startframe) {
        pd_error(0, "%s: bad frame boundaries to internal_fadein: %ld and %ld",OBJECT_NAME,startframe, endframe);
        return;
    }
    fadeframes = endframe - startframe;
    
    for(i = startframe * b_nchans , k = 0; k < fadeframes ; i += b_nchans, k++ ) {
        env = (t_float) k / (t_float) fadeframes;
        for(j = 0; j < b_nchans; j++) {
            b_samples[i + j].w_float *= env;
        }
    }
    buffet_update(x);
}

void buffet_reverse(t_buffet *x)
{
    
    int i,j;
    //  t_float env;
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    t_float tmpsamp;
    long lenm1;
    
    buffet_setbuf(x, x->wavename);
    if(!x->wavebuf->b_valid) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    lenm1 = (b_frames - 1) * b_nchans;
    for( i = 0; i < (b_frames * b_nchans) / 2; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            tmpsamp = b_samples[(lenm1 - i) + j].w_float;
            b_samples[(lenm1 - i) + j].w_float = b_samples[i + j].w_float;
            b_samples[i + j].w_float = tmpsamp;
        }
    }
    
    buffet_update(x);
}

void buffet_update(t_buffet *x)
{
    outlet_bang(x->bang);
    if(x->autoredraw) {
        buffet_redraw(x);
    }
}


void buffet_redraw_named(t_buffet *x, t_symbol *arrayname)
{
    t_garray *a;
    if (!(a = (t_garray *)pd_findbyclass(arrayname, garray_class))) {
        if (*arrayname->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME, arrayname->s_name);
    }
    else  {
        garray_redraw(a);
    }
}

void buffet_redraw(t_buffet *x)
{
    t_garray *a;
    if (!(a = (t_garray *)pd_findbyclass(x->wavename, garray_class))) {
        if (*x->wavename->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME, x->wavename->s_name);
        x->wavebuf->b_valid = 0;
    }
    else  {
        garray_redraw(a);
    }
}

void buffet_copy_to_buffer(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_symbol *destname;
    
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    
    t_word *b_dest_samples;
    long b_dest_nchans;
    long b_dest_frames;
    
    long startframe;
    long endframe;
    long chunkframes;
    
    t_float fadein,fadeout;
    int fadeframes;
    // totalframes;
    t_float env;
    int i,j,k;
    
    buffet_setbuf(x, x->wavename);
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    
    destname = atom_getsymbolarg(0,argc,argv);
    
    if(! buffet_setdestbuf(x, destname)) {
        post("%s: could not find buffer %s",OBJECT_NAME,destname->s_name);
        return;
    }
    b_dest_samples = x->destbuf->b_samples;
    b_dest_nchans = x->destbuf->b_nchans;
    b_dest_frames = x->destbuf->b_frames;
    
    startframe = atom_getfloatarg(1,argc,argv) * .001 * x->sr;
    endframe = atom_getfloatarg(2,argc,argv) * .001 * x->sr;
    chunkframes = endframe - startframe;
    if(chunkframes <= 0) {
        return;
    }
    if(b_nchans != b_dest_nchans) {
        pd_error(0, "%s: channel mismatch with buffer %s",OBJECT_NAME, destname->s_name);
        return;
    }
    if(b_dest_frames < chunkframes ) {
        //    post("%s: %s is too small, truncating copy region",OBJECT_NAME, destname->s_name);
        chunkframes = b_dest_frames;
    }
    if(startframe < 0 || endframe >= b_frames) {
        pd_error(0, "%s: bad frame range for source buffer: %ld %ld",OBJECT_NAME, startframe,endframe);
        return;
    }
    /* first clean out destination */
    for(i = 0; i < b_dest_frames; i++) {
        b_dest_samples[i].w_float = 0.0;
    }
    // memset((char *)b_dest_samples, 0, b_dest_frames * b_dest_nchans * sizeof(t_float));
    
    /* now copy segment */
    for(i = 0; i < chunkframes; i++) {
        b_dest_samples[i].w_float = b_samples[i+startframe].w_float;
    }
    
    /* memcpy(b_dest_samples, b_samples + (startframe * b_nchans),
     chunkframes * b_nchans * sizeof(t_float) ); */
    
    if(argc == 5) {
        //    post("enveloping");
        fadein = atom_getfloatarg(3,argc,argv);
        fadeout = atom_getfloatarg(4,argc,argv);
        if(fadein > 0) {
            //      post("fading in");
            fadeframes = fadein * .001 * x->sr;
            
            if( fadeframes > b_dest_frames) {
                pd_error(0, "%s: fadein is too long",OBJECT_NAME);
                return;
            }
            
            
            for(i = 0 , k = 0; k < fadeframes ; i += b_dest_nchans, k++ ) {
                env = (t_float) k / (t_float) fadeframes;
                for(j = 0; j < b_dest_nchans; j++) {
                    b_dest_samples[i + j].w_float *= env;
                }
            }
        }
        
        if(fadeout > 0) {
            //      post("fading out");
            //      totalframes = chunkframes;
            fadeframes = fadeout * .001 * x->sr;
            startframe = chunkframes - fadeframes;
            endframe = chunkframes;
            
            //    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv); // stays the same
            
            if(startframe < 0) {
                pd_error(0, "%s: bad frame boundaries to internal_fadeout: %ld and %ld",
                      OBJECT_NAME,startframe, endframe);
                return;
            }
            //    fadeframes = endframe - startframe; // we already know this
            
            for(i = (chunkframes-1) * b_dest_nchans , k = 0; k < fadeframes ; i -= b_dest_nchans, k++ ) {
                env = (t_float) k / (t_float) fadeframes;
                for(j = 0; j < b_dest_nchans; j++) {
                    b_dest_samples[i + j].w_float *= env;
                }
            }
        }
    }
    buffet_redraw_named(x, destname);
    outlet_bang(x->bang);
}

int buffet_setdestbuf(t_buffet *x, t_symbol *wavename)
{
    
    t_garray *a;
    int b_frames;
    t_word *b_samples;
    if (!(a = (t_garray *)pd_findbyclass(wavename, garray_class))) {
        if (*wavename->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME,wavename->s_name);
        x->destbuf->b_valid = 0;
    }
    else if (!garray_getfloatwords(a, &b_frames, &b_samples)) {
        pd_error(x, "%s: bad array for %s", wavename->s_name,OBJECT_NAME);
        x->destbuf->b_valid = 0;
    }
    else  {
        x->destbuf->b_valid = 1;
        x->destbuf->b_frames = (long)b_frames;
        x->destbuf->b_nchans = 1;
        x->destbuf->b_samples = b_samples;
        garray_usedindsp(a);
    }
    /*We added this in linux - better check that it works */
    return x->destbuf->b_valid;
}

void buffet_rmschunk(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    
    long bufsamps;
    long i;
    t_float meansq;
    t_float rmsval;
    long bindex;
    
//    t_float buffer_duration;
    long startframe;
    long endframe;
    
    buffet_setbuf(x, x->wavename);
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    // duration in ms.
//    buffer_duration = 1000.0 * (t_float)b_frames / x->sr;
    
    
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    if(startframe < 0 || startframe >= b_frames - 1) {
        pd_error(0, "%s: naughty start frame: %ld",OBJECT_NAME,startframe);
        return;
    }
    if(endframe < 2 || endframe >= b_frames) {
        pd_error(0, "%s: naughty start frame: %ld",OBJECT_NAME,startframe);
        return;
    }
    
    
    bufsamps = (endframe - startframe);
    //  post("analyzing %d samples",bufsamps);
    if(!bufsamps) {
        //  post("%s: start and end are identical!",OBJECT_NAME);
        //  post("instantaneous sample value [ch 1 only]: %f", b_samples[startframe*b_nchans]);
        return;
    }
    meansq = 0.0;
    for(i = startframe; i < endframe; i++) {
        
        bindex = b_nchans * i;    /* only analyze first channel */
        meansq += b_samples[bindex].w_float * b_samples[bindex].w_float;
        
    }
    
    meansq /= (t_float) bufsamps;
    rmsval = sqrt(meansq);
    
    x->rmschunk = rmsval;
    outlet_float(x->floater, (double)x->rmschunk);
}


void buffet_dc_gain(t_buffet *x, t_floatarg f)
{
    x->dc_gain = f;
}

void buffet_dc_coef(t_buffet *x, t_floatarg f)
{
    x->dc_coef = f;
}

void buffet_autoredraw(t_buffet *x, t_floatarg f)
{
    x->autoredraw = (short)f;
}

void buffet_minswap(t_buffet *x, t_floatarg f)
{
    if(f < 2.0 * 1000.0 * x->fade  / x->sr) {
        pd_error(0, "minimum must be at least twice fade time which is %f", 1000. * x->fade / x->sr);
        return;
    }
    x->minframes = f * .001 * x->sr;
    //  post("min set to %f samples",x->minframes);
}

void buffet_maxswap(t_buffet *x, t_floatarg f)
{
    long oldmem;
    //  long framelimit;
    long newframes;
    
    newframes = f * .001 * x->sr;
    if(newframes <= x->minframes) {
        pd_error(0, "max blocksize must exceed minimum blocksize, which is %f", 1000. * x->minframes/ x->sr);
    }
    if(newframes > x->storage_maxframes ) {
        //  post("extending memory");
        oldmem = x->storage_bytes;
        //  post("old memory %d", oldmem);
        x->storage_maxframes = newframes;
        x->storage_bytes = (x->storage_maxframes + 1) * 2 * sizeof(t_float);
        //  post("new memory %d", x->storage_bytes);
        
        x->storage = (t_float *) resizebytes((char *)x->storage, oldmem, x->storage_bytes);
    }
    x->maxframes = newframes;
}

void buffet_erase(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    long startframe, endframe;
    //int slab;
    //int offset;
    int i;
    
    if(argc < 2) {
        post("%s: erase requires start and end times",OBJECT_NAME);
        return;
    }
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    if(startframe < 0) {
        startframe = 0;
    }
    if(endframe > b_frames - 1) {
        endframe = b_frames - 1;
    }
    if(startframe >= b_frames - 1) {
        pd_error(0, "%s: naughty start frame: %ld",OBJECT_NAME,startframe);
        return;
    }
    if(endframe < 2 || endframe <= startframe) {
        pd_error(0, "%s: naughty end frame: %ld",OBJECT_NAME,endframe);
        return;
    }
    
    for(i = startframe * b_nchans; i < endframe * b_nchans; i++) {
        b_samples[i].w_float = 0.0;
    }
    
    buffet_update(x);
}


void buffet_rotatetozero(t_buffet *x, t_floatarg f)
{
    int i;
    t_float target = (t_float) f;
    long shiftframes = (long) (target * 0.001 * x->sr);
    
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    t_float *tmpmem;
    
    
    
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    
    
    if(shiftframes <= 0 || shiftframes >= b_frames) {
        pd_error(0, "%s: shift target %f is out of range",OBJECT_NAME,target);
        return;
    }
    
    tmpmem = (t_float *) getbytes(shiftframes * b_nchans * sizeof(t_float));
    
    /* copy shift block to tmp */
    for(i = 0; i < shiftframes; i++) {
        tmpmem[i] = b_samples[i].w_float;
    }
    
    /* now shift the rest to the top */
    for(i = 0; i < b_frames - shiftframes; i++) {
        b_samples[i].w_float =b_samples[i+shiftframes].w_float;
    }

    /* finally copy tmp to the tail */
    
    for(i = 0; i < shiftframes; i++) {
        b_samples[(b_frames - shiftframes)+i].w_float = tmpmem[i];
    }

    freebytes(tmpmem, shiftframes * b_nchans * sizeof(t_float));
    buffet_update(x);
}

void buffet_normalize(t_buffet *x, t_floatarg f)
{
    
    t_float target = (t_float) f;
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    long i;
    t_float maxamp = 0.0;
    t_float amptest;
    t_float rescale;
    
    if(target <= 0.0) {
        pd_error(0, "%s: normalize target %f is too low",OBJECT_NAME,target);
        return;
    }
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    for(i = 0; i < b_frames * b_nchans; i++) {
        amptest = fabs(b_samples[i].w_float);
        if(maxamp < amptest)
            maxamp = amptest;
    }
    
    if(maxamp < .000000001) {
        post("%s: amplitude zero or too low to normalize in \"%s\"",OBJECT_NAME,x->wavename->s_name);
        return;
    }
    rescale = target / maxamp;
    if(rescale > .99 && rescale < 1.01) {
        post("%s: \"%s\" already normalized to %f",OBJECT_NAME,x->wavename->s_name,target);
    }
    else {
        for(i = 0; i < b_frames * b_nchans; i++) {
            b_samples[i].w_float *= rescale;
        }
    }
    buffet_update(x);
}


void buffet_fadein(t_buffet *x, t_floatarg f)
{
    long fadeframes;
    long totalframes;
    int i,j,k;
    t_float env;
    t_word *b_samples;
    long b_nchans;
    
    
    if( ! x->sr) {
        pd_error(0, "zero sample rate!");
        return;
    }
    fadeframes = f * .001 * x->sr;
    //post("fading in %d frames",fadeframes);
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    
    
    if( fadeframes > totalframes) {
        pd_error(0, "fadein is too long");
        return;
    }
    for(i = 0, k = 0;i < fadeframes * b_nchans; i += b_nchans, k++ ) {
        env = (t_float) k / (t_float) fadeframes;
        for(j = 0; j < b_nchans; j++) {
            b_samples[i + j].w_float *= env;
        }
    }
    buffet_update(x);
}

void buffet_fadeout(t_buffet *x, t_floatarg f)
{
    long fadeframes;
    long totalframes;
    int i,j,k;
    t_float env;
    t_word *b_samples;
    long b_nchans;
    
    
    if( ! x->sr) {
        pd_error(0, "zero sample rate!");
        return;
    }
    fadeframes = f * .001 * x->sr;
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    
    
    if( fadeframes > totalframes) {
        pd_error(0, "%s: fadein is too long",OBJECT_NAME);
        return;
    }
    for(i = (totalframes-1) * b_nchans , k = 0; k < fadeframes ; i -= b_nchans, k++ ) {
        env = (t_float) k / (t_float) fadeframes;
        for(j = 0; j < b_nchans; j++) {
            b_samples[i + j].w_float *= env;
        }
    }
    buffet_update(x);
}



void buffet_killdc(t_buffet *x)
{
    //  long fadeframes;
    long totalframes;
    int i,j;
    t_word *b_samples;
    long b_nchans;
    t_float dc_coef = x->dc_coef;
    t_float a0[MAX_CHANNELS];
    t_float a1[MAX_CHANNELS];
    t_float b0[MAX_CHANNELS];
    t_float b1[MAX_CHANNELS];
    
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(b_nchans > MAX_CHANNELS) {
        pd_error(0, "buffer has too many channels");
        return;
    }
    for(j = 0; j < b_nchans; j++) {
        a0[j] = a1[j] = b0[j] = b1[j] = 0.0;
    }
    
    for(i = 0; i < totalframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            a0[j] = b_samples[i + j].w_float;
            b0[j] = a0[j] - a1[j] + dc_coef * b1[j];
            b_samples[i + j].w_float = b0[j];
            b1[j] = b0[j];
            a1[j] = a0[j];
        }
    }
    buffet_update(x);
}


void *buffet_new(t_symbol *msg, int argc, t_atom *argv)
{
    srand(clock());
    t_buffet *x = (t_buffet *)pd_new(buffet_class);
    x->bang = outlet_new(&x->x_obj, gensym("bang"));
    x->list = outlet_new(&x->x_obj, gensym("list"));
    x->floater = outlet_new(&x->x_obj, gensym("float"));
    x->sr = sys_getsr();
    if(! x->sr )
        x->sr = 44100;
    
    x->minframes = 0;
    x->maxframes = 0;
    if(argc < 1) {
        pd_error(0, "%s: you must provide a valid buffer name",OBJECT_NAME);
        x->wavename = &s_;
    }
    atom_arg_getsym(&x->wavename,0,argc,argv);
    atom_arg_getfloat(&x->minframes,1,argc,argv);
    atom_arg_getfloat(&x->maxframes,2,argc,argv);
    if(!x->minframes)
        x->minframes = 100;
    if(!x->maxframes)
        x->maxframes = x->minframes + 10;
    buffet_init(x,0);
    
    return x;
}

void buffet_init(t_buffet *x, short initialized)
{
    if(x->minframes <= 0)
        x->minframes = 250;
    if(x->maxframes <= 0)
        x->maxframes = 1000;
    
    if(!initialized) {
        x->minframes *= .001 * x->sr;
        x->storage_maxframes = x->maxframes *= .001 * x->sr;
        x->fade = .001 * 20 * x->sr; // 20 ms fadetime to start
        x->storage_bytes = (x->maxframes + 1) * 2 * sizeof(t_float); // stereo storage frames
        x->storage = (t_float *) getbytes(x->storage_bytes);
        x->dc_coef = .995; // for dc blocker
        x->dc_gain = 4.0;
        x->autoredraw = 1;
        x->rmsbuf = (t_float *) getbytes(MAX_RMS_BUFFER * x->sr * sizeof(t_float));
        memset((char *)x->rmsbuf, 0, MAX_RMS_BUFFER * x->sr * sizeof(t_float));
        x->listdata = (t_atom *)getbytes(MAX_EVENTS * sizeof(t_atom));// lots of events
        x->analbuf = (t_float *) getbytes(MAX_RMS_FRAMES * sizeof(t_float));
        memset((char *)x->analbuf, 0, MAX_RMS_FRAMES * sizeof(t_float));
        x->onset = (t_float *) getbytes(MAX_EVENTS * sizeof(t_float));
        x->wavebuf = (t_guffer *) getbytes(1 * sizeof(t_guffer));
        x->destbuf = (t_guffer *) getbytes(1 * sizeof(t_guffer));
    } else {
        x->minframes *= .001 * x->sr;
        x->storage_maxframes = x->maxframes *= .001 * x->sr;
        x->fade = .001 * 20 * x->sr; // 20 ms fadetime to start
        x->storage_bytes = (x->maxframes + 1) * 2 * sizeof(t_float); // stereo storage frames
        x->storage = (t_float *) resizebytes((char *)x->storage_bytes,0,x->storage_bytes);
        x->rmsbuf = (t_float *)resizebytes((char *)x->rmsbuf,0,MAX_RMS_BUFFER * x->sr * sizeof(t_float));
        memset((char *)x->rmsbuf, 0, MAX_RMS_BUFFER * x->sr * sizeof(t_float));
    }
}


void buffet_nosync_setswap(t_buffet *x)
{
    long totalframes = x->wavebuf->b_frames;
    
    t_float minframes = x->minframes;
    t_float maxframes = x->maxframes;
    long swapframes = x->swapframes;
    long r1startframe = x->r1startframe;
    long r2startframe = x->r2startframe;
    long r1endframe;
    long r2endframe;
    long region1;
    long region2;
    
    
    swapframes = buffet_boundrand(minframes, maxframes);
    r1startframe = buffet_boundrand(0.0, (t_float)(totalframes-swapframes));
    r1endframe = r1startframe + swapframes;
    region1 = r1startframe;
    region2 = totalframes - r1endframe;
    if(swapframes > region1) {
        r2startframe = buffet_boundrand((t_float)r1endframe,(t_float)(totalframes-swapframes));
    } else if(swapframes > region2) {
        r2startframe = buffet_boundrand(0.0,(t_float)(r1startframe-swapframes));
    } else { // either region ok
        if(buffet_boundrand(0.0,1.0) > 0.5) {
            r2startframe = buffet_boundrand(0.0,(t_float)(r1startframe-swapframes));
        } else {
            r2startframe = buffet_boundrand((t_float)r1endframe,(t_float)(totalframes-swapframes));
        }
    }
    r2endframe = r2startframe + swapframes;
    
    if(r2startframe < 0 || r1startframe < 0) {
        pd_error(0, "start frame less than zero!");
        return;
    }
    if(r2endframe >= totalframes || r1endframe >= totalframes) {
        pd_error(0, "end frame reads beyond buffer!");
        return;
    }
    x->swapframes = swapframes;
    x->r1startframe = r1startframe;
    x->r2startframe = r2startframe;
    outlet_bang(x->bang);
}


void buffet_swap(t_buffet *x)
{
    
    t_float maxframes = x->maxframes;
    t_float fade = x->fade;
    long totalframes;
    long swapframes;
    long start_sample1;
    long start_sample2;
    t_float *storage = x->storage;
    t_word *b_samples;
    t_float mix_sample;
    t_float fadein_gain;
    t_float fadeout_gain;
    t_float fadephase;
    long b_nchans;
    int i,j,k;
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(totalframes < maxframes * 2 + 1) {
        pd_error(0, "buffer must contain at least twice as many samples as the maximum swap size");
        return;
    }
    if(b_nchans > 2) {
        pd_error(0, "buffet~ only accepts mono or stereo buffers");
        return;
    }
    
    buffet_nosync_setswap(x);
    
    
    start_sample1 = x->r1startframe * b_nchans;
    start_sample2 = x->r2startframe * b_nchans;
    swapframes = x->swapframes;
    
    // store block1 samples
    for(i = 0; i < swapframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            storage[i+j] = b_samples[start_sample1 + i + j].w_float;
        }
    }
    // swap block2 into block1 location, fadein
    
    for(i = 0, k = 0; i < fade * b_nchans; i += b_nchans, k++) {
        fadephase = ((t_float)k / (t_float) fade) * PIOVERTWO;
        fadein_gain = sin(fadephase);
        fadeout_gain = cos(fadephase);
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadein_gain * b_samples[start_sample2 + i + j].w_float + fadeout_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    // middle part, pure swap
    
    for(i = fade * b_nchans; i < (swapframes-fade) * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample1 + i + j].w_float = b_samples[start_sample2 + i + j].w_float;
        }
    }
    // fade out
    
    for(i = (swapframes-fade) * b_nchans, k = 0; i < swapframes * b_nchans; i += b_nchans, k++) {
        fadephase = ((t_float)k / (t_float) fade) * PIOVERTWO;
        fadein_gain = sin(fadephase);
        fadeout_gain = cos(fadephase);
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadeout_gain * b_samples[start_sample2 + i + j].w_float + fadein_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    // now mix stored block1 into block2 location
    // swap block2 into block1 location, fadein
    
    for(i = 0, k = 0; i < fade * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadein_gain * storage[i + j] + fadeout_gain * b_samples[start_sample2 + i + j].w_float;
            b_samples[start_sample2 + i + j].w_float = mix_sample;
        }
    }
    // middle part, pure swap
    
    for(i = fade * b_nchans; i < (swapframes-fade) * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample2 + i + j].w_float = storage[i + j];
        }
    }
    // fade out
    
    for(i = (swapframes-fade) * b_nchans, k = 0; i < swapframes * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadeout_gain * storage[i + j] + fadein_gain * b_samples[start_sample2 + i + j].w_float;
            b_samples[start_sample2 + i + j].w_float = mix_sample;
        }
    }
    buffet_update(x);
}


void buffet_specswap(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    //  t_float minframes = x->minframes;
    t_float maxframes = x->maxframes;
    t_float fade = x->fade;
    long totalframes;
    long swapframes;
    long r1startframe;
    long r1endframe;
    long r2startframe;
    long r2endframe;
    long region1;
    long region2;
    long start_sample1;
    long start_sample2;
    t_float *storage = x->storage;
    t_word *b_samples;
    t_float mix_sample;
    t_float fadein_gain;
    t_float fadeout_gain;
    t_float fadephase;
    long b_nchans;
    int i,j,k;
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(totalframes < maxframes * 2 + 1) {
        pd_error(0, "buffer must contain at least twice as many samples as the maximum swap size");
        return;
    }
    if(b_nchans > 2) {
        pd_error(0, "buffet~ only accepts mono or stereo buffers");
        return;
    }
    
    r1startframe = x->sr * .001 * atom_getfloatarg(0,argc,argv);
    r2startframe = x->sr * .001 * atom_getfloatarg(1,argc,argv);
    swapframes = x->sr * .001 * atom_getfloatarg(2,argc,argv);
    if( r1startframe < 0 || r1startframe >= totalframes) {
        pd_error(0, "bad first skip time");
        return;
    }
    if( r2startframe < 0 || r2startframe >= totalframes) {
        pd_error(0, "bad second skip time");
        return;
    }
    //  post("start1 %d start2 %d swaps %d",r1startframe,r2startframe,swapframes );
    
    r1endframe = r1startframe + swapframes;
    r2endframe = r2startframe + swapframes;
    region1 = r1startframe;
    region2 = totalframes - r1endframe;
    /*
     post("min: %.0f, max: %.0f, total: %d, swap: %d",minframes, maxframes,totalframes, swapframes);
     post("r1st %d, r1end %d, region1 %d, region2 %d",r1startframe, r1endframe, region1, region2);
     */
    if(swapframes > x->storage_maxframes) {
        pd_error(0, "swapsize %ld is larger than %ld; reset maximum swap.", swapframes,x->storage_maxframes );
        return;
    }
    if(r1endframe >= totalframes) {
        pd_error(0, "block 1 reads beyond buffer!");
        return;
    }
    if(r2endframe >= totalframes) {
        pd_error(0, "block 2 reads beyond buffer!");
        return;
    }
    
    // store block1 samples
    start_sample1 = r1startframe * b_nchans;
    start_sample2 = r2startframe * b_nchans;
    for(i = 0; i < swapframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            storage[i+j] = b_samples[start_sample1 + i + j].w_float;
        }
    }
    // swap block2 into block1 location, fadein
    
    for(i = 0, k = 0; i < fade * b_nchans; i += b_nchans, k++) {
        fadephase = ((t_float)k / (t_float) fade) * PIOVERTWO;
        fadein_gain = sin(fadephase);
        fadeout_gain = cos(fadephase);
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadein_gain * b_samples[start_sample2 + i + j].w_float + fadeout_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    // middle part, pure swap
    
    for(i = fade * b_nchans; i < (swapframes-fade) * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample1 + i + j].w_float = b_samples[start_sample2 + i + j].w_float;
        }
    }
    // fade out
    
    for(i = (swapframes-fade) * b_nchans, k = 0; i < swapframes * b_nchans; i += b_nchans, k++) {
        fadephase = ((t_float)k / (t_float) fade) * PIOVERTWO;
        fadein_gain = sin(fadephase);
        fadeout_gain = cos(fadephase);
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadeout_gain * b_samples[start_sample2 + i + j].w_float + fadein_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    // now mix stored block1 into block2 location
    // swap block2 into block1 location, fadein
    
    for(i = 0, k = 0; i < fade * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadein_gain * storage[i + j] + fadeout_gain * b_samples[start_sample2 + i + j].w_float;
            b_samples[start_sample2 + i + j].w_float = mix_sample;
        }
    }
    // middle part, pure swap
    
    for(i = fade * b_nchans; i < (swapframes-fade) * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample2 + i + j].w_float = storage[i + j];
        }
    }
    // fade out
    
    for(i = (swapframes-fade) * b_nchans, k = 0; i < swapframes * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadeout_gain * storage[i + j] + fadein_gain * b_samples[start_sample2 + i + j].w_float;
            b_samples[start_sample2 + i + j].w_float = mix_sample;
        }
    }
    buffet_update(x);
}


void buffet_retroblock(t_buffet *x)
{
    t_float minframes = x->minframes;
    t_float maxframes = x->maxframes;
    t_float fade = x->fade;
    long totalframes;
    long swapframes;
    long r1startframe;
    long r1endframe;
    long start_sample1;
    t_float *storage = x->storage;
    t_word *b_samples;
    t_float mix_sample;
    t_float fadein_gain;
    t_float fadeout_gain;
    //  t_float fadephase;
    long ub1, lb1;
    long ub2, lb2;
    long block1size,block2size;
    long b_nchans;
    long syncframe;
    int i,j,k;
    buffet_setbuf(x, x->wavename);
    if(x->hosed) {
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(totalframes < maxframes * b_nchans + 1) {
        pd_error(0, "buffer must contain at least twice as many samples as the maximum swap size");
        return;
    }
    if(b_nchans > 2) {
        pd_error(0, "buffet~ only accepts mono or stereo buffers");
        return;
    }
    syncframe = x->sync * totalframes;
    swapframes = buffet_boundrand(minframes, maxframes);
    if(x->sync <= 0.0 ){ // either no sync signal or it is zero
        r1startframe = buffet_boundrand(0.0, (t_float)(totalframes-swapframes));
        r1endframe = r1startframe + swapframes;
    } else {
        // avoid swapping where we are playing from buffer
        lb1 = 0;
        ub1 = syncframe - CUSHION_FRAMES; // could be reading buffer backwards
        lb2 = ub1 + CUSHION_FRAMES;
        ub2 = totalframes - 1;
        block1size = ub1 - lb1;
        block2size = ub2 - lb2;
        
        if(block1size < maxframes && block2size < maxframes ) {
            pd_error(0, "could not reverse block");
            return;
        }
        if(block1size >= maxframes && block2size >= maxframes) {
            if(buffet_boundrand(0.0,1.0) > 0.5) {
                r1startframe = buffet_boundrand(0.0, (t_float)(ub1-swapframes));
            } else {
                r1startframe = buffet_boundrand((t_float)lb2, (t_float)(ub2-swapframes));
            }
        } else if(block1size < maxframes) {
            r1startframe = buffet_boundrand((t_float)lb2, (t_float)(ub2-swapframes));
        } else {
            r1startframe = buffet_boundrand(0.0, (t_float)(ub1-swapframes));
        }
        r1endframe = r1startframe + swapframes;
    }
    
    
    if(r1endframe >= totalframes) {
        pd_error(0, "%s: retro beyond bounds",OBJECT_NAME);
        return;
    }
    
    // store block1 samples
    start_sample1 = r1startframe * b_nchans;
    // start_sample2 = r2startframe * b_nchans;
    // store block in reversed order
    for(k = 0, i = (swapframes-1) * b_nchans; i > 0; i -= b_nchans, k += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            storage[i+j] = b_samples[start_sample1 + k + j].w_float;
        }
    }
    // swap block2 into block1 location, fadein
    
    for(i = 0, k = 0; i < fade * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadein_gain * storage[i + j] + fadeout_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    // middle part, pure swap
    
    
    for(i = fade * b_nchans; i < (swapframes-fade) * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample1 + i + j].w_float = storage[i + j];
        }
    }
    // fade out
    
    for(i = (swapframes-fade) * b_nchans, k = 0; i < swapframes * b_nchans; i += b_nchans, k++) {
        fadein_gain = (t_float)k / (t_float) fade;
        fadeout_gain = 1.0 - fadein_gain;
        for(j = 0; j < b_nchans; j++) {
            mix_sample = fadeout_gain * storage[i + j] + fadein_gain * b_samples[start_sample1 + i + j].w_float;
            b_samples[start_sample1 + i + j].w_float = mix_sample;
        }
    }
    buffet_update(x);
}

// clicky version
void buffet_nakedswap(t_buffet *x)
{
    //long
    long minframes = x->minframes;
    long maxframes = x->maxframes;
    long totalframes;
    long swapframes;
    long r1startframe;
    long r1endframe;
    long r2startframe;
    long r2endframe;
    long region1;
    long region2;
    long start_sample1;
    long start_sample2;
    t_float *storage = x->storage;
    t_word *b_samples;
    long b_nchans;
    int i,j;
    buffet_setbuf(x, x->wavename);
    
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    totalframes = x->wavebuf->b_frames;
    
    if(totalframes < maxframes * 2 + 1) {
        pd_error(0, "buffer must contain at least twice as many samples as the maximum swap size");
        return;
    }
    if(b_nchans != 2) {
        pd_error(0, "buffet~ only accepts stereo buffers");
        return;
    }
    swapframes = buffet_boundrand((t_float)minframes, (t_float)maxframes);
    r1startframe = buffet_boundrand(0.0, (t_float)(totalframes-swapframes));
    r1endframe = r1startframe + swapframes;
    region1 = r1startframe;
    region2 = totalframes - r1endframe;
    
    if(swapframes > region1) {
        r2startframe = buffet_boundrand((t_float)r1endframe,(t_float)(totalframes-swapframes));
    } else if(swapframes > region2) {
        r2startframe = buffet_boundrand(0.0,(t_float)(r1startframe-swapframes));
    } else {
        if(buffet_boundrand(0.0,1.0) > 0.5) {
            r2startframe = buffet_boundrand(0.0,(t_float)(r1startframe-swapframes));
        } else {
            r2startframe = buffet_boundrand((t_float)r1endframe,(t_float)(totalframes-swapframes));
        }
    }
    r2endframe = r2startframe + swapframes;
    if(r2startframe < 0 || r1startframe < 0) {
        pd_error(0, "start frame less than zero!");
        return;
    }
    if(r2endframe >= totalframes || r1endframe >= totalframes) {
        pd_error(0, "end frame reads beyond buffer!");
        return;
    }
    // Now swap the samples. For now no cross fade
    start_sample1 = r1startframe * b_nchans;
    start_sample2 = r2startframe * b_nchans;
    for(i = 0; i < swapframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            storage[i+j] = b_samples[start_sample1 + i + j].w_float;
        }
    }
    // swap block1 into block2
    for(i = 0; i < swapframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample1 + i + j].w_float = b_samples[start_sample2 + i + j].w_float;
        }
    }
    // swap stored block into block1
    for(i = 0; i < swapframes * b_nchans; i += b_nchans) {
        for(j = 0; j < b_nchans; j++) {
            b_samples[start_sample2 + i + j].w_float = storage[i + j];
        }
    }
    buffet_update(x);
}

void buffet_setbuf(t_buffet *x, t_symbol *wavename)
{
    
    t_garray *a;
    int b_frames;
    t_word *b_samples;
    if (!(a = (t_garray *)pd_findbyclass(wavename, garray_class))) {
        if (*wavename->s_name) pd_error(x, "%s: %s: no such array",OBJECT_NAME,wavename->s_name);
        x->wavebuf->b_samples = 0;
        x->wavebuf->b_valid = 0;
    }
    else if (!garray_getfloatwords(a, &b_frames, &b_samples)) {
        pd_error(x, "%s: bad array for %s", wavename->s_name,OBJECT_NAME);
        x->wavebuf->b_valid = 0;
    }
    else  {
        x->wavebuf->b_valid = 1;
        x->wavebuf->b_frames = b_frames;
        x->wavebuf->b_nchans = 1;
        x->wavebuf->b_samples = b_samples;
        garray_usedindsp(a);
    }
    
}


t_int *buffet_perform(t_int *w)
{
    t_buffet *x = (t_buffet *) (w[1]);
    t_float *sync = (t_float *)(w[2]);
    int n = (int) w[3];
    
    while(n--) {
        x->sync = *sync++;
    }
    
    return (w+4);
}


t_float buffet_boundrand(t_float min, t_float max)
{
    return min + (max-min) * ((t_float) (rand() % RAND_MAX)/ (t_float) RAND_MAX);
}

void buffet_detect_onsets(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    t_atom *listdata = x->listdata;
    long startframe;
    long currentframe;
    long endframe;
    long sample_frames;
    long fft_frames;
    int R = x->sr;
    int N = 1024;
    int N2 = N / 2;
    int Nw = N;
    int Nw2 = Nw / 2;
    int D = N / 2;
    int i,j,k,m;
    //  t_float sb1,sb2,sb3,sb4,energy=0,lastenergy;
    t_float freqcenter=0;
    t_float tmp;
    t_float fave; // freq average
    short first = 1;
    t_float dsum;
    short suppress = 0;
    t_float suppressrt = 0;
    
    t_float rt, rtadv;
    int inCount;
    t_float *Wanal;
    t_float *Wsyn;
    t_float *input;
    t_float *Hwin;
    t_float *buffer;
    t_float *channel;
    t_float *output;
    t_float *input_vec;
    t_float *specdiff;
    t_float *freqdiff;
    t_float *hfc;
    t_float *freqs, *amps;
    t_float lastfreqs[5];
    t_float fdiffs[5];
    t_float *onsets;
    t_float **loveboat;
    t_float *c_lastphase_in;
    t_float *c_lastphase_out;
    t_float c_fundamental;
    t_float c_factor_in;
    t_float *trigland;
    int *bitshuffle;
    int topanalbin;
    t_float threshold;
    int bytesize;
    
    int MAX_ONSETS = 2048;
    int onset_count = 0;
    
    buffet_setbuf(x, x->wavename);
    if( x->hosed ) {
        pd_error(0, "buffet~ needs a valid buffer");
        return;
    }
    if(argc < 4) {
        post("%s: detect_onsets requires start and end times,threshold and FFTsize",OBJECT_NAME);
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    threshold = atom_getfloatarg(2,argc,argv);
    N = atom_getfloatarg(3,argc,argv);
    N2 = N / 2;
    Nw = N;
    Nw2 = Nw / 2;
    D = N / 2;
    
    
    sample_frames = endframe - startframe;
    fft_frames = sample_frames / D;
    
    c_fundamental = (t_float)R / (t_float)N;
    inCount = -Nw;
    c_factor_in = (t_float)R/((t_float)D * TWOPI);
    
    rt = startframe * 1000.0 / x->sr;
    rtadv = 1000. * (t_float)D/(t_float)R;
    
    post("we will analyze %d FFT frames",fft_frames);
    
    Wanal = (t_float *) getbytes(Nw * sizeof(t_float));
    Wsyn = (t_float *) getbytes(Nw * sizeof(t_float));
    Hwin = (t_float *) getbytes(Nw * sizeof(t_float));
    input = (t_float *) getbytes(Nw * sizeof(t_float) );
    output = (t_float *) getbytes(Nw * sizeof(t_float) );
    buffer = (t_float *) getbytes(N * sizeof(t_float));
    channel = (t_float *) getbytes((N+2) * sizeof(t_float) );
    bitshuffle = (int *) getbytes(N * 2 * sizeof(int));
    bytesize = N * 2 * sizeof(t_float);
    //  post("allocing %d bytes to trigland",N * 2 * sizeof(t_float));
    
    trigland = (t_float *) getbytes(N * 2 * sizeof(t_float));
    //    post("cleaning %d bytes in trigland",N * 2 * sizeof(t_float));
    
    
    c_lastphase_in = (t_float *) getbytes((N2+1) * sizeof(t_float));
    c_lastphase_out = (t_float *) getbytes((N2+1) * sizeof(t_float));
    
    input_vec = (t_float *) getbytes(D * sizeof(t_float));
    onsets = (t_float *) getbytes(MAX_ONSETS * sizeof(t_float));
    specdiff = (t_float *) getbytes(fft_frames * sizeof(t_float));
    freqdiff = (t_float *) getbytes(fft_frames * sizeof(t_float));
    hfc = (t_float *) getbytes(fft_frames * sizeof(t_float));
    loveboat = (t_float **) getbytes(fft_frames * sizeof(t_float *));
    
    freqs = (t_float *) getbytes(N * sizeof(t_float));
    amps = (t_float *) getbytes(N * sizeof(t_float));
    for(i = 0; i < fft_frames; i++) {
        loveboat[i] = (t_float *) getbytes((N+2) * sizeof(t_float));
        if(loveboat[i] == NULL) {
            pd_error(0, "memory error");
            return;
        }
        // memset((char *)loveboat[i],0,(N+2)*sizeof(t_float));
    }
    
    memset((char *)trigland,0, N * 2 * sizeof(t_float));
    memset((char *)input,0,Nw * sizeof(t_float));
    memset((char *)output,0,Nw * sizeof(t_float));
    memset((char *)c_lastphase_in,0,(N2+1) * sizeof(t_float));
    memset((char *)c_lastphase_out,0,(N2+1) * sizeof(t_float));
    memset((char *)bitshuffle,0, 2 * N * sizeof(int));
    
    
    lpp_makewindows(Hwin, Wanal, Wsyn, Nw, N, D);
    
    
    lpp_init_rdft(N, bitshuffle, trigland);
    
    currentframe = startframe;
    
    for(k = 0; k < fft_frames; k++) {
        for(j = 0; j < D; j++) {
            input_vec[j] = b_samples[(currentframe + j) * b_nchans].w_float;// only read first channel
        }
        currentframe += D;
        
        for ( j = 0 ; j < Nw - D ; j++ ) {
            input[j] = input[j+D];
        }
        for (i = 0, j = Nw - D; j < Nw; j++, i++) {
            input[j] = input_vec[i];
        }
        lpp_fold(input, Wanal, Nw, buffer, N, inCount);
        lpp_rdft(N, 1, buffer, bitshuffle, trigland);
        lpp_convert(buffer, loveboat[k], N2, c_lastphase_in, c_fundamental, c_factor_in);
    }
    
    topanalbin = 11000.0 / c_fundamental;
    //  post("top bin %d",topanalbin);
    
    for(j = 0; j < 5; j++)
    lastfreqs[j] = 0;
    
    for(i = 0, k = 1; i < fft_frames - 1; i++, k++) {
        
        for(j = 0; j < 5; j++)
        lastfreqs[j] = freqs[j];
        
        for(j = 2, m = 0; j < topanalbin * 2; j += 2, m++){  // skip lowest bin
            freqcenter += fabs(loveboat[i][j]) * fabs(loveboat[i][j+1]);
            amps[m] = fabs(loveboat[i][j]);
            freqs[m] = loveboat[i][j+1];
        }
        /* swap top freqs to the top - but note there is a problem if some lower partials
         vacillate over which is the weakest since it will swap places and generate spurious
         change of frequency info. */
        
        for(j = 0; j < topanalbin-1; j++) {
            for(m = j + 1 ; m < topanalbin; m++) {
                if(amps[m]> amps[j]) {
                    tmp = amps[j];
                    amps[j] = amps[m];
                    amps[m] = tmp;
                    tmp = freqs[j];
                    freqs[j] = freqs[m];
                    freqs[m] = tmp;
                }
            }
        }
        
        if(first) {
            first = 0;
            for(j = 0; j < 5; j++) {
                lastfreqs[j] = freqs[j];
            }
        }
        for(j = 0; j < 5; j++) {
            fave = (lastfreqs[j] + freqs[j]) * 0.5;
            if(fave > 0)
                fdiffs[j] = fabs(lastfreqs[j] - freqs[j])/fave;
            else
                fdiffs[j] = 0;
        }
        dsum = fdiffs[0] + .5 * fdiffs[1] + .25 * fdiffs[2] + .125 * fdiffs[3] + .0625 * fdiffs[4];
        if(dsum > threshold && ! suppress) {
            suppress = 1;
            suppressrt = 0;
            onsets[onset_count++] = rt;
            /*
             post("%.2f: f1 %.2f d1 %.2f f2 %.2f d2 %.2f f3 %.2f d3 %.2f f4 %.2f d4 %.2f f5 %.2f d5 %.2f, dsum: %f",
             rt, freqs[0],fdiffs[0], freqs[1],fdiffs[1], freqs[2],fdiffs[2], freqs[3],fdiffs[3], freqs[4],fdiffs[4], dsum);
             */
            
        } else {
            suppressrt += rtadv;
            if(suppressrt > 100){ // 100 ms latency
                suppress = 0;
            }
            
        }
        rt += rtadv;
    }
    
    post("%s: %d onsets detected",OBJECT_NAME,onset_count);
    for(i = 0; i < onset_count; i++) {
        SETFLOAT(listdata + i, onsets[i]);
    }
    
    outlet_list(x->list, 0, onset_count, listdata);
    outlet_bang(x->bang);
    // free memory at end of routine
}

void buffet_detect_subband_onsets(t_buffet *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_word *b_samples;
    long b_nchans;
    long b_frames;
    long startframe;
    long currentframe;
    long endframe;
    long sample_frames;
    long fft_frames;
    int R = x->sr;
    int N = 1024;
    int N2 = N / 2;
    int Nw = N;
    int Nw2 = Nw / 2;
    int D = N / 2;
    int i,j,k,l;
    int minbin,maxbin;
    
    
    t_float rt, rtadv;
    int inCount;
    t_float *Wanal;
    t_float *Wsyn;
    t_float *input;
    t_float *Hwin;
    t_float *buffer;
 //   t_float *channel;
    t_float *output;
    t_float *input_vec;
 //   t_float *specdiff;
 //   t_float *freqdiff;
    t_float *hfc;
    t_float *freqs, *amps;
    
    t_float subband_energy[5];
    t_float subband_freqdiff[5];
    int subband_bincuts[6];
    
    t_float *onsets;
    t_float **loveboat;
    
    t_float *c_lastphase_in;
    t_float *c_lastphase_out;
    t_float c_fundamental;
    t_float c_factor_in;
    
    t_float *trigland;
    int *bitshuffle;
    
    
    t_float threshold;
    
    
    int MAX_ONSETS = 2048;
    
    buffet_setbuf(x, x->wavename);
    if( x->hosed ) {
        pd_error(0, "buffet~ needs a valid buffer");
        return;
    }
    if(argc < 4) {
        post("%s: detect_onsets requires start and end times,threshold and FFTsize",OBJECT_NAME);
        return;
    }
    b_samples = x->wavebuf->b_samples;
    b_nchans = x->wavebuf->b_nchans;
    b_frames = x->wavebuf->b_frames;
    startframe = .001 * x->sr * atom_getfloatarg(0,argc,argv);
    endframe = .001 * x->sr * atom_getfloatarg(1,argc,argv);
    threshold = atom_getfloatarg(2,argc,argv);
    N = atom_getfloatarg(3,argc,argv);
    N2 = N / 2;
    Nw = N;
    Nw2 = Nw / 2;
    D = N / 2;
    
    if(!endframe)
        endframe = b_frames;
    
    sample_frames = endframe - startframe;
    fft_frames = sample_frames / D;
    
    c_fundamental = (t_float)R / (t_float)N;
    inCount = -Nw;
    c_factor_in = (t_float)R/((t_float)D * TWOPI);
    
    rt = startframe * 1000.0 / x->sr;
    rtadv = 1000. * (t_float)D/(t_float)R;
    
    post("we will analyze %d FFT frames",fft_frames);
    subband_bincuts[0] = 1;
    subband_bincuts[1] = (int) (500.0 / c_fundamental);
    subband_bincuts[2] = (int) (1000.0 / c_fundamental);
    subband_bincuts[3] = (int) (2000.0 / c_fundamental);
    subband_bincuts[4] = (int) (4000.0 / c_fundamental);
    subband_bincuts[5] = (int) (8000.0 / c_fundamental);
    
    Wanal = (t_float *) getbytes(Nw * sizeof(t_float));
    Wsyn = (t_float *) getbytes(Nw * sizeof(t_float));
    Hwin = (t_float *) getbytes(Nw * sizeof(t_float));
    input = (t_float *) getbytes(Nw * sizeof(t_float) );
    output = (t_float *) getbytes(Nw * sizeof(t_float) );
    buffer = (t_float *) getbytes(N * sizeof(t_float));
//    channel = (t_float *) getbytes((N+2) * sizeof(t_float) );
    bitshuffle = (int *) getbytes(N * 2 * sizeof(int));
    trigland = (t_float *) getbytes(N * 2 * sizeof(t_float));
    c_lastphase_in = (t_float *) getbytes((N2+1) * sizeof(t_float));
    c_lastphase_out = (t_float *) getbytes((N2+1) * sizeof(t_float));
    
    input_vec = (t_float *) getbytes(D * sizeof(t_float));
    onsets = (t_float *) getbytes(MAX_ONSETS * sizeof(t_float));
//    specdiff = (t_float *) getbytes(fft_frames * sizeof(t_float));
//    freqdiff = (t_float *) getbytes(fft_frames * sizeof(t_float));
    hfc = (t_float *) getbytes(fft_frames * sizeof(t_float));
    loveboat = (t_float **) getbytes(fft_frames * sizeof(t_float *));
    
    freqs = (t_float *) getbytes(N * sizeof(t_float));
    amps = (t_float *) getbytes(N * sizeof(t_float));
    for(i = 0; i < fft_frames; i++) {
        loveboat[i] = (t_float *) getbytes((N+2) * sizeof(t_float));
        if(loveboat[i] == NULL) {
            pd_error(0, "memory error");
            return;
        }
        // memset((char *)loveboat[i],0,(N+2)*sizeof(t_float));
    }
    /*
    memset((char *)trigland,0, N * 2 * sizeof(t_float));
    memset((char *)input,0,Nw * sizeof(t_float));
    memset((char *)output,0,Nw * sizeof(t_float));
    memset((char *)c_lastphase_in,0,(N2+1) * sizeof(t_float));
    memset((char *)c_lastphase_out,0,(N2+1) * sizeof(t_float));
    memset((char *)bitshuffle,0, 2 * N * sizeof(int));
    */
    
    lpp_makewindows(Hwin, Wanal, Wsyn, Nw, N, D);
    
    
    lpp_init_rdft(N, bitshuffle, trigland);
    
    currentframe = startframe;
    
    for(k= 0; k < 5; k++) {
        minbin = subband_bincuts[k]; maxbin = subband_bincuts[k+1];
        post("%d: minbin %d maxbin %d",k, minbin,maxbin);
    }
    
    for(k = 0; k < fft_frames; k++) {
        for(j = 0; j < D; j++) {
            input_vec[j] = b_samples[(currentframe + j) * b_nchans].w_float;// only read first channel
        }
        currentframe += D;
        
        for ( j = 0 ; j < Nw - D ; j++ ) {
            input[j] = input[j+D];
        }
        for (i = 0, j = Nw - D; j < Nw; j++, i++) {
            input[j] = input_vec[i];
        }
        lpp_fold(input, Wanal, Nw, buffer, N, inCount);
        lpp_rdft(N, 1, buffer, bitshuffle, trigland);
        lpp_convert(buffer, loveboat[k], N2, c_lastphase_in, c_fundamental, c_factor_in);
    }
    
    rt = rtadv;
    for(i = 0, j = 1; i < fft_frames - 1; i++, j++) {
        
        
        for(k= 0; k < 5; k++) {
            minbin = subband_bincuts[k]; maxbin = subband_bincuts[k+1];
            subband_energy[k] = 0;
            subband_freqdiff[k] = 0;
            /* maxamp = 0;
             for(l = minbin * 2; l < maxbin * 2; l += 2) {
             } */
            for(l = minbin * 2; l < maxbin * 2; l += 2) {
                //    subband_energy[k] += (loveboat[j][l] * loveboat[j][l])  - (loveboat[i][l] * loveboat[i][l]);
                subband_freqdiff[k] += (fabs(loveboat[j][l+1] - loveboat[i][l+1])/loveboat[i][l+1]) *
                loveboat[j][l] * loveboat[j][l] ;// scale by energy
            }
        }
        /*
         post("rt: %.2f sb1 %f sb2 %f sb3 %f sb4 %f sb5 %f",
         rt, subband_energy[0],subband_energy[1],subband_energy[2],subband_energy[3],subband_energy[4]
         );
         */
        post("rt:%f fb1 %f fb2 %f fb3 %f fb4 %f fb5 %f",
             rt, subband_freqdiff[0],subband_freqdiff[1],subband_freqdiff[2],subband_freqdiff[3],subband_freqdiff[4]
             );
        /*
         post("%.2f: f1 %.2f d1 %.2f f2 %.2f d2 %.2f f3 %.2f d3 %.2f f4 %.2f d4 %.2f f5 %.2f d5 %.2f, dsum: %f",
         rt, freqs[0],fdiffs[0], freqs[1],fdiffs[1], freqs[2],fdiffs[2], freqs[3],fdiffs[3], freqs[4],fdiffs[4], dsum);
         */
        
        
        rt += rtadv;
    }
    /*
     post("%s: %d onsets detected",OBJECT_NAME,onset_count);
     for(i = 0; i < onset_count; i++) {
     SETFLOAT(listdata + i, onsets[i]);
     }
     
     outlet_list(x->list, 0, onset_count, listdata);
     
     outlet_bang(x->bang);
     */
    // free memory at end of routine
}


void buffet_dsp_free(t_buffet *x)
{
    
    freebytes(x->storage,0);
    freebytes(x->listdata,0);
    freebytes(x->rmsbuf,0);
    freebytes(x->analbuf,0);
    freebytes(x->onset,0);
    freebytes(x->wavebuf,0);
    freebytes(x->destbuf,0);
}

void buffet_dsp(t_buffet *x, t_signal **sp)
{
    
    buffet_setbuf(x, x->wavename);
    
    if( x->hosed ) {
        pd_error(0, "buffet~ needs a valid buffer");
        return;
    }
    if( x->sr != sp[0]->s_sr) {
        x->sr = sp[0]->s_sr;
        
        if(!x->sr) {
            post("%s: warning: zero sampling rate!",OBJECT_NAME);
            x->sr = 44100;
        }
        buffet_init(x,1);
    }
    
    dsp_add(buffet_perform, 3, x, sp[0]->s_vec, (t_int)sp[0]->s_n);
    
}
