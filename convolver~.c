#include "m_pd.h"
#include "fftease.h"
#include "MSPd.h"

#define OBJECT_NAME "convolver~"
#define DENORM_WANT_FIX   1

#define COMPILE_DATE "9.4.22"
#define OBJECT_VERSION "2.03"

#define FIX_DENORM_FLOAT(v) (v=(fabs(v) < 0.000001 ? 0.f : (v)))
static t_class *convolver_class;

#define CBUF_SIZE 32768
#define NCMAX 52428800

typedef struct _buffy {
    t_word *b_samples;
    long b_frames;
    long b_nchans;
    long b_valid;
    t_garray *b;
    t_symbol *myname;
} t_buffy;


typedef struct _convolver
{
    t_object x_obj;
    t_float x_f;
    t_buffy *impulse; // impulse buffer
    t_buffy *source; // source buffer
    t_buffy *dest; // output buffer
    void *bang; // completion bang
    t_float sr;
    // convolution stuff
    t_float *tbuf;
    t_float *sbuf;
    t_float *filt;
    long N;
    long N2;
    long last_N;
    // for fast fft
    t_float mult;
    t_float *trigland;
    int *bitshuffle;
    short static_memory; // flag to avoid dynamic memory manipulation
} t_convolver;


static t_float boundrand(t_float min, t_float max);
static void convolver_setbuf(t_buffy *trybuf);
static void *convolver_new(t_symbol *msg, int argc, t_atom *argv);
static t_int *convolver_perform(t_int *w);
static void convolver_dsp_free(t_convolver *x);
static void convolver_attach_buffers(t_convolver *x) ;
static void convolver_spikeimp(t_convolver *x, t_floatarg density);
static void convolver_convolve(t_convolver *x);
static void convolver_convolvechans(t_convolver *x, t_symbol *msg, int argc, t_atom *argv);
static void convolver_noiseimp(t_convolver *x, t_floatarg curve);
static void rfft( t_float *x, int N, int forward );
static void cfft( t_float *x, int NC, int forward );
static void rdft(int n, int isgn, t_float *a, int *ip, t_float *w);
static void bitreverse( t_float *x, int N );
static void init_rdft(int n, int *ip, t_float *w);
static void convolver_static_memory(t_convolver *x, t_floatarg toggle);

void convolver_tilde_setup(void) {
    convolver_class = class_new(gensym("convolver~"), (t_newmethod)convolver_new,
                                (t_method)convolver_dsp_free,sizeof(t_convolver), 0,A_GIMME,0);
    CLASS_MAINSIGNALIN(convolver_class, t_convolver, x_f);
    class_addmethod(convolver_class,(t_method)convolver_spikeimp,gensym("spikeimp"), A_FLOAT, 0);
    class_addmethod(convolver_class,(t_method)convolver_noiseimp, gensym("noiseimp"),A_FLOAT, 0);
    class_addmethod(convolver_class,(t_method)convolver_convolve, gensym("convolve"),0);
    class_addmethod(convolver_class,(t_method)convolver_static_memory, gensym("static_memory"),A_FLOAT, 0);
    potpourri_announce(OBJECT_NAME);
}

void convolver_static_memory(t_convolver *x, t_floatarg toggle)
{
    
    long memcount = 0;
    t_float *tbuf = x->tbuf;
    t_float *sbuf = x->sbuf;
    t_float *filt = x->filt;
    int *bitshuffle = x->bitshuffle;
    t_float *trigland = x->trigland;
    t_buffy *impulse = x->impulse;
    long N, N2;
    
    x->static_memory = (short) toggle;
    
    if( x->static_memory ) {
        convolver_attach_buffers( x );
        
        for( N2 = 2; N2 < NCMAX; N2 *= 2) {
            if( N2 >= impulse->b_frames ) {
                // post("%s: Exceeded Impulse Maximum: %d",OBJECT_NAME, NCMAX);
                break;
            }
        }
        N = 2 * N2;
        
        post("%s: memory is now static - do not reload your impulse buffer",OBJECT_NAME);
        
        if ((sbuf = (t_float *) getbytes((N+2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory", OBJECT_NAME);
        memcount += (N+2) * sizeof(t_float);
        if ((tbuf = (t_float *) getbytes(N2 * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        memcount += (N2) * sizeof(t_float);
        if ((filt = (t_float *) getbytes((N+2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        memcount += (N+2) * sizeof(t_float);
        if( (bitshuffle = (int *) getbytes((N*2) * sizeof(int))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        memcount += (N2) * sizeof(t_float);
        if( (trigland = (t_float *) getbytes((N*2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        memcount += (N2) * sizeof(t_float);
        post("%s: allocated %f Megabytes for %s", OBJECT_NAME, (t_float)memcount / 1000000.0, impulse->myname->s_name);
        x->N = N;
        x->N2 = N2;
    }
}


/*
void convolver_seed(t_convolver *x, t_floatarg seed)
{
    (void)x;
    srand((long)seed);
}
*/

void convolver_convolve(t_convolver *x)
{
    int i;
    t_symbol *mymsg;
    short myargc = 3;
    t_atom data[3];
    mymsg = (t_symbol *) getbytes(1 * sizeof(t_symbol));
    convolver_attach_buffers( x );
    if(x->source->b_nchans == x->impulse->b_nchans && x->impulse->b_nchans == x->dest->b_nchans) {
        //post("case 1");
        for(i = 0; i < x->source->b_nchans; i++) {
            SETFLOAT(data, i+1); // source
            SETFLOAT(data+1, i+1); // impulse
            SETFLOAT(data+2, i+1); // destination
            convolver_convolvechans(x, mymsg, myargc, data);
        }
    }
    else if(x->source->b_nchans == 1 && x->impulse->b_nchans == x->dest->b_nchans) {
        //post("case 2");
        for(i = 0; i < x->impulse->b_nchans; i++) {
            SETFLOAT(data, 1); // source
            SETFLOAT(data+1, i+1); // impulse
            SETFLOAT(data+2, i+1); // destination
            convolver_convolvechans(x, mymsg, myargc, data);
        }
    }
    else if(x->impulse->b_nchans == 1 && x->source->b_nchans == x->dest->b_nchans) {
        //post("case 3");
        for(i = 0; i < x->impulse->b_nchans; i++) {
            SETFLOAT(data, i+1); // source
            SETFLOAT(data+1, 1); // impulse
            SETFLOAT(data+2, i+1); // destination
            convolver_convolvechans(x, mymsg, myargc, data);
        }
    } else {
        post("%s: \"convolve\" is not smart enough to figure out what you want to do. Try \"convolvechans\"",OBJECT_NAME);
        post("source chans: %d, impulse chans: %d, dest chans: %d",x->source->b_nchans, x->impulse->b_nchans, x->dest->b_nchans );
    }
    outlet_bang(x->bang);
}

void convolver_convolvechans(t_convolver *x, t_symbol *msg, int argc, t_atom *argv)
{
    t_float *tbuf = x->tbuf;
    t_float *sbuf = x->sbuf;
    t_float *filt = x->filt;
    long N = x->N;
    long N2 = x->N2;
    long i, j, ip, ip1;
    long ifr_cnt = 0, ofr_cnt = 0;
    int target_frames = 2;
    short copacetic; // loop enabler
    t_float a,b,temp,max=0.0,gain=1.0; //,thresh=.0000000001,fmag;
    int readframes, writeframes;
    t_buffy *impulse = x->impulse;
    t_buffy *source = x->source;
    t_buffy *dest = x->dest;
    int *bitshuffle = x->bitshuffle;
    t_float *trigland = x->trigland;
    long source_chan, impulse_chan, dest_chan;
    t_float rescale = 0.5 / (t_float) N;
    (void)msg;
    //  t_atom newsize;
    
    convolver_attach_buffers( x );
    
    source_chan = atom_getfloatarg(0,argc,argv);
    impulse_chan = atom_getfloatarg(1,argc,argv);
    dest_chan = atom_getfloatarg(2,argc,argv);
    // post("chans %d %d %d", source_chan, impulse_chan, dest_chan);
    if( source_chan <= 0 || impulse_chan <= 0 || dest_chan <= 0) {
        pd_error(0, "%s: channels are counted starting from 1",OBJECT_NAME);
        return;
    }
    if( source_chan > source->b_nchans ) {
        pd_error(0, "%s: source channel %ld out of range", OBJECT_NAME, source_chan);
        return;
    }
    if( impulse_chan > impulse->b_nchans ) {
        pd_error(0, "%s: impulse channel %ld out of range", OBJECT_NAME, impulse_chan);
        return;
    }
    if( dest_chan > dest->b_nchans ) {
        pd_error(0, "%s: dest channel %ld out of range", OBJECT_NAME, dest_chan);
        return;
    }
    --source_chan;
    --impulse_chan;
    --dest_chan;
    //inframes = source->b_frames;
    //outframes = dest->b_frames;
    // initialization routine (move out and only do once)
    for( N2 = 2; N2 < NCMAX; N2 *= 2) {
        if( N2 >= impulse->b_frames )
            break;
    }
    N = 2 * N2;
    // be more careful with memory
    // also be sure to clear destination buffer
    
    // post("size of N for convolution is %d", N);
    if(! x->static_memory ) {
        if ((sbuf = (t_float *) getbytes((N+2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory", OBJECT_NAME);
        if ((tbuf = (t_float *) getbytes(N2 * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        if ((filt = (t_float *) getbytes((N+2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        if( (bitshuffle = (int *) getbytes((N*2) * sizeof(int))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
        if( (trigland = (t_float *) getbytes((N*2) * sizeof(t_float))) == NULL)
            pd_error(0, "%s: insufficient memory",OBJECT_NAME);
    }
    
    x->mult = 1. / (t_float) N;
    x->last_N = N;
    init_rdft(N, bitshuffle, trigland);
    
    for(i = 0, j = 0; i < impulse->b_frames; i+= impulse->b_nchans, j++) {
        filt[j] = impulse->b_samples[i + impulse_chan].w_float;
    }
    
    rdft( N, 1, filt, bitshuffle, trigland );
    
    for (i=0; i <= N; i += 2) {
        
        a = filt[i];
        b = filt[i + 1];
        temp = a*a + b*b;
        if (temp > max)
            max = temp;
    }
    
    if (max != 0.) {
        max = gain/(sqrt(max));
    }
    else {
        pd_error(0, "%s: impulse response is all zeros",OBJECT_NAME);
        return;
    }
    // make normalization optional
    for (i=0; i< N+2; i++)
    filt[i] *= max;
    
    ifr_cnt = ofr_cnt = 0;
    
    if(source->b_frames - ifr_cnt >= N2)
        readframes = N2;
    else readframes = source->b_frames - ifr_cnt;
    // read desired channel from multichannel source buffer into sbuf
    for(i = 0, j = ifr_cnt * dest->b_nchans; i < readframes; i++, ifr_cnt++, j += source->b_nchans)
    sbuf[i] = source->b_samples[j + source_chan].w_float;
    // zero pad source buffer
    for (i = readframes; i<N+2; i++)
    sbuf[i] = 0.;
    copacetic = 1;
    
    while( target_frames < source->b_frames + impulse->b_frames ) {
        target_frames *= 2;
    }
    //post("src frames + imp frames %d dest frames %d",source->b_frames + impulse->b_frames, dest->b_frames);
    
    if( dest->b_frames < target_frames) {
        
        
        //SETFLOAT(&newsize, (t_float) target_frames);
        // typedmess((void *) x->dest->b, gensym("sizeinsamps"),1, &newsize);
        garray_resize_long(x->dest->b,target_frames );
        post("%s: destination buffer was too small and has been resized",OBJECT_NAME);
        convolver_attach_buffers( x );
    }
    
    while(copacetic && ofr_cnt < ifr_cnt + N2) {
        // post("ofr %d ifr %d, N %d",ofr_cnt, ifr_cnt, N);
        
        // convolve source buffer with filter buffer
        rdft( N, 1, sbuf, bitshuffle, trigland );
        for (i=0; i<=N2; i++) {
            ip = 2*i;
            ip1 = ip + 1;
            a = sbuf[ip] * filt[ip] - sbuf[ip1] * filt[ip1];
            b = sbuf[ip] * filt[ip1] + sbuf[ip1] * filt[ip];
            sbuf[ip] = a;
            sbuf[ip1] = b;
        }
        // inverse fft
        
        rdft( N, -1, sbuf, bitshuffle, trigland );
        
        //accumulate to output buffer
        // denormals fix is in
        for (i=0; i<N2; i++) {
            FIX_DENORM_FLOAT(sbuf[i]);
            tbuf[i] += sbuf[i];
        }
        
        // write to msp buffer
        if(dest->b_frames - ofr_cnt >= N2)
            writeframes = N2;
        else {
            writeframes = source->b_frames - ofr_cnt;
            // post("cutting off with N2 %d dest frames - ofr %d", N2, dest->b_frames - ofr_cnt );
            copacetic = 0; // reached end of dest buffer
        }
        //shift samples to desired channel of multichannel output buffer
        for(i = 0, j = ofr_cnt * dest->b_nchans; i < writeframes; i++, ofr_cnt++, j += dest->b_nchans)
        dest->b_samples[j + dest_chan].w_float = tbuf[i];
        
        // shift over remaining convolved samples
        for (i=0; i<N2; i++) {
            FIX_DENORM_FLOAT(sbuf[N2 + i]);
            tbuf[i] = sbuf[N2+i];
        }
        // read in next batch
        if(source->b_frames - ifr_cnt >= N2)
            readframes = N2;
        else {
            readframes = source->b_frames - ifr_cnt;
        }
        // arbitrary input channel read
        for(i = 0, j = ifr_cnt * source->b_nchans; i < readframes; i++, ifr_cnt++, j += source->b_nchans)
        sbuf[i] = source->b_samples[j + source_chan].w_float;
        // zero pad
        for (i = readframes; i<N+2; i++)
        sbuf[i] = 0.;
    }
    // now normalize output buffer
    
    
    // OK
    //  post("first rescale: %f", rescale);
    max = 0.0;
    for(i = 0, j = 0; i < dest->b_frames; i++, j += dest->b_nchans) {
        if(max < fabs(dest->b_samples[j + dest_chan].w_float) )
            max = fabs(dest->b_samples[j + dest_chan].w_float);
    }
    if(max <= 0.0) {
        post("convolvesf: zero output");
        return;
    }
    rescale = 1.0 / max;
    // post("max: %f, second rescale: %f", max, rescale);
    
    for(i = 0, j = 0; i < dest->b_frames; i++, j+= dest->b_nchans) {
        dest->b_samples[j + dest_chan].w_float *= rescale;
    }
    // FAILED BY HERE
    // post("rescale done");
    //  return;
    
    if(! x->static_memory ) {
        freebytes(sbuf,(N+2) * sizeof(t_float));
        freebytes(tbuf,N2 * sizeof(t_float));
        freebytes(filt,(N+2) * sizeof(t_float));
        freebytes(bitshuffle,(N*2) * sizeof(int));
        freebytes(trigland,(N*2) * sizeof(t_float));
    } else {
        x->N = N;
        x->N2 =  N2;
    }
    outlet_bang(x->bang);
    garray_redraw(x->dest->b);

}

void convolver_noiseimp(t_convolver *x, t_floatarg curve)
{
    long b_frames;
    t_word *b_samples;
    t_float sr = x->sr;
    int i;
    int count;
    //  int position;
    t_float guess;
    t_float level = 1.0, endLevel = 0.001;
    t_float grow, a1, a2, b1;
    
    if(fabs(curve) < 0.001) {
        curve = 0.001;
    }
    // let's be current
    convolver_attach_buffers(x);
    //b_nchans = x->impulse->b_nchans;
    b_frames = x->impulse->b_frames;
    b_samples = x->impulse->b_samples;
    // chan test
    if( sr == 0. ) {
        pd_error(0, "zero sample rate");
        return;
    }
    // zero out buffer
    //dur = (t_float) b_frames / sr;
    count = b_frames;
    if(b_frames < 20) {
        post("impulse buffer too small!");
        return;
    }
    
    
    //  memset((char *)b_samples, 0, b_nchans * b_frames * sizeof(t_float));
    // return;
    
    level = 1.0;
    endLevel = 0.001;
    grow = exp(curve / (count - 1) );
    a1 = (endLevel - level) / (1.0 - exp(curve));
    a2 = level + a1;
    b1 = a1;
    for( i = 0; i < b_frames; i++ ) {
        guess = boundrand(-1.0, 1.0);
        //gain = 1. - guess;
        
        b1 = b1 * grow;
        level = a2 - b1;
        b_samples[i].w_float = level * guess;
    }
    
    garray_redraw(x->impulse->b);
    outlet_bang(x->bang);
    
}

void convolver_spikeimp(t_convolver *x, t_floatarg density)
{
    long b_nchans;
    long b_frames;
    t_word *b_samples;
    t_float sr = x->sr;
    int i, j;
    int count;
    int position;
    t_float gain, guess;
    t_float dur;
    
    // let's be current
    convolver_attach_buffers(x);
    b_nchans = x->impulse->b_nchans;
    b_frames = x->impulse->b_frames;
    b_samples = x->impulse->b_samples;
    // chan test
    if( sr == 0. ) {
        pd_error(0, "zero sample rate");
        return;
    }
    // zero out buffer
    dur = (t_float) b_frames / sr;
    count = density * dur;
    // memset((char *)b_samples, 0, b_nchans * b_frames * sizeof(t_float));
    // assume mono forever
    for(i = 0; i < b_frames; i++){
        b_samples[i].w_float = 0.0;
    }
    
    // return;
    for( j = 0; j < b_nchans; j++ ) {
        for( i = 0; i < count; i++ ) {
            guess = boundrand(0., 1.);
            gain = 1. - guess;
            gain = gain * gain;
            if( boundrand(0.0,1.0) > 0.5 ) {
                gain = gain * -1.0; // randomly invert signs to remove DC
            }
            position = (int) (dur * guess * guess * sr) * b_nchans + j;
            if( position >= b_frames * b_nchans ) {
                pd_error(0, "%d exceeds %ld",position, b_frames * b_nchans);
            } else{
                b_samples[ position ].w_float = gain;
            }
        }
    }
    garray_redraw(x->impulse->b);
    outlet_bang(x->bang);
}

t_float boundrand(t_float min, t_float max)
{
    return min + (max-min) * ((t_float) (rand() % RAND_MAX)/ (t_float) RAND_MAX);
}


void *convolver_new(t_symbol *msg, int argc, t_atom *argv)
{
    
    t_convolver *x = (t_convolver *)pd_new(convolver_class);
    (void)msg;
    x->bang = outlet_new(&x->x_obj, gensym("bang"));
    srand(time(0)); //need "seed" message
    x->impulse = (t_buffy *)getbytes(sizeof(t_buffy));
    x->source = (t_buffy *)getbytes(sizeof(t_buffy));
    x->dest = (t_buffy *)getbytes(sizeof(t_buffy));
    x->static_memory = 0;
    
    // default names
    x->impulse->myname = gensym("impulse_buf");
    x->source->myname = gensym("source_buf");
    x->dest->myname = gensym("dest_buf");
    
    x->last_N = -1;
    
    x->source->myname = atom_getsymbolarg(0,argc,argv);
    x->impulse->myname = atom_getsymbolarg(1,argc,argv);
    x->dest->myname = atom_getsymbolarg(2,argc,argv);
    
    x->sr = sys_getsr();
    return x;
}


void convolver_attach_buffers(t_convolver *x)
{
    convolver_setbuf(x->source);
    convolver_setbuf(x->impulse);
    convolver_setbuf(x->dest);
}


void convolver_setbuf(t_buffy *trybuf)
{
    t_garray *a;
    int b_frames;
    /* load up sample array */
    if (!(a = (t_garray *)pd_findbyclass(trybuf->myname, garray_class))) {
        if (*trybuf->myname->s_name) pd_error(0, "%s: %s: no such array", OBJECT_NAME, trybuf->myname->s_name);
    }
    else if (!garray_getfloatwords(a, &b_frames, &trybuf->b_samples)) { // possible crash worry?
        pd_error(0, "%s: bad template for %s", trybuf->myname->s_name,OBJECT_NAME);
        trybuf->b_valid = 0;
    }
    else  {
        trybuf->b_frames = b_frames;
        trybuf->b_nchans = 1; // Pd buffers are always mono (so far)
        trybuf->b = a; // link to array
        garray_usedindsp(a);
    }
}

/*
t_int *convolver_perform(t_int *w)
{
    return w + 4; // maybe we don't need this guy at all
}
*/

void convolver_dsp_free(t_convolver *x)
{
    int N = x->N;
    int N2 = x->N2;
    freebytes(x->impulse, sizeof(t_buffy));
    freebytes(x->source, sizeof(t_buffy));
    freebytes(x->dest, sizeof(t_buffy));
    if( x->static_memory ) {
        freebytes(x->sbuf, (N+2) * sizeof(t_float));
        freebytes(x->tbuf, N2 * sizeof(t_float));
        freebytes(x->filt, (N+2) * sizeof(t_float));
        freebytes(x->bitshuffle, (N*2) * sizeof(int));
        freebytes(x->trigland, (N*2) * sizeof(t_float));
       //  outlet_bang(x->bang);
    }
}

/*
 void convolver_dsp(t_convolver *x, t_signal **sp, short *count)
 {
 // never actually do anything here
 }
 */


// old FFT stuff, soon to be replaced




void cfft( t_float *x, int NC, int forward )

{
    t_float   wr,wi,
    wpr,wpi,
    theta,
    scale;
    int     mmax,
    ND,
    m,
    i,j,
    delta;
    
    // void bitreverse();
    
    ND = NC<<1;
    bitreverse( x, ND );
    for ( mmax = 2; mmax < ND; mmax = delta ) {
        delta = mmax<<1;
        theta = TWOPI/( forward? mmax : -mmax );
        wpr = -2.*pow( sin( 0.5*theta ), 2. );
        wpi = sin( theta );
        wr = 1.;
        wi = 0.;
        for ( m = 0; m < mmax; m += 2 ) {
            register t_float rtemp, itemp;
            for ( i = m; i < ND; i += delta ) {
                j = i + mmax;
                rtemp = wr*x[j] - wi*x[j+1];
                itemp = wr*x[j+1] + wi*x[j];
                x[j] = x[i] - rtemp;
                x[j+1] = x[i+1] - itemp;
                x[i] += rtemp;
                x[i+1] += itemp;
            }
            wr = (rtemp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + rtemp*wpi + wi;
        }
    }
    
    /* scale output */
    
    scale = forward ? 1./ND : 2.;
    { register t_float *xi=x, *xe=x+ND;
        while ( xi < xe )
            *xi++ *= scale;
    }
}

/* bitreverse places t_float array x containing N/2 complex values
 into bit-reversed order */

void bitreverse( t_float *x, int N )

{
    t_float   rtemp,itemp;
    int     i,j,
    m;
    
    for ( i = j = 0; i < N; i += 2, j += m ) {
        if ( j > i ) {
            rtemp = x[j]; itemp = x[j+1]; /* complex exchange */
            x[j] = x[i]; x[j+1] = x[i+1];
            x[i] = rtemp; x[i+1] = itemp;
        }
        for ( m = N>>1; m >= 2 && j >= m; m >>= 1 )
        j -= m;
    }
}

void init_rdft(int n, int *ip, t_float *w)
{
    
    int nw,
    nc;
    
    nw = n >> 2;
    lpp_makewt(nw, ip, w);
    
    nc = n >> 2;
    lpp_makect(nc, ip, w + nw);
    
    return;
}


void rdft(int n, int isgn, t_float *a, int *ip, t_float *w)
{
    
    int   j,
    nw,
    nc;
    
    t_float   xi;
    
    // void   bitrv2(int n, int *ip, t_float *a),
    //    cftsub(int n, t_float *a, t_float *w),
    //    rftsub(int n, t_float *a, int nc, t_float *c);
    
    
    nw = ip[0];
    nc = ip[1];
    
    if (isgn < 0) {
        a[1] = 0.5 * (a[1] - a[0]);
        a[0] += a[1];
        
        for (j = 3; j <= n - 1; j += 2) {
            a[j] = -a[j];
        }
        
        if (n > 4) {
            lpp_rftsub(n, a, nc, w + nw);
            lpp_bitrv2(n, ip + 2, a);
        }
        
        lpp_cftsub(n, a, w);
        
        for (j = 1; j <= n - 1; j += 2) {
            a[j] = -a[j];
        }
    }
    
    else {
        
        if (n > 4) {
            lpp_bitrv2(n, ip + 2, a);
        }
        
        lpp_cftsub(n, a, w);
        
        if (n > 4) {
            lpp_rftsub(n, a, nc, w + nw);
        }
        
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    }
}


void bitrv2(int n, int *ip, t_float *a)
{
    int j, jj1, k, k1, l, m, m2;
    t_float xr, xi;
    
    ip[0] = 0;
    l = n;
    m = 1;
    
    while ((m << 2) < l) {
        l >>= 1;
        for (j = 0; j <= m - 1; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    
    if ((m << 2) > l) {
        
        for (k = 1; k <= m - 1; k++) {
            
            for (j = 0; j <= k - 1; j++) {
                jj1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[jj1];
                xi = a[jj1 + 1];
                a[jj1] = a[k1];
                a[jj1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
    
    else {
        m2 = m << 1;
        
        for (k = 1; k <= m - 1; k++) {
            
            for (j = 0; j <= k - 1; j++) {
                jj1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[jj1];
                xi = a[jj1 + 1];
                a[jj1] = a[k1];
                a[jj1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
                jj1 += m2;
                k1 += m2;
                xr = a[jj1];
                xi = a[jj1 + 1];
                a[jj1] = a[k1];
                a[jj1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}


void cftsub(int n, t_float *a, t_float *w)
{
    int j, jj1, j2, j3, k, k1, ks, l, m;
    t_float wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    t_float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    l = 2;
    
    while ((l << 1) < n) {
        m = l << 2;
        
        for (j = 0; j <= l - 2; j += 2) {
            jj1 = j + l;
            j2 = jj1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[jj1];
            x0i = a[j + 1] + a[jj1 + 1];
            x1r = a[j] - a[jj1];
            x1i = a[j + 1] - a[jj1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[jj1] = x1r - x3i;
            a[jj1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
        
        if (m < n) {
            wk1r = w[2];
            
            for (j = m; j <= l + m - 2; j += 2) {
                jj1 = j + l;
                j2 = jj1 + l;
                j3 = j2 + l;
                x0r = a[j] + a[jj1];
                x0i = a[j + 1] + a[jj1 + 1];
                x1r = a[j] - a[jj1];
                x1i = a[j + 1] - a[jj1 + 1];
                x2r = a[j2] + a[j3];
                x2i = a[j2 + 1] + a[j3 + 1];
                x3r = a[j2] - a[j3];
                x3i = a[j2 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                a[j2] = x2i - x0i;
                a[j2 + 1] = x0r - x2r;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[jj1] = wk1r * (x0r - x0i);
                a[jj1 + 1] = wk1r * (x0r + x0i);
                x0r = x3i + x1r;
                x0i = x3r - x1i;
                a[j3] = wk1r * (x0i - x0r);
                a[j3 + 1] = wk1r * (x0i + x0r);
            }
            
            k1 = 1;
            ks = -1;
            
            for (k = (m << 1); k <= n - m; k += m) {
                k1++;
                ks = -ks;
                wk1r = w[k1 << 1];
                wk1i = w[(k1 << 1) + 1];
                wk2r = ks * w[k1];
                wk2i = w[k1 + ks];
                wk3r = wk1r - 2 * wk2i * wk1i;
                wk3i = 2 * wk2i * wk1r - wk1i;
                
                for (j = k; j <= l + k - 2; j += 2) {
                    jj1 = j + l;
                    j2 = jj1 + l;
                    j3 = j2 + l;
                    x0r = a[j] + a[jj1];
                    x0i = a[j + 1] + a[jj1 + 1];
                    x1r = a[j] - a[jj1];
                    x1i = a[j + 1] - a[jj1 + 1];
                    x2r = a[j2] + a[j3];
                    x2i = a[j2 + 1] + a[j3 + 1];
                    x3r = a[j2] - a[j3];
                    x3i = a[j2 + 1] - a[j3 + 1];
                    a[j] = x0r + x2r;
                    a[j + 1] = x0i + x2i;
                    x0r -= x2r;
                    x0i -= x2i;
                    a[j2] = wk2r * x0r - wk2i * x0i;
                    a[j2 + 1] = wk2r * x0i + wk2i * x0r;
                    x0r = x1r - x3i;
                    x0i = x1i + x3r;
                    a[jj1] = wk1r * x0r - wk1i * x0i;
                    a[jj1 + 1] = wk1r * x0i + wk1i * x0r;
                    x0r = x1r + x3i;
                    x0i = x1i - x3r;
                    a[j3] = wk3r * x0r - wk3i * x0i;
                    a[j3 + 1] = wk3r * x0i + wk3i * x0r;
                }
            }
        }
        
        l = m;
    }
    
    if (l < n) {
        
        for (j = 0; j <= l - 2; j += 2) {
            jj1 = j + l;
            x0r = a[j] - a[jj1];
            x0i = a[j + 1] - a[jj1 + 1];
            a[j] += a[jj1];
            a[j + 1] += a[jj1 + 1];
            a[jj1] = x0r;
            a[jj1 + 1] = x0i;
        }
    }
}


void rftsub(int n, t_float *a, int nc, t_float *c)
{
    int j, k, kk, ks;
    t_float wkr, wki, xr, xi, yr, yi;
    
    ks = (nc << 2) / n;
    kk = 0;
    
    for (k = (n >> 1) - 2; k >= 2; k -= 2) {
        j = n - k;
        kk += ks;
        wkr = 0.5 - c[kk];
        wki = c[nc - kk];
        xr = a[k] - a[j];
        xi = a[k + 1] + a[j + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        a[k] -= yr;
        a[k + 1] -= yi;
        a[j] += yr;
        a[j + 1] -= yi;
    }
}


void lpp_makewt(int nw, int *ip, t_float *w)
{
    //    void bitrv2(int n, int *ip, t_float *a);
    int nwh, j;
    t_float delta, x, y;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j <= nwh - 2; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}


void lpp_makect(int nc, int *ip, t_float *c)
{
    int nch, j;
    t_float delta;
    
    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = 0.5;
        c[nch] = 0.5 * cos(delta * nch);
        for (j = 1; j <= nch - 1; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}
