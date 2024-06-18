#include "bashfest.h"
#include "stdlib.h"

void lpp_putsine (t_float *arr, int len);
t_float lpp_boundrand(t_float min, t_float max);


void lpp_putsine (t_float *arr, int len)
{
    int i;
    double twopi;
    twopi = 8.0 * atan2(1.,1.);
    
    for ( i = 0; i < len ; i++) {
        *(arr + i) = sin( twopi * i / len);
    }
}

t_float lpp_boundrand(t_float min, t_float max)
{
    return min + (max-min) * ((t_float)rand()/MY_MAX);
}


void lpp_mycombset(t_float loopt,t_float rvt,int init,t_float *a,t_float srate)
{
    int j;
    
    a[0] =  (3.0 + (loopt * srate + .5));
    a[1] = rvt;
    if(!init) {
        for(j=3; j<(int)*a; j++)
        a[j] = 0;
        a[2] = 3;
    }
}

t_float lpp_mycomb(t_float samp,t_float *a)
{
    t_float temp,*aptr;
    if ( a[2] >= (int) a[0])
        a[2] = 3;
    aptr = a + (int)a[2];
    a[2]++;
    temp = *aptr;
    *aptr = *aptr * a[1] + samp;
    return(temp);
}

void lpp_setweights(t_float *a, int len)
{
    t_float sum = 0.0;
    int i;
    for(i=0;i<len;i++)
    sum += a[i];
    if(sum == 0.0) {
        pd_error(0, "zero odds sum");
    }
    for(i=0;i<len;i++)
    a[i] /= sum;
    for(i=1;i<len;i++)
    a[i] += a[i-1];
}

void  lpp_delset2(t_float *a,int *l,t_float xmax, t_float srate)
{
    /* delay initialization.  a is address of t_float array, l is size-2 int
     * array for bookkeeping variables, xmax, is maximum expected delay */
    
    int i;
    *l = 0;
    *(l+1) = (int)(xmax * srate + .5);
    for(i = 0; i < *(l+1); i++) *(a+i) = 0;
}

void lpp_delput2(t_float x,t_float *a,int *l)
{
    
    /* put value in delay line. See delset. x is t_float */
    
    *(a + (*l)++) = x;
    if(*(l) >= *(l+1)) *l -= *(l+1);
}

t_float lpp_dliget2(t_float *a,t_float wait,int *l,t_float srate)
{
    /* get interpolated value from delay line, wait seconds old */
    register int im1;
    t_float x = wait * srate;
    register int i = x;
    t_float frac = x - i;
    i = *l - i;
    im1 = i - 1;
    if(i <= 0) {
        if(i < 0) i += *(l+1);
        if(i < 0) return(0.);
        if(im1 < 0) im1 += *(l+1);
    }
    return(*(a+i) + frac * (*(a+im1) - *(a+i)));
}

void lpp_butterLopass( t_float *in, t_float *out, t_float cutoff, int frames, int channels, t_float srate)

{
    int channel_to_compute;
    t_float data[8];
    
    for( channel_to_compute = 0; channel_to_compute < channels; channel_to_compute++) {
        lpp_butset( data );
        lpp_lobut(data, cutoff, srate);
        lpp_butter_filter( in, out, data, frames, channels, channel_to_compute);
    }
    
}

void lpp_butterBandpass(t_float *in, t_float *out, t_float center, t_float bandwidth, int frames,int  channels, t_float srate)
{
    int channel_to_compute;
    t_float data[8];
    
    for( channel_to_compute = 0; channel_to_compute < channels; channel_to_compute++) {
        lpp_butset( data );
        lpp_bpbut(data, center, bandwidth, srate);
        lpp_butter_filter( in, out, data, frames, channels, channel_to_compute);
    }
    
}


void lpp_butterHipass(t_float *in, t_float *out, t_float cutoff, int frames,int channels, t_float srate)
{
    int channel_to_compute;
    t_float data[8];
    
    for( channel_to_compute = 0; channel_to_compute < channels; channel_to_compute++) {
        lpp_butset( data );
        lpp_hibut(data, cutoff, srate);
        lpp_butter_filter( in, out, data, frames, channels, channel_to_compute);
    }
    
}

void lpp_butset(t_float *a)
{
    a[6] = a[7] = 0.0;
}

void lpp_lobut(t_float *a, t_float cutoff,t_float SR)
{
    register t_float   c;
    
    c = 1.0 / tan( PI * cutoff / SR);
    a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
    a[2] = a[1] + a[1];
    a[3] = a[1];
    a[4] = 2.0 * ( 1.0 - c*c) * a[1];
    a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    
    
}

void lpp_hibut(t_float *a, t_float cutoff, t_float SR)
{
    
    register t_float  c;
    
    c = tan( PI * cutoff / SR);
    a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
    a[2] = -2.0 * a[1];
    a[3] = a[1];
    a[4] = 2.0 * ( c*c - 1.0) * a[1];
    a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    
}

void lpp_bpbut(t_float *a, t_float formant, t_float bandwidth,t_float  SR)
{
    register t_float  c, d;
    
    c = 1.0 / tan( PI * bandwidth / SR);
    d = 2.0 * cos( 2.0 * PI * formant / SR);
    a[1] = 1.0 / ( 1.0 + c);
    a[2] = 0.0;
    a[3] = -a[1];
    a[4] = - c * d * a[1];
    a[5] = ( c - 1.0) * a[1];
    
}
/* in array can == out array */

void lpp_butter_filter(t_float *in,t_float *out,t_float *a, int frames, int channels, int channel)
{
    
    int i;
    t_float t ,y ;
    
    for( i = channel ; i < frames * channels; i+= channels )
    {
        t = *(in + i) - a[4] * a[6] - a[5] * a[7];
        y = t * a[1] + a[2] * a[6] + a[3] * a[7];
        a[7] = a[6];
        a[6] = t;
        *(out + i) = y;
    }
}

void lpp_rsnset2(t_float cf,t_float bw,t_float scl,t_float xinit,t_float *a,t_float srate)
{
    //  double exp(),cos(),sqrt();
    t_float c,temp;
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

t_float lpp_reson(t_float x,t_float *a)
{
    t_float temp;
    temp = *a * x + *(a+1) * *(a+3) - *(a+2) * *(a+4);
    *(a+4) = *(a+3);
    *(a+3) = temp;
    return(temp);
}

t_float lpp_allpass(t_float samp,t_float *a)
{
    t_float temp,*aptr;
    if ( a[STARTM1] >= (int) a[0]) a[STARTM1] = START;
    aptr = a + (int)a[STARTM1];
    a[STARTM1] ++;
    temp = *aptr;
    *aptr = *aptr * a[1] + samp;
    return(temp - a[1] * *aptr);
}

void lpp_init_reverb_data(t_float *a)
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

void lpp_reverb1me(t_float *in, t_float *out, int inFrames, int out_frames, int nchans,
                   int channel, t_float revtime, t_float dry, t_bashfest *x)
{
    t_float dels[4];// stick into main structure
    t_float **alpo = x->mini_delay ;
    t_float a1,a2,a3,a4;
    int i;
    //  int alsmp ;
    t_float *fltdata = x->reverb_ellipse_data;
    
    int nsects;
    t_float xnorm;
    LSTRUCT *eel = x->eel;
    
    t_float wet;
    //  t_float max;
    t_float srate = x->sr;
    //  t_float max_del = x->max_mini_delay ;
    
    wet = cos(1.570796 * dry);
    dry = sin(1.570796 * dry);
    
    /* combset uses reverb time , mycombset uses feedback */
    for( i = 0; i < 4; i++ ) {
        dels[i] = lpp_boundrand(.005, .1 );
        if(dels[i] < .005 || dels[i] > 0.1) {
            post("reverb1: bad random delay time: %f",dels[i]);
            dels[i] = .05;
        }
        lpp_mycombset(dels[i], revtime, 0, alpo[i], srate);
    }
    
    lpp_ellipset(fltdata,eel,&nsects,&xnorm);
    
    for( i = channel ; i < inFrames * nchans; i += nchans ) {
        
        a1 = lpp_allpass(in[i], alpo[0]);
        a2 = lpp_allpass(in[i], alpo[1]);
        a3 = lpp_allpass(in[i], alpo[2]);
        a4 = lpp_allpass(in[i], alpo[3]);
        
        out[i] = in[i] * dry + lpp_ellipse((a1+a2+a3+a4), eel, nsects,xnorm) * wet;
    }
    
    for( i = channel + inFrames * nchans; i < out_frames * nchans; i += nchans ) {
        
        a1 = lpp_allpass(0.0, alpo[0]);
        a2 = lpp_allpass(0.0, alpo[1]);
        a3 = lpp_allpass(0.0, alpo[2]);
        a4 = lpp_allpass(0.0, alpo[3]);
        
        out[i] =  lpp_ellipse((a1+a2+a3+a4), eel, nsects,xnorm) * wet;
        
    }
    
}

void lpp_feed1(t_float *inbuf, t_float *outbuf, int in_frames, int out_frames,int channels, t_float *functab1,
               t_float *functab2,t_float *functab3,t_float *functab4,int funclen,
               t_float duration, t_float maxDelay, t_bashfest *x)
{
    int i;
    t_float srate = x->sr;
    t_float *delayLine1a = x->mini_delay[0];
    t_float *delayLine2a = x->mini_delay[1];
    t_float *delayLine1b = x->mini_delay[2];
    t_float *delayLine2b = x->mini_delay[3];
    int dv1a[2], dv2a[2];   /* cmix bookkeeping */
    int dv1b[2], dv2b[2];   /* cmix bookkeeping */
    t_float delsamp1a=0, delsamp2a=0 ;
    t_float delsamp1b=0, delsamp2b=0 ;
    t_float delay1, delay2, feedback1, feedback2;
    t_float funcSi, funcPhs;
    t_float putsamp;
    
    /***************************/
    
    funcPhs = 0.;
    
    // read once during note
    
    funcSi = ((t_float) funclen / srate) / duration ;
    
    
    lpp_delset2(delayLine1a, dv1a, maxDelay,srate);
    lpp_delset2(delayLine2a, dv2a, maxDelay,srate);
    
    if( channels == 2 ) {
        lpp_delset2(delayLine1b, dv1b, maxDelay,srate);
        lpp_delset2(delayLine2b, dv2b, maxDelay,srate);
    }
    
    
    for(i = 0; i < out_frames*channels; i += channels ) {
        // buffer loop
        
        delay1 = functab1[ (int) funcPhs ];
        delay2 = functab2[ (int) funcPhs ];
        feedback1 = functab3[ (int) funcPhs ];
        feedback2 = functab4[ (int) funcPhs ];
        
        funcPhs += funcSi;
        if( funcPhs >= (t_float) funclen )
            funcPhs = 0;
        
        putsamp = i < in_frames * channels ? inbuf[i] + delsamp1a*feedback1 : 0.0;
        outbuf[i] = putsamp; // zero instead ??
        
        lpp_delput2( putsamp, delayLine1a, dv1a);
        delsamp1a = lpp_dliget2(delayLine1a, delay1, dv1a,srate);
        
        putsamp = delsamp1a+delsamp2a*feedback2 ;
        
        lpp_delput2( putsamp, delayLine2a, dv2a);
        delsamp2a = lpp_dliget2(delayLine2a, delay2, dv2a, srate);
        outbuf[i] += delsamp2a;
        
        
        if( channels == 2 ) {
            putsamp = i < in_frames * channels ? inbuf[i+1] + delsamp1a*feedback1 : 0.0;
            outbuf[i+1] = putsamp;
            lpp_delput2( putsamp, delayLine1b, dv1b);
            delsamp1b = lpp_dliget2(delayLine1b, delay1, dv1b, srate);
            putsamp = delsamp1b+delsamp2b*feedback2;
            lpp_delput2( putsamp, delayLine2b, dv2b);
            delsamp2b = lpp_dliget2(delayLine2b, delay2, dv2b, srate);
            outbuf[i+1] += delsamp2b;
        }
    }
}

void lpp_setflamfunc1(t_float *arr, int flen)
{
    int i;
    t_float x;
    for ( i = 0; i < flen; i++) {
        x = (t_float)i / (t_float) flen ;
        *(arr + i) = ((x - 1) / (x + 1)) * -1.  ;
    }
}


void lpp_setExpFlamFunc(t_float *arr, int flen, t_float v1,t_float v2,t_float alpha)
{
    int i;
    
    if( alpha == 0 )
        alpha = .00000001 ;
    
    for ( i = 0; i < flen; i++) {
        *(arr + i) = v1 + (v2-v1) * ((1-exp((t_float)i*alpha/((t_float)flen-1.)))/(1-exp(alpha)));
    }
}

void lpp_funcgen1(t_float *outArray, int outlen, t_float duration, t_float outMin, t_float outMax,
                  t_float speed1, t_float speed2, t_float gain1, t_float gain2, t_float *phs1, t_float *phs2,
                  t_float *sine, int sinelen)
{
    t_float si1, si2;
    t_float localSR;
    int i;
    
    localSR = duration * (t_float) outlen ;
    *phs1 *= (t_float) sinelen;
    *phs2 *= (t_float) sinelen;
    si1 = ((t_float)sinelen/localSR)  * speed1;
    si2 = ((t_float)sinelen/localSR)  * speed2;
    
    for( i = 0; i < outlen; i++ ) {
        *(outArray + i) = lpp_oscil(gain1, si1, sine, sinelen, phs1) ;
        *(outArray + i) += lpp_oscil(gain2, si2, sine, sinelen, phs2) ;
    }
    lpp_normtab( outArray, outArray, outMin, outMax, outlen);
}


void lpp_normtab(t_float *inarr,t_float *outarr, t_float min, t_float max, int len)
{
    int i;
    
    t_float imin=9999999999., imax=-9999999999.;
    
    for(i = 0; i < len ; i++) {
        if( imin > inarr[i] )
            imin = inarr[i];
        if( imax < inarr[i] )
            imax = inarr[i];
    }
    for(i = 0; i < len; i++ )
    outarr[i] = lpp_mapp(inarr[i], imin, imax, min, max);
}

t_float lpp_mapp(t_float in,t_float imin,t_float imax,t_float omin,t_float omax)
{
    if( imax == 0.0 )
    {
        return 0.0 ;
    }
    return( omin+((omax-omin)*((in-imin)/(imax-imin))) );
}

t_float lpp_oscil(t_float amp,t_float si,t_float *farray,int len,t_float *phs)
{
    register int i =  *phs;
    *phs += si;
    while(*phs >= len)
        *phs -= len;
    return(*(farray+i) * amp);
}

void lpp_killdc( t_float *inbuf, int in_frames, int channels, t_bashfest *x)
{
    int i,j=1;
    LSTRUCT *eel = x->eel;
    int nsects;
    t_float xnorm;
    t_float *dcflt = x->dcflt;
    
    /* t_float dcflt[64] =
     {3, -1.9999924    , -1.9992482    ,  1.0000000
     ,  .99928019    ,
     -1.9999956    , -1.9964080    ,  1.0000000    ,  .99645999    ,
     -1.9999994    , -1.9805074    ,  1.0000000    ,  .98069401    ,
     .98817413E+00};*/
    
    for( j = 0; j < channels; j++) {
        lpp_ellipset(dcflt,eel,&nsects,&xnorm);
        
        for( i = j; i < in_frames * channels ; i += channels ) {
            inbuf[i] = lpp_ellipse(inbuf[i], eel, nsects,xnorm);
        }
    }
}

void lpp_set_dcflt(t_float *a)
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

void lpp_set_distortion_table(t_float *arr, t_float cut, t_float max, int len)
{
    int i, len2;
    t_float samp;
    
    len2 = len>>1 ;
    for( i = len2; i < len; i++ ) {
        samp = (t_float)(i - len2) / (t_float) len2 ;
        if( samp > cut )
            samp = lpp_mapp( samp, cut, 1.0,  cut, max );
        *(arr + i) = samp;
    }
    for( i = 0; i < len2; i++ )
    *(arr + i) = - *(arr + len - (i+1));
}

t_float lpp_dlookup(t_float samp,t_float *arr,int len)
{
    return arr[(int) (((samp+1.0)/2.0) * (t_float) len)];
}

void lpp_do_compdist(t_float *in,t_float *out,int sampFrames,int nchans,int channel,
                     t_float cutoff,t_float maxmult,int lookupflag,t_float *table,int range,t_float bufMaxamp)
{
    int i;
    t_float rectsamp;
    
    for( i = channel ; i < sampFrames * nchans; i+= nchans ) {
        if( lookupflag) {
            *(out + i) = lpp_dlookup( *(in + i)/bufMaxamp, table, range );
        } else {
            rectsamp = fabs( *(in + i) ) / bufMaxamp;
            if( rectsamp > cutoff ) {
                *(in + i) = *(out + i) *
                lpp_mapp( rectsamp, cutoff, 1.0, cutoff, maxmult);
            }
        }
    }
}

t_float lpp_getmaxamp(t_float *arr, int len)
{
    int i;
    t_float max = 0;
    
    for(i = 0; i < len; i++ ) {
        if( fabs(arr[i]) > max )
            max = fabs(arr[i]);
    }
    return max;
}

void lpp_buildadsr(CMIXADSR *a)
{
    t_float A = a->a;
    t_float D = a->d;
    t_float S = a->s;
    t_float R = a->r;
    t_float f1 = a->v1;
    t_float f2 = a->v2;
    t_float f3 = a->v3;
    t_float f4 = a->v4;
    
    int funclen = a->len;
    t_float *func = a->func;
    t_float total;
    int ipoint = 0;
    int i;
    int segs[4];
    t_float m1,m2;
    total = A + D + S + R ;
    
    segs[0] = (A/total) * funclen;
    segs[1] = (D/total) * funclen;
    segs[2] = (S/total) * funclen;
    segs[3] = funclen - (segs[0]+segs[1]+segs[2]);
    
    if( f1 > 20000. || f1 < -20000. ) {
        f1 = 250.0;
    }
    if( f2 > 20000. || f2 < -20000. ) {
        f2 = 1250.0;
    }
    if( f3 > 20000. || f3 < -20000. ) {
        f3 = 950.0;
    }
    if( f4 > 20000. || f4 < -20000. ) {
        f4 = f1;
    }
    
    if( segs[0] <= 0 || segs[1] <= 0 || segs[2] <= 0 || segs[3] <= 0 ) {
        for( i = 0; i < 4; i++ ) {
            segs[i] = funclen / 4;
        }
    }
    
    for( i = 0 ; i < segs[0]; i++ ) {
        m1 = 1.-(t_float)i/(t_float)(segs[0]);
        m2 = 1. - m1;
        *(func +i ) = f1 * m1 + f2 * m2;
    }
    ipoint = i;
    
    for( i = 0 ; i < segs[1]; i++ ) {
        m1 = 1.-(t_float)i/(t_float)(segs[1]);
        m2 = 1. - m1;
        *(func + i + ipoint) = f2 * m1 + f3 * m2;
    }
    ipoint += i;
    
    for( i = 0 ; i < segs[2]; i++ ) {
        m1 = 1.-(t_float)i/(t_float)(segs[2]);
        m2 = 1. - m1;
        *(func + i + ipoint) = f3;
    }
    ipoint += i;
    
    for( i = 0 ; i < segs[3]; i++ ) {
        m1 = 1.-(t_float)i/(t_float)(segs[3]);
        m2 = 1. - m1;
        *(func + ipoint + i) = f3 * m1 + f4 * m2;
    }
    ipoint += i;
}
