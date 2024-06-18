/* THIS IS A TOTALLY HACKED HEADER - NO LONGER ANY GOOD FOR CMIX */

#define MAXSECTS 20
#define RESON_NO_SCL (0.)
#define START 3
#define STARTM1 2 /* for start of comb memory in a array */
#define NCOMBS  6 /* for reverb */
#define NALPASSES 2 /* for reverb */

typedef struct {
  t_float ps0;
  t_float ps1;
  t_float ps2;
  t_float ps3;
  t_float c0;
  t_float c1;
  t_float c2;
  t_float c3;
} LSTRUCT ;

typedef struct {
  int len;
  t_float *func;
  t_float amp;
  t_float phs;
  t_float si;
} CMIXOSC ;

typedef struct {
  t_float *arr;
  t_float lpt;
  t_float rvbt;
  int len;
  int status;
} CMIXCOMB ;

typedef struct {
  t_float cf;
  t_float bw;
  t_float scl;
  t_float q[5];
} CMIXRESON ;

typedef struct {
  t_float a;
  t_float d;
  t_float s;
  t_float r;
  t_float v1;
  t_float v2;
  t_float v3;
  t_float v4;
  t_float v5;
  t_float *func;
  int len;
} CMIXADSR ;
