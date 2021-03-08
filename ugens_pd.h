/* THIS IS A TOTALLY HACKED HEADER - NO LONGER ANY GOOD FOR CMIX */

#define MAXSECTS 20
#define RESON_NO_SCL (0.)
#define START 3
#define STARTM1 2 /* for start of comb memory in a array */
#define NCOMBS  6 /* for reverb */
#define NALPASSES 2 /* for reverb */

typedef struct {
    double ps0;
    double ps1;
    double ps2;
    double ps3;
    double c0;
    double c1;
    double c2;
    double c3;
} LSTRUCT ;

typedef struct {
    int len;
    double *func;
    double amp;
    double phs;
    double si;
} CMIXOSC ;

typedef struct {
    double *arr;
    double lpt;
    double rvbt;
    int len;
    int status;
} CMIXCOMB ;

typedef struct {
    double cf;
    double bw;
    double scl;
    double q[5];
} CMIXRESON ;

typedef struct {
    double a;
    double d;
    double s;
    double r;
    double v1;
    double v2;
    double v3;
    double v4;
    double v5;
    double *func;
    int len;
} CMIXADSR ;
