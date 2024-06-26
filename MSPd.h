#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <strings.h>
#include <string.h>


#ifndef LPP_VERSION
# define LPP_VERSION "3.0"
#endif
#ifndef LPP_COMPILE_DATE
# define LPP_COMPILE_DATE __DATE__
#endif
#define LPP_ANNOUNCEMENT "-< LyonPotpourri "LPP_VERSION " >-   |  "
#define LPP_VERSION4PD "LyonPotpourri " LPP_VERSION " for Pd"


#define lpp_version(objectname) post("%s: version %s compiled %s",objectname,LPP_VERSION4PDLPP_COMPILE_DATE);

// #define LYONPOTPOURRI_MSG "-<LyonPotpourri 3.0>-"

#define NO_FREE_FUNCTION 0


/* because Max and Pd have different ideas of what A_FLOAT is, use t_floatarg
   to force consistency. Otherwise functions that look good will fail on some
   hardware. Also note that Pd messages cannot accept arguments of type A_LONG. */


#include "m_pd.h"
// #define t_floatarg t_float
#define t_double double

#define atom_getsymarg atom_getsymbolarg


#ifndef PIOVERTWO
#define PIOVERTWO 1.5707963268
#endif
#ifndef TWOPI
#define TWOPI 6.2831853072
#endif
#ifndef PI
#define PI 3.14159265358979
#endif

/*** MSP helper functions **/
void atom_arg_getfloat(t_float *c, long idx, long ac, t_atom *av);
void atom_arg_getsym(t_symbol **c, long idx, long ac, t_atom *av);

// #define potpourri_announce(objname)  post("( %s )\t%s",objname,LYONPOTPOURRI_MSG)

// #define potpourri_announce(objname)  post("%s (  %s  )",LYONPOTPOURRI_MSG,objname)
#define potpourri_announce(objname)  post("%s ( %s )",LPP_ANNOUNCEMENT,objname)
