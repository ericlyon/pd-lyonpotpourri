#include "PenroseRand.h"

t_float rrand(int *seed)
{
  int i = ((*seed = *seed * 1103515245 + 12345)>>16) & 077777;
  return((t_float)i/16384. - 1.);
}

t_float prand(int *seed)
{
  int i = ((*seed = *seed * 1103515245 + 12345)>>16) & 077777;
  return((t_float)i/32768.);
}
