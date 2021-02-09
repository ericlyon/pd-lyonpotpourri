#include "fftease.h"
int power_of_two(int test)
{
  int limit = MAX_N;
  int compare = 1;
  //  post("testing what we thing is an int:%d",test);
  do {
    if(test == compare) {
      //      post("good power of 2 found!");
      return 1;
    }
    compare *= 2;
  } while (compare <= limit);

  return 0;
}

