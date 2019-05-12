#include <stdlib.h>  // malloc, alloc, rand, NULL

#include "utils.h"


double rand01()
// generate random number from 0 to 1
{
   return double(rand() % 10000)/10000;
}
