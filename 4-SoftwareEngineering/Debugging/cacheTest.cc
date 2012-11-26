#include <stdio.h>
#include <stdlib.h>

/**
Tests cache misses.
*/

int main(int argc, char **argv)
{
  if (argc < 4){
    printf("Usage: cacheTest fast sizeI sizeJ\nIn first input, use 1 for fast code (cache-smart) and anything else for slow (cache-poor) execution.\n");
    return 1;
  }
  int runFast = atoi(argv[1]);
  long sI = atoi(argv[2]);
  long sJ = atoi(argv[3]);

  printf("Operating on matrix of size %d by %d\n",sI,sJ);

  long *arr = new long[sI*sJ]; // double array.

  // El mejor!
  //for (int i=0; i < sI*sJ; ++i) arr[i] = i;

  if (runFast == 1){

    for (long i=0; i < sI; ++i)
      for (long j=0; j < sJ; ++j)
         arr[(i * (sJ)) + j ] = i;

  }else{

    for (long i=0; i < sI; ++i)
      for (long j=0; j < sJ; ++j)
         arr[(j * (sI)) + i ] = i;

  }

  return 1;
}
