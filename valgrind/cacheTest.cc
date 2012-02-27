#include <stdio.h>
#include <stdlib.h>

/**
Tests cache misses.
*/

int main(int argc, char **argv)
{
  if (argc < 2){
    printf("Usage: cacheTest sizeI sizeJ\n");
    return 1;
  }
  int sI = atoi(argv[1]);
  int sJ = atoi(argv[2]);
  int *arr = new int[sI*sJ]; // double array.

  //for (int i=0; i < sI*sJ; ++i) arr[i] = i;

  for (int i=0; i < sI; ++i)
    for (int j=0; j < sJ; ++j)
       arr[(i * (sJ)) + j ] = i;

  for (int i=0; i < sI; ++i)
    for (int j=0; j < sJ; ++j)
       arr[(j * (sI)) + i ] = i;

  return 1;
}
