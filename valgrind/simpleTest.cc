#include <stdio.h>

/**
Borrows from 
http://cs.ecs.baylor.edu/~donahoo/tools/valgrind/
*/

int main()
{
  char *p;
  int sz;
  p = new char[sz];

  p = new char[12];
  delete(p);

  p = new char[16];
  printf("Size is %d",sz);
  return 0;
}
