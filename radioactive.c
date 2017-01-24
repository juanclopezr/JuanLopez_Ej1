#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define N 100

int main()
{
  time_t t;
  float par;
  float rands;
  srand((unsigned) time(&t));
  for(int i=0;i<N;i++)
    {
      rands = rand()/(1.0 + RAND_MAX);
      printf("%f\n",rands);
    }
  return 0;
}
