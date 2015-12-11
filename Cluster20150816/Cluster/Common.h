#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <string.h>
#include <math.h>
using namespace std;

#define RAND (rand()/(RAND_MAX+1.0))
typedef double (*PEnergy) ( double* ,int );
typedef double (*PForce) (double *, double *, double *, int );
