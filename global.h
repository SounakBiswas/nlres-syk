#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <assert.h>
#define N 18
#define FNUM 1
#define NR 100
#define SEED 115
#define DELTA 0.01
#define dcomplex double complex
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
int nstates,seed,nrz,nspins,nmaj;
dcomplex *H,*chi;
double *eval;
dcomplex *warray;
double delta;
