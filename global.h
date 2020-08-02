#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <assert.h>
#define N 16
#define FNUM 1
#define NR 10
#define SEED 115
#define DELTA 0.01
#define dcomplex double complex
#define TEMP 0.1
#define FINITE_TEMP 0
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
int nstates,nstates2,seed,nrz,nspins,nmaj;
dcomplex *H,*chi,*chi1,*chi2,*chi3,*darr;
double *eval;
dcomplex *warray;
double delta;
double beta;
