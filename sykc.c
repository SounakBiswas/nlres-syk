#include "global.h"
void diagonalize_wrapper(double *eval, dcomplex * H, int nstates);
void addSykQ(int q, double var);
double make_P(dcomplex *A);
double make_sharpP(dcomplex *A,int eignp);
void makeomat(int nop,dcomplex *opmat);
double calc_osus(int nop,int rz);
double calc_osus2(int nop,int rz);
int combination(int n1,int n2);
void resetH(){
  int i;
  for(i=0; i<nstates*nstates; i++)
    H[i]=0;

}
void resetchi(){
  int i;
  for(i=0; i<nstates*nstates; i++)
    chi[i]=0;

}
void init_genrand64(int);
double overlap(dcomplex *P,dcomplex *O,int q);
double overlap_mc(dcomplex *P,dcomplex *O,int q,int samples);
void makechi(dcomplex *chi,int midx);
void toEigenBasis(dcomplex *chi,dcomplex *temp, int nstates);
dcomplex calc_gE( double eps, double *eval, dcomplex *rotchi );
dcomplex calc_gt( double eps, double *eval, dcomplex *rotchi );
extern void   zcopy(int *, dcomplex *, int*, dcomplex*, int *);
#define TSTEPS 32
#define TOT_TIME 10

void main() {
  int rz;
  char fname1[300],fname2[300], fname3[300];
  char fname4[300];
  FILE *fp;
  int i,j;
  double tgrid[TSTEPS];
  double Tgrid[TSTEPS];
  double dt=TOT_TIME/TSTEPS;
  double dT=TOT_TIME/TSTEPS;
  int nop=combination(N,4)+combination(N,2);




  nmaj=N;
  seed=SEED;
  nrz=NR;
  delta=DELTA;
  init_genrand64(seed);
  nspins=nmaj/2;
  nstates=1<<nspins;
  H=(double complex *)malloc(nstates*nstates*sizeof(double complex));
  eval=(double  *)malloc(nstates*sizeof(double));
  dcomplex *O=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  dcomplex *P=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));


  // dcomplex * evec;
  // evec= (dcomplex *)malloc(totelem*sizeof(dcomplex));
  long long np1;
  int np1mc;
  int neigenp=20; // number of eigenprojectors considered.
  int *eigenps=(int *)malloc(neigenp*sizeof(int));
  for(i=0; i<neigenp; i++){
    eigenps[i]= (int)(1.0*i*nstates/neigenp);
    printf("%d\n",eigenps[i]);
  }
  for(rz=0; rz<nrz; rz++){
    printf("rz=%d\n",rz);
    resetH();
    addSykQ(4,0.66666);
    printf("Hamiltonian created.\n");

    diagonalize_wrapper(eval,H,nstates);
    int idx[4]={1,0,0,0};
    makeKick(chi,idx);
    toEigenbasis(chi,temp,nstates);//rotate chi to diagonal basis

    for(tidx=0; tidx)


      fclose(fp);
    }


  }
  free(eigenps);
  free(H);
  free(eval);
  free(O);
  free(P);


}

