#include "global.h"
void diagonalize_wrapper(double *eval, dcomplex * H, int nstates);
void addSykQ(int q, double var);
void multiply_quick(dcomplex *A,dcomplex*B,dcomplex*C);
void muldiag_fromL(dcomplex *arr,dcomplex *diag);
void muldiag_fromR(dcomplex *arr,dcomplex *diag);

int combination(int n1,int n2);

double makeKick(dcomplex *chi,int pos);
void toEigenBasis(dcomplex *chi,dcomplex *temp, int nstates);
void construct_heisenberg_op(dcomplex *op0, double t);
void commutator(dcomplex *arr1,dcomplex *arr2, dcomplex *res);
double calc_term1_simp(double t, double T);
double calc_term2_simp(double t, double T);
double zt_calc_term1_simp(double t, double T);
double zt_calc_term2_simp(double t, double T);



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

//construct exp(iHt)
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
  double t,T;
  int tidx,Tidx;
  int nop=combination(N,4)+combination(N,2);




  nmaj=N;
  beta=1/TEMP;
  seed=SEED;
  nrz=NR;
  delta=DELTA;
  init_genrand64(seed);
  nspins=nmaj/2;
  nstates=1<<nspins;
  H=(double complex *)malloc(nstates*nstates*sizeof(double complex));
  eval=(double  *)malloc(nstates*sizeof(double));
  chi1=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  //if(FINITE_TEMP)
  // chi2=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  //else
  // chi2=(dcomplex*)malloc(nstates*sizeof(dcomplex));
  chi2=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  darr=(dcomplex*)malloc(nstates*sizeof(dcomplex));


  // dcomplex * evec;
  // evec= (dcomplex *)malloc(totelem*sizeof(dcomplex));
  long long np1;
  int one=1;
  int np1mc;
  int neigenp=20; // number of eigenprojectors considered.
  nstates2=nstates*nstates;
  double pref[6];
  double htime[6][4];
  chi3=H;
  for(rz=0; rz<nrz; rz++){
    printf("rz=%d\n",rz);
    resetH();
    addSykQ(4,0.66666);
    printf("Hamiltonian created.\n");
    t=20;
    T=20;

    diagonalize_wrapper(eval,H,nstates);
    printf("diagonalised\n");
    //fp=fopen("testHc.dat","w");
    //for(j=0; j<nstates*nstates; j++)
    //  fprintf(fp,"%f  %f\n",creal(H[j]),cimag(H[j]));
    //fclose(fp);
    int idx[1]={0};
    double matelem;
    makeKick(chi1,0); //set chi1=g in computational basis
    toEigenBasis(chi1,chi2,nstates);//rotate chi1=g to eigenbasis,chi2 wksp

    for(tidx=0; tidx<TSTEPS; tidx++){
      for(Tidx=0; Tidx<TSTEPS; Tidx++){
        //printf("tidx %d Tidx %d\n",tidx,Tidx);

        t=tidx*dt;
        T=Tidx*dT;
        if(FINITE_TEMP)
          matelem=calc_term1_simp(t,T);
        else
          matelem=zt_calc_term1_simp(t,T);

        if(FINITE_TEMP)
          matelem=calc_term2_simp(t,T);
        else
          matelem=zt_calc_term2_simp(t,T);

      }
    }
  }
  free(H);
  free(eval);
  free(chi1);
  free(chi2);
  free(darr);
  mkl_free_buffers();
}

