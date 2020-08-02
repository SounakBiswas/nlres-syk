#include "global.h"
void diagonalize_wrapper(double *eval, dcomplex * H, int nstates);
void addSykQ(int q, double var);
void multiply_quick(dcomplex *A,dcomplex*B,dcomplex*C);
void muldiag_fromL(dcomplex *arr,dcomplex *diag);
void muldiag_fromR(dcomplex *arr,dcomplex *diag);
void mkl_free_buffers();

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
#define TSTEPS 64
#define TOT_TIME 20

//construct exp(iHt)
void main() {
  int rz;
  FILE *fp;
  char fname1[300],fname2[300], fname3[300];
  char fname4[300];
  char fname[300];
  int i,j;
  double tgrid[TSTEPS];
  double Tgrid[TSTEPS];
  double dt=(TOT_TIME/(1.0*TSTEPS));
  double dT=(TOT_TIME/(1.0*TSTEPS));
  double t,T;
  double *data1,*data2;
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
  data1=(double*)malloc(TSTEPS*TSTEPS*sizeof(double));
  data2=(double*)malloc(TSTEPS*TSTEPS*sizeof(double));
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

        t=tidx*dt;
        T=Tidx*dT;
        if(FINITE_TEMP)
          matelem=calc_term1_simp(t,T);
        else
          matelem=zt_calc_term1_simp(t,T);
        data1[tidx+Tidx*TSTEPS]=matelem;

        if(FINITE_TEMP)
          matelem=calc_term2_simp(t,T);
        else
          matelem=zt_calc_term2_simp(t,T);
        data2[tidx+Tidx*TSTEPS]=matelem;

      }
    }
    sprintf(fname,"./outfiles/pump_probe_N%d_seed%d_rz%d.dat",N,FNUM,rz);
    fp=fopen(fname,"w");
    for(i=0; i<TSTEPS*TSTEPS; i++){
      fprintf(fp,"%.16f \n",data1[i]);
    }
    fclose(fp);

    sprintf(fname,"./outfiles/rephasing_N%d_seed%d_rz%d.dat",N,FNUM,rz);
    fp=fopen(fname,"w");
    for(i=0; i<TSTEPS*TSTEPS; i++){
      fprintf(fp,"%.16f \n",data2[i]);
    }
    fclose(fp);
  }
  free(H);
  free(eval);
  free(chi1);
  free(chi2);
  free(darr);
  mkl_free_buffers();
}

