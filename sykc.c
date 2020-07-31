#include "global.h"
void diagonalize_wrapper(double *eval, dcomplex * H, int nstates);
void addSykQ(int q, double var);

int combination(int n1,int n2);

double makeKick(dcomplex *chi,int pos);
void toEigenBasis(dcomplex *chi,dcomplex *temp, int nstates);
void construct_heisenberg_op(dcomplex *op0, double t);
void commutator(dcomplex *arr1,dcomplex *arr2, dcomplex *res);


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
  double t,T;
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
  dcomplex *chi1=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  dcomplex *chi2=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));
  dcomplex *temp=(dcomplex*)malloc(nstates*nstates*sizeof(dcomplex));


  // dcomplex * evec;
  // evec= (dcomplex *)malloc(totelem*sizeof(dcomplex));
  long long np1;
  int one=1;
  int np1mc;
  int neigenp=20; // number of eigenprojectors considered.
  int *eigenps=(int *)malloc(neigenp*sizeof(int));
  int nstates2=nstates*nstates;
  for(i=0; i<neigenp; i++){
    eigenps[i]= (int)(1.0*i*nstates/neigenp);
    printf("%d\n",eigenps[i]);
  }
  for(rz=0; rz<nrz; rz++){
    printf("rz=%d\n",rz);
    resetH();
    addSykQ(4,0.66666);
    printf("Hamiltonian created.\n");
    t=20;
    T=20;

    diagonalize_wrapper(eval,H,nstates);
    printf("diagonalised\n");
    int idx[1]={0};
    makeKick(chi1,0); //set chi1=g in computational basis
    toEigenBasis(chi1,chi2,nstates);//rotate chi1=g to eigenbasis,chi2 wksp
    zcopy(&nstates2,chi1,&one,chi2,&one); // set chi2=g,eigenbasis
    construct_heisenberg_op(chi1,t+T); //chi1= g_{T+t}
    construct_heisenberg_op(chi2,t); //set chi2=g_{t}
    commutator(chi1,chi2,temp); //temp= [g_{T+t},g_t]
    commutator(temp,chi2,chi1); // chi1= [temp,g_t]
    makeKick(chi2,0); //set chi2=g in computational basis
    toEigenBasis(chi2,temp,nstates);//rotate chi2=g to eigenbasis,temp wksp
    commutator(chi1,chi2,temp); // chi1= 
    printf("comm=%f\n",creal(temp[0]));

    makeKick(chi1,0); //set chi1=g in computational basis
    toEigenBasis(chi1,chi2,nstates);//rotate chi1=g to eigenbasis,chi2 wksp
    zcopy(&nstates2,chi1,&one,chi2,&one); // set chi2=g,eigenbasis
    construct_heisenberg_op(chi1,t+T); //chi1= g_{T+t}
    construct_heisenberg_op(chi2,t); //set chi2=g_{t}
    commutator(chi1,chi2,temp); //temp= [g_{T+t},g_t]
    makeKick(chi2,0); //set chi2=g in computational basis
    toEigenBasis(chi2,temp,nstates);//rotate chi2=g to eigenbasis,temp wksp
    commutator(temp,chi2,chi1); // chi1= [temp,g_t]
    commutator(chi1,chi2,temp); // chi1= 
    printf("comm=%f\n",creal(temp[0]));


  }
  free(H);
  free(eval);
  free(chi1);
  free(chi2);
  free(temp);


}

