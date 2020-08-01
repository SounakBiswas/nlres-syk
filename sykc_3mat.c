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
#define TSTEPS 10
#define TOT_TIME 10

//construct exp(iHt)
void init_darr(dcomplex *darr,double t){
  int i;
  for(i=0; i<nstates; i++)
    darr[i]=cexp(1j*eval[i]*t);
}
void init_prefs(double t, double T, double (*h)[4], double *pref){
  pref[0]=1;
  h[0][0]=0;
  h[0][1]=0;
  h[0][2]=t;
  h[0][3]=t+T;

  pref[1] =-1;
  h[1][0]=T+t;
  h[1][1]=t;
  h[1][2]=0;
  h[1][3]=0;

  pref[2]=-1;
  h[2][0]=0;
  h[2][1]=0;
  h[2][2]=T+t;
  h[2][3]=t;

  pref[3]=-2;
  h[3][0]=0;
  h[3][1]=t;
  h[3][2]=T+t;
  h[3][3]=0;

  pref[4]=1;
  h[4][0]=t;
  h[4][1]=T+t;
  h[4][2]=0;
  h[4][3]=0;

  pref[5]=2;
  h[5][0]=0;
  h[5][1]=T+t;
  h[5][2]=t;
  h[5][3]=0;

}

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
  dcomplex *darr=(dcomplex*)malloc(nstates*sizeof(dcomplex));


  // dcomplex * evec;
  // evec= (dcomplex *)malloc(totelem*sizeof(dcomplex));
  long long np1;
  int one=1;
  int np1mc;
  int neigenp=20; // number of eigenprojectors considered.
  int *eigenps=(int *)malloc(neigenp*sizeof(int));
  int nstates2=nstates*nstates;
  double pref[6];
  double htime[6][4];
  dcomplex *chi3;
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

    for(tidx=0; tidx<TSTEPS; tidx++){
      for(Tidx=0; Tidx<TSTEPS; Tidx++){
        //printf("tidx %d Tidx %d\n",tidx,Tidx);
        dcomplex matelem=0;

        t=tidx*dt;
        T=Tidx*dT;
        t=20;
        T=20;
        init_prefs(t,T,htime,pref);
        makeKick(chi1,0); //set chi1=g in computational basis
        printf("chi1 init %f\n",creal(chi1[0]));
        toEigenBasis(chi1,chi2,nstates);//rotate chi1=g to eigenbasis,chi2 wksp
        printf("chi1 init %f\n",creal(chi1[0]));
        zcopy(&nstates2,chi1,&one,chi2,&one); // set chi2=g,eigenbasis
        printf("chi2 init %f\n",creal(chi2[0]));


        int term;
        double tol_eps=1e-5;
        double targ;
        
         
        for(term=0; term<6; term++){
          zcopy(&nstates2,chi1,&one,chi2,&one); // set chi2=g,eigenbasis
          targ=-htime[term][0];

          if(fabs(targ)>tol_eps){
            init_darr(darr,targ);
            muldiag_fromR(chi2,darr);
          }
          targ=htime[term][0]-htime[term][1];

          if(fabs(targ)>tol_eps){
            init_darr(darr,targ);
            muldiag_fromL(chi2,darr);
          }
          multiply_quick(chi1,chi2,chi3);

          targ=htime[term][1]-htime[term][2];
          if(fabs(targ)>tol_eps){
            init_darr(darr,targ);
            muldiag_fromL(chi3,darr);
          }

          multiply_quick(chi1,chi3,chi2);

          targ=htime[term][2]-htime[term][3];
          if(fabs(targ)>tol_eps){
            init_darr(darr,targ);
            muldiag_fromL(chi2,darr);
          }
          multiply_quick(chi1,chi2,chi3);


          targ=htime[term][3];
          if(fabs(targ)>tol_eps){
            init_darr(darr,targ);
            muldiag_fromL(chi3,darr);
          }
          matelem+=chi3[0]*pref[term];
          printf("term=%d, mat=%f +i%f\n",term,creal(chi3[0]),cimag(chi3[0]));

        }
        printf("t=%f T=%f, whatve we got:%f +I%f\n",t,T,creal(matelem),cimag(matelem));

      }
    }
  }
  free(H);
  free(eval);
  free(chi1);
  free(chi2);
  free(temp);
}

