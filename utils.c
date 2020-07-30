#include"global.h"
extern void zcopy(int *nt,dcomplex *a,int *ix,dcomplex *b, int *iy);
extern void zdotc(dcomplex *res,int *nx,dcomplex *x,int *incx, dcomplex *y, int *incy);
extern void zscal(int *nt,dcomplex *a,dcomplex *b, int *iy);
extern void zheev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda,
    double* w, dcomplex* work, int* lwork, double* rwork, int* info );
extern void  zgemm(char *, char *, int *, int *, int *, dcomplex *, dcomplex *, int *, dcomplex *, int *, dcomplex *, dcomplex *, int *);
double genrand64_real2();
double random_normal(double mu,double sigma){
  static int rcall=1;
  static double rand1,rand2;
  rcall=!rcall;

  if(rcall)
    return mu+sigma*rand2;
  else{
    double u1,u2;
    u1=genrand64_real2();
    u2=genrand64_real2();
    rand1=sqrt(-2.0*logl(u1))*cos(2*M_PI*u2);
    rand2=sqrt(-2.0*logl(u1))*sin(2*M_PI*u2);
    return mu+sigma*rand1;

  }

}
int compare (const void * a, const void * b){
  return ( *(int*)a - *(int*)b );
}
//double random_normal(double mu,double sigma){
//  return 1;
//}
int factorial(int num){
  int res=1;
  while(num>1)
    res*=(num--);
  return res;
}

unsigned long long combination(int num1,int num2){
  unsigned long long res=1;
  int k=1;
  if(num2>num1-num2)
    num2=num1-num2;
  while(k<=num2){
    res=res*num1/(k);
    num1--;
    k++;
  }
  //printf("%lld \n",res);
  return res;
}

void addSykQ(int q,double Jq){
  int *idx;
  int i,j,p;
  int count=0;
  int mparity,midx,sidx,state,tstate;
  double complex coeff;
  double sigma,coupling;
  idx=(int *)malloc((q+1)*sizeof(int));
  for(i=0;i<q;i++)
    idx[i]=i;
  sigma=Jq*sqrt(factorial(q-1)/pow(N,q-1));
  p=0;
  idx[q]=N;
  double complex mul=1;
  double pref=pow(0.5,q/2);
  if((q/2)%2)
    mul*= I;
  FILE *fp;
  FILE *gp;

  while(idx[q]==N){ // loop over combinations of NCq sets of majoranas
    count++;


    coupling= random_normal(0.0,sigma);
    //coupling=1.0;
    for (state=0; state<nstates; state++){
      coeff=1 +0*I;
      //get sigma_z products
      for(i=0;i<q; i+=2){
        for(j=idx[i]/2; j< idx[i+1]/2; j++) {
          if(state & (1<<j)){
            coeff=coeff*(-1+0*I);
          }
        }
      }

      tstate=state;
      //compute effect of spinflip terms.
      for(i=0; i<q; i++){
        midx= idx[i];
        sidx= idx[i]/2;
        mparity= idx[i]%2;
        if(idx[i]%2) {
          coeff=coeff*(0+1*I);
          if(tstate &(1<<sidx))
            coeff=coeff*(-1+0*I);
        }
        tstate=tstate^(1<<sidx);

      }
      //      printf("%d %d %f %f\n",state,tstate,creal(coupling*coeff),cimag(coupling*coeff));
      H[state+tstate*nstates] += mul*pref*coupling*coeff;


    }
    idx[0]++;
    while(idx[p]==idx[p+1]){
      idx[p]=p;
      idx[++p]++;
      if(p==q)
        break;

    }
    p=0;
  }
  free(idx);
}

double makeKick(dcomplex *chi,int *idx){
  int i,j,p;
  int count=0;
  int mparity,midx,sidx,state,tstate;
  double complex coeff;
  double sigma,coupling;
  p=0;
  double complex mul=1;
  double pref=pow(0.5,q/2);
  pref=pow(0.5,q);
  pref=1/sqrt(nstates);
  //pref=1.0;
  int nstates2=nstates*nstates;
  int one=1;
  if((q/2)%2)
    mul*= I;
  for(i=0;i<nstates;i++)
    O[i]=0;
  char fname[300];

    for(i=0;i<nstates*nstates;i++)
      chi[i]=0;
    for (state=0; state<nstates; state++){
      coeff=1 +0*I;
      for(i=0;i<q; i+=2){
        for(j=idx[i]/2; j< idx[i+1]/2; j++) {
          if(state & (1<<j)){
            coeff=coeff*(-1+0*I);
          }
        }
      }

      tstate=state;
      //compute effect of spinflip terms.
      for(i=0; i<q; i++){
        midx= idx[i];
        sidx= idx[i]/2;
        mparity= idx[i]%2;
        if(idx[i]%2) {
          coeff=coeff*(0+1*I);
          if(tstate &(1<<sidx))
            coeff=coeff*(-1+0*I);
        }
        tstate=tstate^(1<<sidx);

      }
      chi[state+tstate*nstates] += mul*pref*coeff;

    }
}
void commutator(dcomplex *arr1,dcomplex *arr2, double *res){
  char trans1[1],trans2[1];
  trans1[0]='N';
  trans2[0]='N';
  dcomplex alpha=1;
  dcomplex beta=0;
  zgemm(trans1, trans2, &nstates, &nstates, &nstates, &alpha, arr1, &nstates, arr2, &nstates, &beta, res, &nstates);
  alpha=-1;
  beta=1;
  zgemm(trans1, trans2, &nstates, &nstates, &nstates, &alpha, arr2, &nstates, arr1, &nstates, &beta, res, &nstates);
}

//in-place construction o heisenberg operator
void construct_heisenberg_op(dcomplex *op0, dcomplex *opt,double t){
  int i;
  for(i=0; i<nstates;i++){
    larr[i]=cexp(I*t*eval[i])
  }
  muldiag_fromL(op0,larr);
  for(i=0; i<nstates;i++){
    larr[i]=cexp(-I*t*eval[i])
  }
  muldiag_fromR(op0,larr);
}
void muldiag_fromL(dcomplex *arr,double *diag){
  int i;
  //zcopy()
  int inc=1;
  dcomplex fac;
  for(i=0;i<nstatess; i++){
    fac=diag[i];
    zscal(&nstates,&fac,arr+i,&nstates);
  }
}
void muldiag_fromR(dcomplex *arr,double *diag){
  int i;
  //zcopy()
  int inc=1;
  int one=1;
  for(i=0;i<nstates; i++){
    fac=diag[i];
    zscal(&nstates,&fac,arr+i*n1,&one);
  }
}


void diagonalize_wrapper(double *eval, dcomplex * H, int nstates){
  char jobz[1],uplo[1];
  jobz[0]='V';
  uplo[0]='L';
  int info=0;
  int lwork=-1;
  int i;
  dcomplex *work;
  dcomplex wopt[1];

  double *rwork= (double *) malloc((3*nstates-2)*sizeof(double));
  zheev(jobz,uplo,&nstates, H, &nstates, eval, wopt, &lwork, rwork, &info);
  lwork=(int)creal(wopt[0]);
  work=(dcomplex *)malloc(lwork*sizeof(dcomplex));
  zheev(jobz,uplo,&nstates, H, &nstates, eval, work, &lwork, rwork, &info);
  free(work);
  free(rwork);
}
//rotates to eigenbasis, needs temp;
void toEigenBasis(dcomplex *chi,dcomplex *temp, int nstates){
  char trans1[1];
  char trans2[1];
  temp=(dcomplex *)malloc(nstates*nstates*sizeof(dcomplex));
  trans1[0]='N';
  trans2[0]='C';
  double epszero=1e-8;
  dcomplex alpha=1.0;
  dcomplex beta=0.0;
  int i;
  for(i=0;i<nstates;i++){
   assert(creal(chi[nstates*i]-chi[i])<epszero);
   assert(cimag(chi[nstates*i]-chi[i])<epszero);
  }
//  printf("this is okay \n");

  zgemm(trans1, trans1, &nstates, &nstates, &nstates, &alpha, chi, &nstates, H, &nstates, &beta, temp, &nstates);
  zgemm(trans2, trans1, &nstates, &nstates, &nstates, &alpha, H, &nstates, temp, &nstates, &beta, chi, &nstates);
  free(temp);
  for(i=0;i<nstates;i++){
   assert(fabs(creal(chi[nstates*i]-chi[i]))<epszero);
   assert(fabs(cimag(chi[nstates*i]+chi[i]))<epszero);
  }

}
