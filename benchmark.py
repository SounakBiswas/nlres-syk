import numpy as np
import scipy.linalg as LA
import numpy.linalg as nLA
import math
import itertools
import time;
import sys
N=8 #N of syk model
J2=1;
J4=1;
nstates=2**(N//2)
H=np.zeros([nstates,nstates],dtype='complex')
majtime=0;
#print("enter seed")
#seedr=int(input())
seedr=11
np.random.seed(seedr)
print(nstates)
coup=np.loadtxt("couplings.dat")
ind=np.loadtxt("idxs.dat",dtype=int)



def addSYK(q,Jq) :
    global H
    #var=Jq*(math.factorial(q-1)*2**(q-1)/(q*N**(q-1)))
    var=1.0*Jq*np.sqrt(1.0*math.factorial(q-1)/N**(q-1))
    mul=1j**(q/2)   #imposes hermiticity for odd q/2
    t=time.time()
    interactions=list(itertools.combinations(range(N),q))
    print(len(interactions))
    nint=len(interactions)
    print (time.time()-t, "interactions evaluated",len(interactions))
    st=0
    pref=0.5**(q/2.0)
    #for qtuple in interactions : 
    for qidx in range(nint) :
      qtuple=ind[qidx,:]
      #print("qtup",st)
      st+=1
      #coupling=np.random.normal(0,var);
      coupling=coup[st-1]
      #coupling=1.0
      for state in range(nstates) :
          coeff=1;
          #accumulate -1*sigma^{3} for strings
          for i in range(0,q,2) :
              for j in range(qtuple[i]//2, qtuple[i+1]//2) :
                  if(state&(2**j)):
                      coeff *= -1;
          tstate=state;
          temp=0;
          for i in range(0,q) :
              midx=qtuple[i];
              sidx=qtuple[i]//2;
              mparity=qtuple[i]%2;
              if(mparity==1) :
                  coeff*=1j;
                  if( tstate&(2**sidx)>0) :
                     coeff*=-1;
              tstate=tstate^(2**sidx)
              

          H[state,tstate]+=pref*coupling*coeff;

def makeG(pos) :
    M=np.zeros([nstates,nstates])
    for state in range(nstates) :
        coeff=1;
        #accumulate -1*sigma^{3} for strings
        for j in range(0, pos//2) :
                if(state&(2**j)):
                    coeff *= -1;
        tstate=state;
        temp=0;
        midx=pos;
        sidx=pos//2;
        mparity=pos%2;
        if(mparity==1) :
            coeff*=1j;
            if( tstate&(2**sidx)>0) :
               coeff*=-1;
        tstate=tstate^(2**sidx)
            

        M[tstate,state]+= coeff;
        #print tstate,st,"state=",state,"Hs=",gstate[state],coeff*gstate[state]
        #sys.stdin.read(1)
    return M

addSYK(4,1)
e,v=nLA.eigh(H)
#v=np.loadtxt('testHc.dat',usecols=0)+1j*np.loadtxt('testHc.dat',usecols=1)
#v=np.reshape(v,[nstates,nstates],order='F');
#Hd=np.dot(np.transpose(np.conjugate(v)),np.dot(H,v))
#print "norm"
#print LA.norm(np.diagonal(Hd)-e)
#print Hd[0,0]
#G=makeG(0)
#print np.shape(G)
#temp=np.reshape(G,[nstates*nstates],order='F')
#np.savetxt('testGpyo.dat',temp)
#G=np.dot(np.transpose(np.conjugate(v)),np.dot(G,v))
#temp=np.reshape(G,[nstates*nstates],order='F')
#temp2=np.loadtxt('testGc.dat',usecols=0)+1j*np.loadtxt('testGc.dat',usecols=2)
#print LA.norm(temp-temp2)
#np.savetxt("testGpy.dat",temp)
#print np.diagonal(G)
T=20
t=20
G=makeG(0)
G=np.dot(np.transpose(np.conjugate(v)),np.dot(G,v))
diag1=np.diag(np.exp(1j*(T+t)*e))
diag1c=np.conjugate(diag1)
GTpt=np.dot(diag1,np.dot(G,diag1c))

diag1=np.diag(np.exp(1j*(t)*e))
diag1c=np.conjugate(diag1)
Gt=np.dot(diag1,np.dot(G,diag1c))

temp=np.dot(GTpt,Gt)-np.dot(Gt,GTpt)
temp=np.dot(temp,Gt)-np.dot(Gt,temp)
comm=np.dot(temp,G) - np.dot(G,temp)
print comm[0,0]

temp=np.dot(GTpt,Gt)-np.dot(Gt,GTpt)
temp=np.dot(temp,G)-np.dot(G,temp)
comm=np.dot(temp,G) - np.dot(G,temp)
print comm[0,0]

print e[0:5]
sys.stdin.read(1)
