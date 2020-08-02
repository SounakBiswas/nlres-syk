MKLROOT=/opt/intel/compilers_and_libraries_2018/linux/mkl/
gcc -O0 sykc_simp.c mt19937-64.c utils.c  -o  ./N24SEED1.out -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -g  -m64 -I${MKLROOT}/include
#gcc -O0 sykc_3mat.c mt19937-64.c utils.c  -o  ./N24SEED1.out -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -g  -m64 -I${MKLROOT}/include
#gcc -O0 sykc.c mt19937-64.c utils.c  /usr/local/lib/librefblas.a /usr/local/lib/liblapack.a -o  ./N24SEED1.out   -lm
