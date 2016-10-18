FC=mpif90
MKL_ROOT=/home/junteng/intel/mkl/lib/intel64
LAPACKBLAS = -L$(MKL_ROOT) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
FFLAGS =  -O3 -ip -check bounds -check uninit -check pointers -traceback -g -fpe0

install :
	$(FC) $(FFLAGS) -o DFT_MBD_AT_rsSCS.x UTILS.F90 MBD_AT_rsSCS.F90 $(LAPACKBLAS) 

clean: 
	rm -f *.o  *.mod *.x
