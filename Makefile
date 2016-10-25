FC = mpif90
LD = mpif90
MKLROOT = /home/junteng/intel/mkl
FFLAGS  = -g -O3 -qopt-report=5 -qopt-report-phase=vec -qopt-report-phase=par -align array64byte -assume byterecl -D__INTEL -fpp -xCORE-AVX2 -fpic -ip -check bounds -check uninit -check pointers -traceback -g -fpe0
LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

exec: UTILS.o MBD_rsSCS.o 
	${LD} -o pmbd.x UTILS.o MBD_rsSCS.o ${LDFLAGS}

MBD_rsSCS.o: MBD_rsSCS.f90
	${FC} -o MBD_rsSCS.o -c MBD_rsSCS.f90 ${FFLAGS}

UTILS.o: UTILS.f90
	${FC} -o UTILS.o -c UTILS.f90 ${FFLAGS}

test_C6H6:
	mpiexec -n 4 ./pmbd.x C6H6.in settings.in > out

test_graphene:
	mpiexec -n 4 ./pmbd.x graphene.in settings.in > out

clean: 
	rm -f *.o  *.mod *.x *.optrpt out
