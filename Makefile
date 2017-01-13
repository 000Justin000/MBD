include Makefile.in

exec: UTILS.o MBD_rsSCS.o 
	${LD} -o pmbd.x UTILS.o MBD_rsSCS.o ${LDFLAGS}

MBD_rsSCS.o: MBD_rsSCS.f90
	${FC} -o MBD_rsSCS.o -c MBD_rsSCS.f90 ${FFLAGS}

UTILS.o: UTILS.f90
	${FC} -o UTILS.o -c UTILS.f90 ${FFLAGS}

test_C6H6:
	mpiexec -n 4 ./pmbd.x examples/C6H6.in settings.in > out

test_graphene:
	mpiexec -n 4 ./pmbd.x examples/graphene_111.in settings.in > out

test_PBBA:
	mpiexec -n 4 ./pmbd.x examples/PBBA.in settings.in > out

test_TiOCHF:
	mpiexec -n 4 ./pmbd.x examples/TiOCHF.in settings.in > out

clean: 
	rm -f *.o  *.mod *.x *.optrpt out
