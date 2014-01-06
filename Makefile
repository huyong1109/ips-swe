defaults: swe

CFLAGS	         = -O3 -Wall
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

swe: main.o ips_swe.o  chkopts
	-${CLINKER} -o swe  main.o ips_swe.o  ${PETSC_LIB}
run: swe
	-${MPIEXEC} -np 1 ./swe  -px 1 -py 1 -meshsize 50 -tsize 0.1 -tfinal 15 -tsmax 10  -log_summary 
#	-${MPIEXEC} -np 24 ./swe.13d -preload 0 -meshsize 64 -tsize 2.5e-3 -tfinal 15 -tsmax 6000 -tsback 400 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 54 ./swe.13d -preload 0 -meshsize 64 -tsize 2.5e-3 -tfinal 15 -tsmax 6000 -tsback 400 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 54 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 24 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 6 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary

