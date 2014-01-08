defaults: swe

CFLAGS	         = -O3 -Wall
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

swe: main.o ips_swe.o misc.o  chkopts
	-${CLINKER} -o swe  main.o ips_swe.o misc.o  ${PETSC_LIB}
run: swe
	-${MPIEXEC} -np 2 ./swe  -px 2 -py 1 -meshsize 10 -tsize 0.1 -tfinal 15 -tsmax 3  -log_summary -vec_view_matlab -draw_pause 1
#	-${MPIEXEC} -np 24 ./swe.13d -preload 0 -meshsize 64 -tsize 2.5e-3 -tfinal 15 -tsmax 6000 -tsback 400 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 54 ./swe.13d -preload 0 -meshsize 64 -tsize 2.5e-3 -tfinal 15 -tsmax 6000 -tsback 400 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 54 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 24 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary
#	-${MPIEXEC} -np 6 ./swe.13d -preload 0 -meshsize 1024 -tsize 1.6e-4 -tfinal 10 -tsmax 10 -tscomp 1 -log_summary

