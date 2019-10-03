#makefile 1
#hiya : hello.o hello_isla.o
#	#gfortran -o hiya hello.o hello_isla.o

#hello.o : hello.f90
#	#gfortran  hello.f90
#hello_isla.o : hello_isla.f90
#	#gfortran hello_isla.f90
#FFLAGS = -ffpe-trap=invalid,zero -g -fcheck=all -fbacktrace

objects = packet.o\
	grid.o\
	search_bisec.o\
	get_dim.o\
	load.o\
	interpolate.o\
	load_spec2.o\
	get_cdf.o\
	op_prop_set.o\
	iarray.o\
	density.o\
	gridset.o\
	mc_sample.o\
	sourceph.o\
	reflect.o\
	tauint2.o\
	stokes.o\
	pl_estimators.o\
	mcpolar.o

#Comp= gfortran

mcgrid : $(objects)
	gfortran -o mcgrid $(objects)

%.o : %.f90
	gfortran -c -g -fcheck=all -ffpe-trap=invalid,zero -fbacktrace $<
