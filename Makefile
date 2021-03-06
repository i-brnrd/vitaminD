#makefile 1
#hiya : hello.o hello_isla.o
#	#gfortran -o hiya hello.o hello_isla.o

#hello.o : hello.f90
#	#gfortran  hello.f90
#hello_isla.o : hello_isla.f90
#	#gfortran hello_isla.f90
#FFLAGS = -ffpe-trap=invalid,zero -g -fcheck=all -fbacktrace

objects = constants.o\
	optical_properties.o\
	packet.o\
	grid.o\
	search_bisec.o\
	get_dim.o\
	load.o\
	interpolate.o\
	load_spec2.o\
	get_cdf.o\
	optical_properties_init.o\
	jacques_verification.o\
	iarray.o\
	density.o\
	gridset.o\
	mc_sample.o\
	sourceph.o\
	n_interface.o\
	tauint3.o\
	stokes.o\
	pl_estimators.o\
	mcpolar.o

#Comp= gfortran

mcgrid : $(objects)
	gfortran -o mcgrid $(objects)

%.o : %.f90
	gfortran -c -g -fcheck=all -ffpe-trap=invalid,zero -fbacktrace $<
