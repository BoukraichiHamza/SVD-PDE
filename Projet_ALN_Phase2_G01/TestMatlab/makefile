include Make.inc

# Machines salle TP
OBJS  = main.o 
MODS  = m_subspace_iter.o

all: lib main

main:  $(MODS) $(OBJS)
	$(FC) -o $@ $^  $(LIBS)

lib:
	(cd matgen; make)

matlabsetup:
	/applications/matlab/bin/mex -setup
	/applications/matlab/bin/mex -setup FORTRAN
	cp mex_FORTRAN_glnxa64.xml mex_C_glnxa64.xml $(HOME)/.matlab/R2016a/

mexfile:
	/applications/matlab/bin/mex -Dmex -largeArrayDims -c m_subspace_iter.F90
	/applications/matlab/bin/mex -Dmex -largeArrayDims mex_subspace_iter_sv.c m_subspace_iter.o $(PLIBS) -lgfortran

main.o: matgen/libaux.a

clean:
	(rm -f *.o *.mod main; cd matgen; make clean)

%.o: %.f90
	$(FC) -c $<  

%.o: %.F90
	$(FC) -c $<  

%.o: %.f
	$(FC) -c $<  

