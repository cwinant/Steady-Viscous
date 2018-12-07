$(info have you set "export OMP_NUM_THREADS=XXX"   ?)
FC = gfortran
PROGRAMS = SV
FCFLAGS = -O3 -ffast-math -march=corei7 -fopenmp
#FCFLAGS +=  -g -fbounds-check ##for debug
LDFLAGS = -lm -fopenmp -L/usr/local/lib /home/abraca/Software/SuperLU_MT_3.1/lib/libsuperlu_mt_OPENMP.a /home/abraca/Software/SuperLU_MT_3.1/lib/libblas_OPENMP.a -llapack -lblas
# flags forall (e.g. look for system .mod files, required in gfortran)
# "make" builds all
all: $(PROGRAMS)
# Using Fortran MODULES:

SV.o: mod_SparseHB.o mod_BuildMatrices.o
SV: mod_SparseHB.o mod_BuildMatrices.o

# Linking
%: %.o
	$(FC) $(FCFLAGS) -o  $@  $^  /home/abraca/Software/SuperLU_MT_3.1/EXAMPLE/c_bridge_pdgssv.o $(LDFLAGS)

# Compiling (fortran)
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<


# Utility targets

clean:
	rm -f *.o *.mod 
