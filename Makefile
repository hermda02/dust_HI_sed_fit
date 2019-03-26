# -*- Makefile -*-

FC      = ifort
OPTIM   = -g -C

FITSDIR = -L/mn/stornext/u3/hke/local/lib -lcfitsio
LAPACK  = -L/mn/stornext/u3/hke/local/lib -llapack -lblas
HEALPIX = -L/mn/stornext/u3/hke/local/lib -lhealpix
HEALINC = -I/mn/stornext/u3/hke/local/include
OUTPUT  = fit_dust_hi

OBJS    = hi_dust_fit.o

fit_dust_hi: $(OBJS)
	$(FC) $(OBJS) $(HEALPIX) $(FITSDIR) -fopenmp -o $(OUTPUT)

# Compilation stage
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) $(LAPACK) $(CFITSIO) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *~ fit_dust_hi
