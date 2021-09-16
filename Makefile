# This will help you compile everything here.
FC       = gfortran
CPPFLAGS = -I/usr/include
LDFLAGS  = -L/usr/lib
FFLAGS   = -g
COMPILE  = $(FC) $(FFLAGS) $(CPPFLAGS) -c
LINK     = $(FC) $(LDFLAGS) -o
OBJECT   = calormod.o monsindex.o

all: monsindex

monsindex: $(OBJECT)
	$(LINK) bin/monsindex.exe $(OBJECT)

%.o: src/%.f90
	$(COMPILE) $<

clean-all: clean
	rm -f bin/*.exe

clean:
	rm -f *.o *.mod *.mod0
