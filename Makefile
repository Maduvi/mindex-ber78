# This will help you compile everything here.
FC       = gfortran
CPPFLAGS = -I/usr/include
LDFLAGS  = -L/usr/lib
FFLAGS   = -g
COMPILE  = $(FC) $(FFLAGS) $(CPPFLAGS) -c
LINK     = $(FC) $(LDFLAGS) -o
OBJECT1   = calormod.o single.o
OBJECT2   = calormod.o transient.o

all: single transient

transient: $(OBJECT2)
	$(LINK) bin/transient.exe $(OBJECT2)

single: $(OBJECT1)
	$(LINK) bin/single.exe $(OBJECT1)

%.o: src/%.f90
	$(COMPILE) $<

clean-all: clean
	rm -f bin/*.exe

clean:
	rm -f *.o *.mod *.mod0
