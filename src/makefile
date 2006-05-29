TOPDIR=..
SRCDIR=$(TOPDIR)/src
INSTALLDIR=$(HOME)
INSTALLLIBDIR=$(HOME)/lib
INSTALLINCDIR=$(HOME)/include
FC=ifort
FCARGS=-132 -O -check bounds
##INC=-I/home/users/raynaud/lib/lapack95_modules
LIB=-L/usr/local/install/lapack-3.0/lib/ -llapack -lblas -llapack -lblas
TARGET=spanlib


######################################################################

.SUFFIXES: .f90 .a .o .mod 

%.a: %.o
	ar -r $@ $<
##	ld -o $@ $<

%.o: %.mod

%.o: %.f90
	$(FC) $(FCARGS) $(INC) $(LIB) -c $<

all: install

lib: $(TARGET).a $(TARGET).mod

install: $(TARGET).a $(TARGET).mod
	install --mode=755 -d $(INSTALLLIBDIR) $(INSTALLINCDIR)
	install --mode=644 $(TARGET).a $(INSTALLLIBDIR)/lib$(TARGET).a
	install --mode=644 $(TARGET).mod $(INSTALLINCDIR)

archive: $(TARGET).a

clean:
	rm -f *.o *.a *.mod









