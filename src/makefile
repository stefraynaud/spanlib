TARGET=spanlib
TOPDIR=..
SRCDIR=$(TOPDIR)/src
INSTALLDIR=$(HOME)
INSTALLLIBDIR=$(HOME)/lib
INSTALLINCDIR=$(HOME)/include
FC=ifort
FCARGS=-132 -O -check bounds
##INC=-I/home/users/raynaud/lib/lapack95_modules
LIBS=-L/usr/local/install/lapack-3.0/lib/ -llapack -lblas -llapack -lblas
LIB=lib$(TARGET).a
OBJ=$(TARGET).o
MOD=$(TARGET).mod 

######################################################################

.SUFFIXES: .f90 .a .o .mod 

##%.a: %.o
$(LIB): $(OBJ)
	ar -r $@ $<
##	ld -o $@ $<

%.o: %.mod

%.o: %.f90
	$(FC) $(FCARGS) $(INC) $(LIBS) -c $<

all: install

lib: $(LIB) $(MOD)

install: $(LIB) $(MOD)
	install --mode=755 -d $(INSTALLLIBDIR) $(INSTALLINCDIR)
	install --mode=644 $(LIB) $(INSTALLLIBDIR)
	install --mode=644 $(MOD) $(INSTALLINCDIR)

archive: $(LIB)

clean:
	rm -f *.o *.a *.mod









