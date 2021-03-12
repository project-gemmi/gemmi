# includes part of shroud/regression/run/defaults.mk

SUFFIX = -10
CC = gcc$(SUFFIX)
CXX = g++$(SUFFIX)
FC = gfortran$(SUFFIX)

CFLAGS = -O -Wall -Wstrict-prototypes -fno-strict-aliasing -std=c99
FFLAGS = -O -cpp -std=f2003 -ffree-form -Wall -fbounds-check
CXXFLAGS = -O -Wall -std=c++11 -fno-strict-aliasing
INCLUDE = -I../include
FLIBS = -lstdc++

C_OBJS = \
    wrapPosition.o \
    wrapFractional.o \
    wrapUnitCell.o \
    wrapSpaceGroup.o \
    wrapGroupOps.o \
    wrapOp.o \
    wrapCcp4_float.o \
    wrapCcp4_int8_t.o \
    wrapMtz.o \
    wrapgemmi.o \
    utilgemmi.o

F_OBJS = wrapfgemmi.o

all: wrapfgemmi.f ftest fsym

clean:
	rm -f ftest *.o *.mod
	rm -f wrap* types*.h util*.cpp *.json *.log gemmi_types.yaml

wrapfgemmi.f: gemmi.yaml
	shroud gemmi.yaml

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $*.o $<

%.o %.mod: %.f
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $<

%.o %.mod: %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $<

ftest: $(C_OBJS) $(F_OBJS) ftest.o
	$(FC) $(FFLAGS) $^ -o $@ $(FLIBS)

fsym: $(C_OBJS) $(F_OBJS) fsym.o
	$(FC) $(FFLAGS) $^ -o $@ $(FLIBS)

.PHONY: clean all
