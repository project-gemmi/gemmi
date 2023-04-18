# includes part of shroud/regression/run/defaults.mk

SUFFIX = -12
CC = gcc$(SUFFIX)
CXX = g++$(SUFFIX)
FC = gfortran$(SUFFIX)

CFLAGS = -O -Wall -Wstrict-prototypes -fno-strict-aliasing -std=c99
FFLAGS = -O -cpp -std=f2003 -ffree-form -Wall -fbounds-check
CXXFLAGS = -O -Wall -fno-strict-aliasing
INCLUDE = -I../include
FLIBS = -lstdc++ -lz

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

TEST_PROGS = ftest fsym

all: $(TEST_PROGS)

clean:
	rm -f $(TEST_PROGS) *.a *.o *.mod
	rm -f wrap* types*.h util*.cpp *.json *.log gemmi_types.yaml

$(C_OBJS:.o=.cpp) $(F_OBJS:.o=.f90): gemmi.yaml
	shroud gemmi.yaml

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

%.o %.mod: %.f
	$(FC) $(FFLAGS) $(INCLUDE) -c $<

%.o: %.f90 wrapfgemmi.o
	$(FC) $(FFLAGS) $(INCLUDE) -c $<

gemmi.a: $(C_OBJS) $(F_OBJS) write.o
	ar rcs gemmi.a $^

ftest: ftest.o gemmi.a
	$(FC) $(FFLAGS) $^ -o $@ $(FLIBS)

fsym: fsym.o gemmi.a
	$(FC) $(FFLAGS) $^ -o $@ $(FLIBS)

.PHONY: clean all
