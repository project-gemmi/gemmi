# This makefile is redundant. It is useful during development, but users
# are recommended to use the setup.py file from the top-level directory.

.SUFFIXES:

CXX=g++
CXXSTD=-std=c++11
PYTHON=python
PYTHON_CONFIG=$(PYTHON)-config

FLAGS=-O2 -g -pipe $(CXXSTD) -I../include -I../third_party -fPIC \
      -fvisibility=hidden -fwrapv -D_FORTIFY_SOURCE=2 \
       -Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls -Wfloat-conversion $(EXTRA_WFLAGS)

all: gemmi.so

igdir=../include/gemmi

gemmi.o: gemmi.cpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

grid.o: grid.cpp $(igdir)/grid.hpp $(igdir)/unitcell.hpp \
       $(igdir)/symmetry.hpp $(igdir)/util.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

mm.o: mm.cpp $(igdir)/elem.hpp $(igdir)/unitcell.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

sym.o: sym.cpp $(igdir)/symmetry.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

cif.o: cif.cpp $(igdir)/cifdoc.hpp \
       $(igdir)/to_json.hpp $(igdir)/to_cif.hpp $(igdir)/util.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

cif_read.o: cif_read.cpp $(igdir)/cifdoc.hpp $(igdir)/cif.hpp \
            $(igdir)/numb.hpp $(igdir)/gz.hpp $(igdir)/util.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

OBJ=gemmi.o grid.o mm.o sym.o cif.o cif_read.o

gemmi.so: $(OBJ)
	$(CXX) $(FLAGS) -shared \
	-Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro $(OBJ) -lz -o $@

clean:
	rm -f *.so *.o

.PHONY: all clean
