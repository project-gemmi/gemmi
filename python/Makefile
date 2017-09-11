
.SUFFIXES:

#CXX=g++-7
CXXSTD=-std=c++11
PYTHON=python
PYTHON_CONFIG=$(PYTHON)-config

FLAGS=-O2 -g -pipe $(CXXSTD) -I../include -I../third_party -fPIC \
      -fvisibility=hidden -fwrapv -D_FORTIFY_SOURCE=2 \
       -Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls -Wfloat-conversion $(EXTRA_WFLAGS)

py: gemmi.so

igdir=../include/gemmi

gemmi.o: gemmi.cpp $(igdir)/elem.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

cif.o: cif.cpp $(igdir)/cifdoc.hpp $(igdir)/cif.hpp $(igdir)/numb.hpp \
       $(igdir)/gz.hpp $(igdir)/to_json.hpp $(igdir)/to_cif.hpp
	$(CXX) $(FLAGS) `$(PYTHON_CONFIG) --includes` -c $<

gemmi.so: gemmi.o cif.o
	$(CXX) $(FLAGS) -shared \
	-Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro $< cif.o -o $@

clean:
	rm -f *.so *.o

.PHONY: all py clean
