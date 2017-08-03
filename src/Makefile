
.SUFFIXES:

#CXX=g++-7
CXXSTD=-std=c++11
PYTHON=python
PYTHON_CONFIG=$(PYTHON)-config

WFLAGS=-Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls $(EXTRA_WFLAGS)
FLAGS=-O2 -g -pipe $(CXXSTD) $(WFLAGS) -I../include -I../third_party #-DNDEBUG

PYFLAGS=$(FLAGS) -Wno-shadow -fPIC \
       -fvisibility=hidden -fwrapv -D_FORTIFY_SOURCE=2

PROGRAMS=gemmi-mask gemmi-map gemmi-validate gemmi-convert gemmi-grep

all: $(PROGRAMS)

py: gemmi.so

igdir=../include/gemmi

gemmi-validate: validate.cpp options.h $(igdir)/cif.hpp \
                $(igdir)/ddl.hpp $(igdir)/cifgz.hpp $(igdir)/numb.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

input.o: input.cpp input.h $(igdir)/cifgz.hpp $(igdir)/cif.hpp \
         $(igdir)/numb.hpp $(igdir)/cifdoc.hpp $(igdir)/mmcif.hpp \
         $(igdir)/model.hpp $(igdir)/pdb.hpp $(igdir)/pdbgz.hpp \
	 $(igdir)/util.hpp
	$(CXX) $(FLAGS) -c $<

output.o: output.cpp output.h $(igdir)/to_pdb.hpp $(igdir)/model.hpp \
          $(igdir)/elem.hpp $(igdir)/unitcell.hpp $(igdir)/to_mmcif.hpp \
          $(igdir)/cifdoc.hpp $(igdir)/sprintf.hpp
	$(CXX) $(FLAGS) -Wno-strict-aliasing -c $<

convert.o: convert.cpp options.h input.h $(igdir)/to_cif.hpp \
           $(igdir)/to_json.hpp $(igdir)/to_mmcif.hpp $(igdir)/cifdoc.hpp \
           $(igdir)/util.hpp $(igdir)/numb.hpp $(igdir)/model.hpp
	$(CXX) $(FLAGS) -Wno-strict-aliasing -c $<

gemmi-convert: convert.o input.o output.o
	$(CXX) $(FLAGS) $< input.o output.o -o $@ -lz

gemmi-grep: grep.cpp options.h $(igdir)/cif.hpp $(igdir)/cifgz.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

gemmi-mask: mask.cpp options.h input.h input.o $(igdir)/grid.hpp \
            $(igdir)/unitcell.hpp $(igdir)/model.hpp
	$(CXX) $(FLAGS) $< input.o -o $@ -lz

gemmi-map: map.cpp options.h input.h input.o $(igdir)/grid.hpp \
           $(igdir)/unitcell.hpp
	$(CXX) $(FLAGS) $< input.o -o $@ -lz

# for debugging only
trace: validate.cpp options.h $(igdir)/cif.hpp
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@ -lz

pygemmi.o: ../pygemmi.cpp $(igdir)/cif.hpp $(igdir)/to_json.hpp \
           $(igdir)/numb.hpp $(igdir)/to_cif.hpp $(igdir)/elem.hpp
	$(CXX) $(PYFLAGS) `$(PYTHON_CONFIG) --includes` -c $<

gemmi.so: pygemmi.o
	$(CXX) $(PYFLAGS) -shared \
	-Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro $< -o $@

write-help: $(PROGRAMS)
	./gemmi-convert ../tests/misc.cif ../tests/misc.json
	echo '$$ gemmi-validate -h' > ../docs/validate-help.txt
	./gemmi-validate -h >> ../docs/validate-help.txt
	echo '$$ gemmi-convert -h' > ../docs/convert-help.txt
	./gemmi-convert -h >> ../docs/convert-help.txt
	echo '$$ gemmi-grep -h' > ../docs/grep-help.txt
	./gemmi-grep -h >> ../docs/grep-help.txt

clean:
	rm -f $(PROGRAMS) trace *.so *.o

.PHONY: all py clean write-help
