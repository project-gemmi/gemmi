
.SUFFIXES:

#CXX=g++-7
CXXSTD=-std=c++11
PYTHON=python
PYTHON_CONFIG=$(PYTHON)-config

WFLAGS=-Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls -Wshadow $(EXTRA_WFLAGS)
FLAGS=-O2 -g -pipe $(CXXSTD) $(WFLAGS) -Iinclude -Ithird_party #-DNDEBUG

PYFLAGS=$(FLAGS) -Wno-shadow -fPIC \
       -fvisibility=hidden -fwrapv -D_FORTIFY_SOURCE=2

all: gemmi-validate gemmi-convert gemmi-grep

py: gemmi.so

gemmi-validate: validate.cpp options.h include/gemmi/cif.hpp \
                include/gemmi/ddl.hpp include/gemmi/cifgz.hpp \
                include/gemmi/numb.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

gemmi-convert: convert.cpp options.h include/gemmi/*.hpp
	$(CXX) $(FLAGS) -Wno-strict-aliasing $< -o $@ -lz

gemmi-convert-snprintf: convert.cpp options.h include/gemmi/*.hpp
	$(CXX) -DUSE_STD_SNPRINTF $(FLAGS) $< -o $@ -lz

gemmi-grep: grep.cpp options.h include/gemmi/cif.hpp include/gemmi/cifgz.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

gemmi-mask: mask.cpp options.h include/gemmi/grid.hpp include/gemmi/unitcell.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

# for debugging only
trace: validate.cpp options.h include/gemmi/cif.hpp
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@ -lz

pygemmi.o: pygemmi.cpp include/gemmi/cif.hpp include/gemmi/to_json.hpp \
           include/gemmi/numb.hpp include/gemmi/to_cif.hpp \
	   include/gemmi/elem.hpp
	$(CXX) $(PYFLAGS) `$(PYTHON_CONFIG) --includes` -c $<

gemmi.so: pygemmi.o
	$(CXX) $(PYFLAGS) -shared \
	-Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro $< -o $@

write-help: gemmi-validate gemmi-convert
	./gemmi-convert tests/misc.cif tests/misc.json
	echo '$$ gemmi-validate -h' > docs/validate-help.txt
	./gemmi-validate -h >> docs/validate-help.txt
	echo '$$ gemmi-convert -h' > docs/convert-help.txt
	./gemmi-convert -h >> docs/convert-help.txt
	echo '$$ gemmi-grep -h' > docs/grep-help.txt
	./gemmi-grep -h >> docs/grep-help.txt

clean:
	rm -f gemmi-validate gemmi-convert trace gemmi.so pygemmi.o

.PHONY: all py clean write-help
