
CXX=g++

WFLAGS=-Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls -Wshadow
FLAGS=-O2 -g --std=c++11 $(WFLAGS) -Iinclude -Ithird_party #-DNDEBUG

PYFLAGS=$(FLAGS) -Wno-shadow -fPIC \
       -fvisibility=hidden -fwrapv -D_FORTIFY_SOURCE=2 -fstack-protector-strong

all: gemmi-validate gemmi-convert gemmi-grep

py: gemmi.so

gemmi-validate: validate.cpp include/gemmi/cif.hpp include/gemmi/ddl.hpp \
                include/gemmi/cifgz.hpp include/gemmi/numb.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

gemmi-convert: convert.cpp include/gemmi/*.hpp
	$(CXX) $(FLAGS) -Wno-strict-aliasing $< -o $@ -lz

gemmi-convert-snprintf: convert.cpp include/gemmi/*.hpp
	$(CXX) -DUSE_STD_SNPRINTF $(FLAGS) $< -o $@ -lz

gemmi-grep: grep.cpp include/gemmi/cif.hpp include/gemmi/cifgz.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

# for debugging only
trace: validate.cpp include/gemmi/cif.hpp
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@ -lz

pygemmi.o: pygemmi.cpp include/gemmi/cif.hpp include/gemmi/to_json.hpp \
           include/gemmi/numb.hpp include/gemmi/to_cif.hpp \
	   include/gemmi/elem.hpp
	$(CXX) $(PYFLAGS) -I/usr/include/python2.7 -c $<

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
