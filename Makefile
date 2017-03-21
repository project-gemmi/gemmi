CXX=g++
FLAGS=-O3 -g --std=c++11 -Wall -Wextra -Wdisabled-optimization -Wformat=2 \
      -Wredundant-decls -Wshadow -Ithirdparty
PYFLAGS=-O3 -g --std=c++14 -Wall -Wformat=2 -Ithirdparty -fvisibility=hidden \
	-fwrapv -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fPIC

all: validate to_json mmcif gemmi.so

validate: validate.cc cif.hh ddl.hh cifgz.hh numb.hh
	$(CXX) $(FLAGS) $< -o $@ -lz
to_json: to_json.cc cif.hh
	$(CXX) $(FLAGS) $< -o $@

trace: validate.cc cif.hh
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@

mmcif: mmcif.cc mmcif.hh cif.hh cifgz.hh numb.hh
	$(CXX) $(FLAGS) $< -o $@ -lz

pygemmi.o: pygemmi.cc cif.hh
	$(CXX) $(PYFLAGS) -Ithirdparty -I/usr/include/python2.7 -c $<

gemmi.so: pygemmi.o
	$(CXX) $(PYFLAGS) -shared \
	-Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro $< -o $@

clean:
	rm -f validate to_json trace gemmi.so pygemmi.o
