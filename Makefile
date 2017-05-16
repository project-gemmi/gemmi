#CXX=clang++ -stdlib=libc++
#CXX=/home/wojdyr/local/clang40/bin/clang++ -L /home/wojdyr/local/clang40/lib -stdlib=libc++ -Wl,-rpath,/home/wojdyr/local/clang40/lib
CXX=g++
WFLAGS=-Wall -Wextra -Wdisabled-optimization -Wformat=2 -Wredundant-decls
FLAGS=-O2 -g --std=c++11 $(WFLAGS) -Wshadow -Ithird_party #-DNDEBUG
PYFLAGS=-O2 -g --std=c++14 $(WFLAGS) -Ithird_party -fvisibility=hidden \
	-fwrapv -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fPIC

all: gemmi-validate gemmi-convert gemmi.so

gemmi-validate: validate.cc cif.hh ddl.hh cifgz.hh numb.hh
	$(CXX) $(FLAGS) $< -o $@ -lz
gemmi-convert: convert.cc to_json.hh cif.hh write_cif.hh mmcif.hh model.hh \
               to_pdb.hh unitcell.hh read_pdb.hh elem.hh util.hh
	$(CXX) $(FLAGS) -Wno-strict-aliasing $< -o $@ -lz

# for debugging only
trace: validate.cc cif.hh
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@

pygemmi.o: pygemmi.cc cif.hh to_json.hh numb.hh write_cif.hh elem.hh
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

clean:
	rm -f gemmi-validate gemmi-convert trace gemmi.so pygemmi.o
