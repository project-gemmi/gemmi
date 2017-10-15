
.SUFFIXES:

#CXX=g++-7
CXXSTD=-std=c++11

WFLAGS=-Wall -Wextra -Wpedantic -Wdisabled-optimization -Wformat=2 \
       -Wredundant-decls -Wfloat-conversion $(EXTRA_WFLAGS)
FLAGS=-O2 -g -pipe $(CXXSTD) $(WFLAGS) -I../include -I../third_party #-DNDEBUG

PROGRAMS=gemmi-mask gemmi-map gemmi-validate gemmi-convert gemmi-grep gemmi-sg \
	 gemmi-contents

all: $(PROGRAMS)

igdir=../include/gemmi

gemmi-validate: validate.cpp options.h $(igdir)/cif.hpp \
                $(igdir)/ddl.hpp $(igdir)/gz.hpp $(igdir)/numb.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

input.o: input.cpp input.h $(igdir)/gz.hpp $(igdir)/cif.hpp \
         $(igdir)/numb.hpp $(igdir)/cifdoc.hpp $(igdir)/mmcif.hpp \
         $(igdir)/model.hpp $(igdir)/pdb.hpp $(igdir)/mmread.hpp \
	 $(igdir)/util.hpp $(igdir)/json.hpp $(igdir)/unitcell.hpp
	$(CXX) $(FLAGS) -c $<

output.o: output.cpp output.h $(igdir)/to_pdb.hpp $(igdir)/model.hpp \
          $(igdir)/elem.hpp $(igdir)/unitcell.hpp $(igdir)/to_mmcif.hpp \
          $(igdir)/cifdoc.hpp $(igdir)/sprintf.hpp
	$(CXX) $(FLAGS) -c $<

convert.o: convert.cpp options.h input.h $(igdir)/to_cif.hpp \
           $(igdir)/to_json.hpp $(igdir)/to_mmcif.hpp $(igdir)/cifdoc.hpp \
           $(igdir)/util.hpp $(igdir)/numb.hpp $(igdir)/model.hpp
	$(CXX) $(FLAGS) -c $<

gemmi-convert: convert.o input.o output.o
	$(CXX) $(FLAGS) $< input.o output.o -o $@ -lz

gemmi-grep: grep.cpp options.h $(igdir)/cif.hpp $(igdir)/gz.hpp
	$(CXX) $(FLAGS) $< -o $@ -lz

gemmi-mask: mask.cpp options.h input.o $(igdir)/grid.hpp \
            $(igdir)/unitcell.hpp $(igdir)/model.hpp
	$(CXX) $(FLAGS) $< input.o -o $@ -lz

gemmi-map: map.cpp options.h input.o $(igdir)/grid.hpp \
           $(igdir)/unitcell.hpp
	$(CXX) $(FLAGS) $< input.o -o $@ -lz

gemmi-sg: sg.cpp $(igdir)/symmetry.hpp
	$(CXX) $(FLAGS) $< -o $@

gemmi-contents: contents.cpp options.h input.o
	$(CXX) $(FLAGS) $< input.o -o $@ -lz

write-help: $(PROGRAMS)
	../tools/docs-help.sh

clean:
	rm -f $(PROGRAMS) trace *.so *.o

.PHONY: all clean write-help
