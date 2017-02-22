CXX=g++-6
FLAGS=-O3 -g --std=c++11 -Wall -Wextra -Wdisabled-optimization -Wformat=2 \
      -Wredundant-decls -Wshadow -I../PEGTL

all: validate to_json

validate: validate.cc cif.hh ddl.hh
	$(CXX) $(FLAGS) $< -o $@
to_json: to_json.cc cif.hh
	$(CXX) $(FLAGS) $< -o $@

trace: validate.cc cif.hh
	$(CXX) -DCIF_VALIDATE_SHOW_TRACE $(FLAGS) $< -o $@

clean:
	rm -f validate to_json trace
