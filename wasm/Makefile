CXX = em++

INCLUDE = -I../include

FLAGS = \
  -O3 -flto -Wall -Wextra -std=c++20 \
  -DEMSCRIPTEN_HAS_UNBOUND_TYPE_NAMES=0 -fno-rtti \
  -s STRICT=1 \
  -s DISABLE_EXCEPTION_CATCHING=0

LINK_FLAGS = \
  --pre-js pre.js \
  -s EXPORTED_RUNTIME_METHODS=ccall,writeArrayToMemory,getValue,UTF8ToString \
  -s EXPORTED_FUNCTIONS=_malloc \
  -s INCOMING_MODULE_JS_API=print,printErr,setStatus,onRuntimeInitialized \
  -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1

#  -s ASSERTIONS=1 -s SAFE_HEAP=1 \

GEMMI_OBJS = mmcif.o polyheur.o resinfo.o sprintf.o json.o

BINDING_OBJS = mol.o

gemmi.js: $(GEMMI_OBJS) $(BINDING_OBJS)
	$(CXX) --bind $(FLAGS) $(LINK_FLAGS) $(GEMMI_OBJS) $(BINDING_OBJS) -o $@

$(GEMMI_OBJS): %.o: ../src/%.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $< -o $@

$(BINDING_OBJS): %.o: %.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $<

clean:
	rm -f *.o gemmi.js *.wasm

.PHONY: clean
