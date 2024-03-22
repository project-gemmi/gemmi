# currently using Emscripten 3.1.56
CXX = em++

INCLUDE = -I../include

FLAGS = \
  -O3 -flto -Wall -Wextra -std=c++20 \
  -DEMSCRIPTEN_HAS_UNBOUND_TYPE_NAMES=0 -fno-rtti \
  -s STRICT=1 \
  -s DISABLE_EXCEPTION_CATCHING=0

LINK_FLAGS = \
  --pre-js pre.js \
  --post-js post.js \
  --no-entry \
  -lembind \
  -s EXPORTED_RUNTIME_METHODS=writeArrayToMemory,HEAPU8 \
  -s EXPORTED_FUNCTIONS=_malloc \
  -s INCOMING_MODULE_JS_API=print,printErr,setStatus,onRuntimeInitialized \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=Gemmi

#  -s ASSERTIONS=1 -s SAFE_HEAP=1 \

GEMMI_OBJS = mmcif.o polyheur.o resinfo.o sprintf.o json.o

BINDING_OBJS = gemmi.o cell.o mol.o mtz_fft.o

gemmi.js: $(GEMMI_OBJS) $(BINDING_OBJS)
	$(CXX) $(FLAGS) $(LINK_FLAGS) $(GEMMI_OBJS) $(BINDING_OBJS) -o $@

# old mtz module (not built by default)
mtz.js: mtz_fft.cpp cell.o
	$(CXX) -DSTANDALONE_MTZ=1 $(INCLUDE) $(FLAGS) \
	-lembind --no-entry --pre-js=pre.js \
	-s MODULARIZE=1 -s ALLOW_MEMORY_GROWTH=1 \
	-s EXPORTED_FUNCTIONS=_malloc \
	-s EXPORTED_RUNTIME_METHODS=writeArrayToMemory \
	-s EXPORT_NAME=GemmiMtz $^ -o $@

$(GEMMI_OBJS): %.o: ../src/%.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $< -o $@

$(BINDING_OBJS): %.o: %.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(INCLUDE) $(FLAGS) -c $<

clean:
	rm -f *.o gemmi.js mtz.js *.wasm

.PHONY: clean
