cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(gemmi C CXX)

option(USE_FORTRAN "Build Fortran bindings" OFF)
option(USE_PYTHON "Build Python bindings" OFF)
option(INSTALL_EGG_INFO "Install .egg-info via setup.py" ON)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
option(USE_WMAIN "(Windows only) take Unicode arguments in gemmi program" ON)

if (DEFINED ENV{FC} OR CMAKE_Fortran_COMPILER)
  set(USE_FORTRAN ON CACHE BOOL "Build Fortran bindings" FORCE)
endif()

if (INFO)
  set(GEMMI_VERSION_INFO ${INFO} CACHE STRING "Extra text for gemmi -V" FORCE)
endif()

if (USE_FORTRAN)
  enable_language(Fortran)
else()
  message(STATUS "Skipping Fortran bindings. Add -D USE_FORTRAN=1 to build them.")
endif()

if (DEFINED ENV{EXTRA_WFLAGS})
  set(EXTRA_WARNINGS ON CACHE BOOL "Set extra warning flags" FORCE)
endif()

if (NOT CMAKE_CXX_STANDARD)
  if (CMAKE_CXX17_STANDARD_COMPILE_OPTION)
    # Python bindings don't compile as C++17 in VS 2017
    if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
      set(CMAKE_CXX_STANDARD 14)
    else()
      set(CMAKE_CXX_STANDARD 17)
    endif()
  elseif (CMAKE_CXX11_STANDARD_COMPILE_OPTION)
    set(CMAKE_CXX_STANDARD 11)
  endif()
endif()
message(STATUS "Compiling with C++ standard: ${CMAKE_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

# avoid CMake warning about unset policy
if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if (POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
endif()

find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories("${ZLIB_INCLUDE_DIR}")
else()
  message(STATUS "The build will use zlib code from third_party/zlib.")
  include_directories("${CMAKE_SOURCE_DIR}/third_party/zlib")
endif()
find_package(benchmark QUIET)
if (benchmark_FOUND)
  message(STATUS "Found benchmark: ${benchmark_DIR}")
else (NOT benchmark_FOUND)
  message(STATUS "Benchmarks not configured.")
endif()

include_directories("${CMAKE_SOURCE_DIR}/include"
                    "${CMAKE_SOURCE_DIR}/third_party")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 -Wredundant-decls -Wfloat-conversion -Wdisabled-optimization -Wshadow $ENV{EXTRA_WFLAGS}")
  message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" cxx_flags_config)
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS} ${${cxx_flags_config}}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qvec-report:1")
endif()

if (USE_FORTRAN)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2003 -fimplicit-none -Wall -Wextra -pedantic")
    message(STATUS "Fortran flags set to: ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check")
  endif()
endif()


if (ZLIB_FOUND)
  macro(support_gz exe)
    target_link_libraries(${exe} PRIVATE ZLIB::ZLIB)
  endmacro()
else()
  add_library(ungz OBJECT
      "third_party/zlib/adler32.c"
      "third_party/zlib/crc32.c"
      "third_party/zlib/gzlib.c"
      "third_party/zlib/gzread.c"
      "third_party/zlib/inflate.c"
      "third_party/zlib/inftrees.c"
      "third_party/zlib/inffast.c"
      "third_party/zlib/zutil.c")
  check_include_file(unistd.h has_unistd_h)
  target_compile_definitions(ungz PRIVATE NO_GZCOMPRESS=1)
  if (has_unistd_h)
    target_compile_definitions(ungz PRIVATE Z_HAVE_UNISTD_H=1)
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /wd4996 - the POSIX name for this item is deprecated.
    # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
    target_compile_options(ungz PRIVATE /wd4996 /wd4267)
  endif()
  if (USE_PYTHON)
    set_property(TARGET ungz PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif()
  macro(support_gz exe)
    target_sources(${exe} PUBLIC $<TARGET_OBJECTS:ungz>)
  endmacro()
endif()

if (WIN32 AND USE_WMAIN)
  add_definitions(-D_UNICODE=1)
endif()

#add_library(cgemmi STATIC fortran/grid.cpp fortran/symmetry.cpp)
#
#if (USE_FORTRAN)
#  add_library(fgemmi STATIC fortran/gemmi.f90)
#  target_link_libraries(fgemmi PRIVATE cgemmi)
#endif()


### programs from src/ ###

add_library(input OBJECT src/input.cpp)
add_library(options OBJECT src/options.cpp)
if (GEMMI_VERSION_INFO)
  target_compile_definitions(options PRIVATE GEMMI_VERSION_INFO=${GEMMI_VERSION_INFO})
endif()
add_library(output OBJECT src/output.cpp)
add_library(mapcoef OBJECT src/mapcoef.cpp)

add_executable(gemmi-align EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/align.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-align)

add_executable(gemmi-blobs EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/blobs.cpp $<TARGET_OBJECTS:mapcoef> $<TARGET_OBJECTS:input>)
support_gz(gemmi-blobs)

add_executable(gemmi-cif2json EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/cif2json.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-cif2json)

add_executable(gemmi-cif2mtz EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/cif2mtz.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-cif2mtz)

add_executable(gemmi-contact EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/contact.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-contact)

add_executable(gemmi-contents EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/contents.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-contents)

add_executable(gemmi-convert EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/convert.cpp $<TARGET_OBJECTS:output> $<TARGET_OBJECTS:input>)
support_gz(gemmi-convert)

add_executable(gemmi-crdrst EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/crdrst.cpp $<TARGET_OBJECTS:output>
               $<TARGET_OBJECTS:input>)
support_gz(gemmi-crdrst)

add_executable(gemmi-fprime EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/fprime.cpp)

add_executable(gemmi-grep EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/grep.cpp)
support_gz(gemmi-grep)

add_executable(gemmi-h EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/h.cpp $<TARGET_OBJECTS:output> $<TARGET_OBJECTS:input>)
support_gz(gemmi-h)

add_executable(gemmi-json2cif EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/json2cif.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-json2cif)

add_executable(gemmi-map EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/map.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-map)

add_executable(gemmi-map2sf EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/map2sf.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-map2sf)

add_executable(gemmi-mask EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/mask.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-mask)

add_executable(gemmi-mixmtz EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/mixmtz.cpp $<TARGET_OBJECTS:output>)

add_executable(gemmi-merge EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/merge.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-merge)

add_executable(gemmi-mondiff EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/mondiff.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-mondiff)

add_executable(gemmi-mtz EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/mtz.cpp)
support_gz(gemmi-mtz)

add_executable(gemmi-mtz2cif EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/mtz2cif.cpp)
support_gz(gemmi-mtz2cif)

add_executable(gemmi-reindex EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/reindex.cpp)
support_gz(gemmi-reindex)

add_executable(gemmi-residues EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/residues.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-residues)

add_executable(gemmi-rmsz EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/rmsz.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-rmsz)

add_executable(gemmi-sf2map EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/sf2map.cpp $<TARGET_OBJECTS:mapcoef> $<TARGET_OBJECTS:input>)
support_gz(gemmi-sf2map)

add_executable(gemmi-sfcalc EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/sfcalc.cpp $<TARGET_OBJECTS:input> $<TARGET_OBJECTS:output>)
support_gz(gemmi-sfcalc)

add_executable(gemmi-sg EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/sg.cpp)

add_executable(gemmi-tags EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/tags.cpp)
support_gz(gemmi-tags)

add_executable(gemmi-validate EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/validate.cpp src/validate_mon.cpp)
support_gz(gemmi-validate)

add_executable(gemmi-wcn EXCLUDE_FROM_ALL $<TARGET_OBJECTS:options>
               src/wcn.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-wcn)

add_executable(gemmi-prog
               src/align.cpp src/blobs.cpp
               src/cif2mtz.cpp src/cif2json.cpp src/contact.cpp
               src/contents.cpp src/convert.cpp src/fprime.cpp
               src/grep.cpp src/h.cpp src/json2cif.cpp
               src/main.cpp src/map.cpp src/map2sf.cpp src/mask.cpp
               src/merge.cpp src/mondiff.cpp src/mtz.cpp src/mtz2cif.cpp
               src/reindex.cpp src/residues.cpp src/rmsz.cpp
               src/sf2map.cpp src/sfcalc.cpp src/sg.cpp
               src/tags.cpp src/validate.cpp src/validate_mon.cpp src/wcn.cpp
               $<TARGET_OBJECTS:mapcoef>
               $<TARGET_OBJECTS:input>
               $<TARGET_OBJECTS:output>
               $<TARGET_OBJECTS:options>)
support_gz(gemmi-prog)
target_compile_definitions(gemmi-prog PRIVATE GEMMI_ALL_IN_ONE=1)
set_target_properties(gemmi-prog PROPERTIES OUTPUT_NAME gemmi)
if (WIN32 AND USE_WMAIN)
  # _UNICODE=1 is now set globally
  #target_compile_definitions(gemmi-prog PRIVATE _UNICODE=1)
  if(MINGW)
    # target_link_options were added in cmake 3.13
    set_target_properties(gemmi-prog PROPERTIES LINK_FLAGS "-municode")
  endif()
endif()

### tests and examples ###

#add_executable(c_test EXCLUDE_FROM_ALL fortran/c_test.c)
#target_link_libraries(c_test PRIVATE cgemmi)

add_executable(cpptest EXCLUDE_FROM_ALL tests/main.cpp tests/cif.cpp)

add_executable(hello EXCLUDE_FROM_ALL examples/hello.cpp)
add_executable(doc_example EXCLUDE_FROM_ALL
               docs/code/sym.cpp docs/code/elem.cpp docs/code/resinfo.cpp
               docs/code/cell.cpp)
add_executable(doc_example2 EXCLUDE_FROM_ALL docs/code/cif_cc.cpp)
add_executable(doc_maybegz EXCLUDE_FROM_ALL
               docs/code/maybegz.cpp
               docs/code/mutate.cpp)
support_gz(doc_maybegz)
add_executable(doc_newmtz EXCLUDE_FROM_ALL docs/code/newmtz.cpp)

# always compile these tests with assertions enabled
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
  target_compile_options(doc_example PRIVATE -UNDEBUG)
endif()

add_executable(test_disulf EXCLUDE_FROM_ALL tests/disulf.cpp)
support_gz(test_disulf)

# auth_label requires <experimental/filesystem> and -lstdc++fs
add_executable(auth_label EXCLUDE_FROM_ALL examples/auth_label.cpp)
if (NOT MSVC)
  target_link_libraries(auth_label PRIVATE stdc++fs)
endif()
support_gz(auth_label)

add_executable(check_conn EXCLUDE_FROM_ALL examples/check_conn.cpp)
support_gz(check_conn)

enable_testing()

add_custom_target(print-version
  COMMAND gemmi-prog --version
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "gemmi --version"
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION>)
add_test(NAME cpptest COMMAND cpptest)

add_dependencies(check
    cpptest hello doc_example doc_example2 doc_maybegz doc_newmtz
    test_disulf check_conn print-version)

if (USE_FORTRAN)
#  add_executable(ftest EXCLUDE_FROM_ALL fortran/ftest.f90)
#  target_link_libraries(ftest PRIVATE fgemmi)
#  add_test(NAME ftest COMMAND ftest)
#  add_executable(ftest_grid EXCLUDE_FROM_ALL fortran/ftest_grid.f90)
#  target_link_libraries(ftest_grid PRIVATE fgemmi)
#  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
#    set_property(TARGET ftest ftest_grid PROPERTY LINKER_LANGUAGE Fortran)
#  endif()
#  add_test(NAME ftest_grid COMMAND ftest_grid)
#  add_dependencies(check ftest ftest_grid)

endif()

### benchmarks ###

if (benchmark_FOUND)
  foreach(b stoi elem mod pdb resinfo round sym)
    add_executable(${b}-bm EXCLUDE_FROM_ALL benchmarks/${b}.cpp)
    target_link_libraries(${b}-bm PRIVATE benchmark::benchmark)
    set_target_properties(${b}-bm PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                             "${CMAKE_BINARY_DIR}/benchmarks")
    add_dependencies(check ${b}-bm)
  endforeach()
endif()

### Python bindings ###
# Alternatively, the Python module can be built with setup.py.

if (USE_PYTHON)
  message(STATUS "The python module will be built.")
  if (EXISTS ${CMAKE_HOME_DIRECTORY}/pybind11)
    message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/pybind11 (internal copy).")
    add_subdirectory(pybind11)
  else()
    find_package(pybind11 2.6 CONFIG REQUIRED)
    message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
  endif()
  pybind11_add_module(gemmi python/monlib.cpp
                            python/gemmi.cpp python/align.cpp
                            python/ccp4.cpp python/cif.cpp python/custom.cpp
                            python/elem.cpp python/grid.cpp python/hkl.cpp
                            python/meta.cpp python/mol.cpp
                            python/mtz.cpp python/read.cpp python/recgrid.cpp
                            python/scaling.cpp python/search.cpp
                            python/sf.cpp python/sym.cpp python/topo.cpp
                            python/unitcell.cpp python/write.cpp)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(gemmi PRIVATE "-Wno-shadow")
  endif()
  # avoid GCC warning: the ABI of passing structure with 'complex float' member
  # has changed in GCC 4.4
  set_source_files_properties(python/recgrid.cpp python/hkl.cpp python/mtz.cpp
                              src/mtz.cpp
                              PROPERTIES COMPILE_FLAGS
                              $<$<CXX_COMPILER_ID:GNU>:-Wno-psabi>)
  support_gz(gemmi)
  if (NOT DEFINED PYTHON_INSTALL_DIR)
    set(PYTHON_INSTALL_DIR ${PYTHON_SITE_PACKAGES})
  endif()
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()

install(TARGETS gemmi-prog DESTINATION bin)
install(DIRECTORY include/gemmi DESTINATION include)
if (USE_PYTHON AND DEFINED PYTHON_INSTALL_DIR)
  install(TARGETS gemmi DESTINATION ${PYTHON_INSTALL_DIR})
  install(DIRECTORY examples/ DESTINATION ${PYTHON_INSTALL_DIR}/gemmi-examples
          FILES_MATCHING PATTERN "*.py")
  if (INSTALL_EGG_INFO)
    install(CODE
     "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" setup.py install_egg_info --install-dir \"\$ENV{DESTDIR}${PYTHON_INSTALL_DIR}\"
                      WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR}\"
                      RESULT_VARIABLE EGG_STATUS)"
    )
    if (EGG_STATUS AND NOT EGG_STATUS EQUAL 0)
      message(FATAL_ERROR "Failed to install egg-info. Use -D INSTALL_EGG_INFO=OFF to disable.")
    endif()
  endif()
endif()
