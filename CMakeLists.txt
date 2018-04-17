cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(gemmi C CXX)

option(USE_FORTRAN "Build Fortran bindings" OFF)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)

if (DEFINED ENV{FC} OR CMAKE_Fortran_COMPILER)
  set(USE_FORTRAN ON CACHE BOOL "Build Fortran bindings" FORCE)
endif()

if (USE_FORTRAN)
  enable_language(Fortran)
else()
  message(STATUS
          "Skipping Fortran bindings. If you need them add -D USE_FORTRAN=1")
endif()

if (DEFINED ENV{EXTRA_WFLAGS})
    set(EXTRA_WARNINGS ON)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)

# Set default build mode (based on CMake FAQ)
if(NOT CMAKE_BUILD_TYPE AND NOT DEFINED ENV{CXXFLAGS})
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

find_package(ZLIB)
if (NOT ZLIB_FOUND)
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
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
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
    target_link_libraries(${exe} ZLIB::ZLIB)
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
  set(Z_DEF "-DNO_GZCOMPRESS")
  if (has_unistd_h)
    set(Z_DEF "${Z_DEF} -DZ_HAVE_UNISTD_H")
  endif()
  set_target_properties(ungz PROPERTIES COMPILE_FLAGS "${Z_DEF}")
  macro(support_gz exe)
    target_sources(${exe} PUBLIC $<TARGET_OBJECTS:ungz>)
  endmacro()
endif()

add_library(cgemmi STATIC fortran/grid.cpp fortran/symmetry.cpp)

if (USE_FORTRAN)
  add_library(fgemmi STATIC fortran/gemmi.f90)
  target_link_libraries(fgemmi cgemmi)
endif()


### programs from src/ ###

add_library(input OBJECT src/input.cpp)

add_executable(gemmi-contact EXCLUDE_FROM_ALL
               src/contact.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-contact)

add_executable(gemmi-contents EXCLUDE_FROM_ALL
               src/contents.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-contents)

add_executable(gemmi-convert EXCLUDE_FROM_ALL
               src/convert.cpp src/output.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-convert)

add_executable(gemmi-crdrst EXCLUDE_FROM_ALL
               src/crdrst.cpp src/output.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-crdrst)

add_executable(gemmi-grep EXCLUDE_FROM_ALL src/grep.cpp)
support_gz(gemmi-grep)

add_executable(gemmi-map EXCLUDE_FROM_ALL
               src/map.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-map)

add_executable(gemmi-mask EXCLUDE_FROM_ALL
               src/mask.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-mask)

add_executable(gemmi-sg EXCLUDE_FROM_ALL src/sg.cpp)

add_executable(gemmi-validate EXCLUDE_FROM_ALL src/validate.cpp)
support_gz(gemmi-validate)

add_executable(gemmi-prog
               src/input.cpp src/contact.cpp src/contents.cpp src/convert.cpp
               src/output.cpp src/grep.cpp src/main.cpp src/map.cpp
               src/mask.cpp src/sg.cpp src/validate.cpp)
support_gz(gemmi-prog)
target_compile_definitions(gemmi-prog PRIVATE GEMMI_ALL_IN_ONE=1)
set_target_properties(gemmi-prog PROPERTIES OUTPUT_NAME gemmi)

### tests and examples ###

add_executable(ctest EXCLUDE_FROM_ALL fortran/ctest.c)
target_link_libraries(ctest cgemmi)

add_executable(cpptest EXCLUDE_FROM_ALL tests/main.cpp tests/cif.cpp)

add_executable(hello EXCLUDE_FROM_ALL examples/hello.cc)
add_executable(doc_example EXCLUDE_FROM_ALL
               docs/code/sym.cpp docs/code/elem.cpp docs/code/resinfo.cpp
               docs/code/cell.cpp)
add_executable(doc_example2 EXCLUDE_FROM_ALL docs/code/cif_cc.cpp)
add_executable(doc_maybegz EXCLUDE_FROM_ALL
               docs/code/maybegz.cpp docs/code/structure.cpp
               docs/code/mutate.cpp)
support_gz(doc_maybegz)

# always compile these tests with assertions enabled
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
  target_compile_options(doc_example PRIVATE -UNDEBUG)
endif()

# auth_label requires <experimental/filesystem> and -lstdc++fs
add_executable(auth_label EXCLUDE_FROM_ALL examples/auth_label.cc)
target_link_libraries(auth_label stdc++fs)
support_gz(auth_label)

add_executable(check_conn EXCLUDE_FROM_ALL examples/check_conn.cc)
support_gz(check_conn)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION>)
add_test(NAME cpptest COMMAND cpptest)

add_dependencies(check
    ctest cpptest hello doc_example doc_example2 doc_maybegz check_conn)

if (USE_FORTRAN)
  add_executable(ftest EXCLUDE_FROM_ALL fortran/ftest.f90)
  target_link_libraries(ftest fgemmi)
  add_executable(ftest_grid EXCLUDE_FROM_ALL fortran/ftest_grid.f90)
  target_link_libraries(ftest_grid fgemmi)
  add_dependencies(check ftest ftest_grid)
endif()

### benchmarks ###

if (benchmark_FOUND)
  foreach(b bench_elem bench_mod bench_sym bench_round)
    add_executable(${b} EXCLUDE_FROM_ALL tools/${b}.cpp)
    target_link_libraries(${b} benchmark::benchmark)
    add_dependencies(check ${b})
  endforeach()
endif()

### Python bindings ###
# Note: setup.py is the primary way to build Python bindings.
# Here we build it only for convenience and only in certain configuration.
if (EXISTS "${CMAKE_SOURCE_DIR}/third_party/pybind11/CMakeLists.txt")
  message(STATUS "The python module will be built.")
  add_subdirectory(third_party/pybind11)
  pybind11_add_module(gemmi python/gemmi.cpp python/grid.cpp python/cif.cpp
                            python/unitcell.cpp python/mol.cpp
                            python/read.cpp python/sym.cpp)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(gemmi PRIVATE "-Wno-shadow")
  endif()
else()
  message(STATUS "Note: The python module can be built with setup.py.")
endif()
