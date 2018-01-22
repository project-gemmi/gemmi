cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(gemmi C CXX)

option(USE_FORTRAN "Build Fortran bindings" OFF)

if (DEFINED ENV{FC})
  set(USE_FORTRAN ON)
endif()
if (CMAKE_Fortran_COMPILER)
  set(USE_FORTRAN ON)
endif()

if (USE_FORTRAN)
  enable_language(Fortran)
else()
  message(STATUS
          "Skipping Fortran bindings. If you need them add -D USE_FORTRAN=1")
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

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(CMAKE_CXX_FLAGS
   "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 $ENV{EXTRA_WFLAGS}")
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS}")
  # avoid "unused" warnings when compiling with -DNDEBUG
  set(extra_release_flags "-Wno-unused-parameter -Wno-unused-variable")
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


set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${extra_release_flags}")
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} ${extra_release_flags}")

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

add_library(input OBJECT src/input.cpp)

add_executable(gemmi-validate src/validate.cpp)
support_gz(gemmi-validate)

add_executable(gemmi-grep src/grep.cpp)
support_gz(gemmi-grep)

add_executable(gemmi-convert src/convert.cpp src/output.cpp
                             $<TARGET_OBJECTS:input>)
support_gz(gemmi-convert)

add_executable(gemmi-map src/map.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-map)

add_executable(gemmi-mask src/mask.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-mask)

add_executable(gemmi-sg src/sg.cpp)

add_executable(gemmi-contents src/contents.cpp $<TARGET_OBJECTS:input>)
support_gz(gemmi-contents)


### tests and examples ###

add_executable(ctest EXCLUDE_FROM_ALL fortran/ctest.c)
target_link_libraries(ctest cgemmi)

add_executable(cpptest EXCLUDE_FROM_ALL tests/main.cpp tests/cif.cpp)

add_executable(hello EXCLUDE_FROM_ALL examples/hello.cc)
add_executable(doc_example EXCLUDE_FROM_ALL docs/doc_sym.cpp docs/doc_elem.cpp)
add_executable(doc_example2 EXCLUDE_FROM_ALL docs/doc_cif_cc.cpp)
add_executable(doc_maybegz EXCLUDE_FROM_ALL
               docs/doc_maybegz.cpp docs/doc_structure.cpp)
support_gz(doc_maybegz)

# auth_label requires <experimental/filesystem> and -lstdc++fs
add_executable(auth_label EXCLUDE_FROM_ALL examples/auth_label.cc)
target_link_libraries(auth_label stdc++fs)
support_gz(auth_label)

add_executable(ssbond EXCLUDE_FROM_ALL examples/ssbond.cc)
support_gz(ssbond)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION>)
add_test(NAME cpptest COMMAND cpptest)

add_dependencies(check
    ctest cpptest hello doc_example doc_example2 doc_maybegz ssbond)

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

