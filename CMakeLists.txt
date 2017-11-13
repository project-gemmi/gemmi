cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(gemmi C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckIncludeFile)

# Set default build mode (copied from CMake FAQ)
if(NOT CMAKE_BUILD_TYPE)
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2")
  # avoid "unused" warnings when compiling with -DNDEBUG
  set(extra_release_flags "-Wno-unused-parameter -Wno-unused-variable")
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

if (benchmark_FOUND)
  foreach(b bench_elem bench_mod bench_sym)
    add_executable(${b} tools/${b}.cpp)
    target_link_libraries(${b} benchmark::benchmark)
  endforeach()
endif()

add_executable(hello examples/hello.cc)
add_executable(doc_example docs/doc_sym.cpp docs/doc_elem.cpp)
add_executable(doc_example2 docs/doc_cif_cc.cpp)
add_executable(doc_maybegz docs/doc_maybegz.cpp docs/doc_structure.cpp)
support_gz(doc_maybegz)
# auth_label requires <experimental/filesystem> and -lstdc++fs
#add_executable(examples/auth_label examples/auth_label.cc)
#support_gz(examples/auth_label)
