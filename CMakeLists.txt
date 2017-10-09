cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(gemmi C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build mode (copied from CMake FAQ)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

find_package(ZLIB REQUIRED)
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

add_library(input OBJECT src/input.cpp)

add_executable(gemmi-validate src/validate.cpp)
target_link_libraries(gemmi-validate ZLIB::ZLIB)

add_executable(gemmi-grep src/grep.cpp)
target_link_libraries(gemmi-grep ZLIB::ZLIB)

add_executable(gemmi-convert src/convert.cpp src/output.cpp
                             $<TARGET_OBJECTS:input>)
target_link_libraries(gemmi-convert ZLIB::ZLIB)

add_executable(gemmi-map src/map.cpp $<TARGET_OBJECTS:input>)
target_link_libraries(gemmi-map ZLIB::ZLIB)

add_executable(gemmi-mask src/mask.cpp $<TARGET_OBJECTS:input>)
target_link_libraries(gemmi-mask ZLIB::ZLIB)

add_executable(gemmi-sg src/sg.cpp)

if (benchmark_FOUND)
  foreach(b bench_elem bench_mod bench_sym)
    add_executable(${b} tools/${b}.cpp)
    target_link_libraries(${b} benchmark::benchmark)
  endforeach()
endif()

add_executable(hello examples/hello.cc)
add_executable(doc_example docs/doc_sym.cpp docs/doc_elem.cpp)
add_executable(doc_example2 docs/doc_cif_cc.cpp)
# auth_label requires <experimental/filesystem> and -lstdc++fs
#add_executable(examples/auth_label examples/auth_label.cc)
#target_link_libraries(examples/auth_label ZLIB::ZLIB)
