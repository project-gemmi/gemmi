cmake_minimum_required(VERSION 3.15...3.22)

# get version string from version.hpp
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/gemmi/version.hpp"
     gemmi_version_line REGEX "#define GEMMI_VERSION \"")
string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
       gemmi_version_str ${gemmi_version_line})

project(gemmi LANGUAGES C CXX VERSION ${gemmi_version_str})
message(STATUS "Gemmi version ${PROJECT_VERSION}")

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
  #set(CMAKE_INSTALL_BINDIR "${SKBUILD_SCRIPTS_DIR}")
endif()

option(BUILD_SHARED_LIBS "Build shared library" ON)
option(BUILD_GEMMI_PROGRAM "Build gemmi command-line program" ON)
option(INSTALL_DEV_FILES "Install C++ headers, cmake configs, static lib" ON)
option(USE_FORTRAN "Build Fortran bindings" OFF)
option(USE_PYTHON "Build Python bindings" OFF)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
option(USE_WMAIN "(Windows only) take Unicode arguments in gemmi program" ON)
set(GEMMI_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/gemmi" CACHE STRING
    "Install path for gemmi CMake files")
option(STRIP_BINARY "Strip symbols from program" OFF)

# uncomment to show compilation times for each compilation unit
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")

if (DEFINED ENV{FC} OR CMAKE_Fortran_COMPILER)
  set(USE_FORTRAN ON CACHE BOOL "Build Fortran bindings" FORCE)
endif()

if (NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()
if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif()

# We don't need to set policy CMP0069 for this file, because the minimum
# version is 3.15 and for versions >= 3.9 it's automatically set to NEW.
# But a FindPython script in cmake 3.15 have "cmake_policy(VERSION 3.7)"
# so we need to set CMP0069 because of that.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

# Unless it's explicitly specifed, we want LTO for non-debug builds only.
# pybind11Common.cmake has its own logic to set LTO options.
# Having different options in pybind11_add_module than for gemmi_cpp
# may cause problems -- we may need to handle it at some point.
if (NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT is_ipo_supported OUTPUT ipo_error)
  if(is_ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
  else()
    message(STATUS "Interprocedural optimization not supported: ${error}")
  endif()
endif()

if (INFO)
  set(GEMMI_VERSION_INFO ${INFO} CACHE STRING "Extra text for gemmi -V" FORCE)
endif()

if (USE_FORTRAN)
   message(STATUS "Currently, Fortran bindings are not built by cmake.")
   message(STATUS "They can be built by running make in fortran directory.")
  enable_language(Fortran)
else()
  #message(STATUS "Skipping Fortran bindings. Add -D USE_FORTRAN=1 to build them.")
endif()

if (DEFINED ENV{EXTRA_WFLAGS})
  set(EXTRA_WARNINGS ON CACHE BOOL "Set extra warning flags" FORCE)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)  # for check_include_file

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories("${ZLIB_INCLUDE_DIR}")
else()
  message(STATUS "The build will use zlib code from third_party/zlib.")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib")
endif()

if (NOT DEFINED SKBUILD AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(benchmark 1.3 QUIET)
endif()
if (benchmark_FOUND)
  message(STATUS "Found benchmark: ${benchmark_DIR}")
else (NOT benchmark_FOUND)
  message(STATUS "Benchmarks not configured.")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 -Wredundant-decls -Wfloat-conversion -Wdisabled-optimization -Wshadow $ENV{EXTRA_WFLAGS}")
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
  set_property(TARGET ungz PROPERTY POSITION_INDEPENDENT_CODE ON)
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

add_library(gemmi_headers INTERFACE)
add_library(gemmi::headers ALIAS gemmi_headers)
target_include_directories(gemmi_headers INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_compile_features(gemmi_headers INTERFACE cxx_std_11)
set_target_properties(gemmi_headers PROPERTIES EXPORT_NAME headers)

add_library(gemmi_cpp
            src/align.cpp src/assembly.cpp src/calculate.cpp src/crd.cpp
            src/ddl.cpp src/eig3.cpp
            src/mmcif.cpp src/mmread_gz.cpp src/mtz.cpp src/mtz2cif.cpp
            src/polyheur.cpp src/read_cif.cpp src/resinfo.cpp
            src/riding_h.cpp src/sprintf.cpp src/to_mmcif.cpp
            src/to_pdb.cpp src/monlib.cpp src/topo.cpp src/xds_ascii.cpp)
add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
set_property(TARGET gemmi_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET gemmi_cpp PROPERTY CXX_VISIBILITY_PRESET hidden)
target_compile_definitions(gemmi_cpp PRIVATE GEMMI_BUILD)
target_include_directories(gemmi_cpp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

if (BUILD_SHARED_LIBS)
  target_compile_definitions(gemmi_cpp PUBLIC GEMMI_SHARED)
endif()
target_link_libraries(gemmi_cpp PUBLIC gemmi_headers)
support_gz(gemmi_cpp)

### programs from prog/ ###

add_library(options OBJECT prog/options.cpp)
target_link_libraries(options PUBLIC gemmi_headers)
target_include_directories(options PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
if (GEMMI_VERSION_INFO)
  target_compile_definitions(options PRIVATE GEMMI_VERSION_INFO=${GEMMI_VERSION_INFO})
endif()

# Gemmi subcommands compiled as individual binaries.
# Not built by default and used mostly for testing.
function(add_single_prog target)
  cmake_parse_arguments(arg "NOLIB;GZ" "ADD_SRC" "" ${ARGN})
  add_executable(gemmi-${target} EXCLUDE_FROM_ALL prog/${target}.cpp ${arg_ADD_SRC}
                 $<TARGET_OBJECTS:options>)
  if (arg_NOLIB)
    target_link_libraries(gemmi-${target} PRIVATE gemmi_headers)
  else()
    target_link_libraries(gemmi-${target} PRIVATE gemmi_cpp)
  endif()
  target_include_directories(gemmi-${target} PRIVATE
                             "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  if (arg_GZ)
    support_gz(gemmi-${target})
  endif()
endfunction()

add_single_prog(align)
add_single_prog(blobs ADD_SRC prog/mapcoef.cpp GZ)
add_single_prog(cif2json)
add_single_prog(cif2mtz)
add_single_prog(contact)
add_single_prog(contents)
add_single_prog(convert)
add_single_prog(diff GZ)
add_single_prog(ecalc GZ)
add_single_prog(fprime NOLIB)
add_single_prog(grep NOLIB GZ)
add_single_prog(h)
add_single_prog(json2cif)
add_single_prog(map GZ)
add_single_prog(map2sf GZ)
add_single_prog(mask)
add_single_prog(mixmtz)
add_single_prog(merge GZ)
add_single_prog(mondiff)
add_single_prog(mtz GZ)
add_single_prog(mtz2cif GZ)
add_single_prog(prep ADD_SRC prog/monlib_opt.cpp GZ)
add_single_prog(reindex GZ)
add_single_prog(residues)
add_single_prog(rmsz)
add_single_prog(sf2map ADD_SRC prog/mapcoef.cpp GZ)
add_single_prog(sfcalc GZ)
add_single_prog(sg NOLIB)
add_single_prog(tags NOLIB GZ)
add_single_prog(validate ADD_SRC prog/validate_mon.cpp GZ)
add_single_prog(wcn)
add_single_prog(xds2mtz GZ)

# gemmi program with a number of subcommands
add_executable(gemmi_prog
               prog/align.cpp prog/blobs.cpp
               prog/cif2mtz.cpp prog/cif2json.cpp prog/contact.cpp
               prog/contents.cpp prog/convert.cpp
               prog/diff.cpp prog/ecalc.cpp prog/fprime.cpp
               prog/grep.cpp prog/h.cpp prog/json2cif.cpp
               prog/main.cpp prog/map.cpp prog/map2sf.cpp
               prog/mapcoef.cpp prog/mask.cpp prog/merge.cpp
               prog/mondiff.cpp prog/monlib_opt.cpp prog/mtz.cpp prog/mtz2cif.cpp
               prog/prep.cpp prog/reindex.cpp prog/residues.cpp prog/rmsz.cpp
               prog/sf2map.cpp prog/sfcalc.cpp prog/sg.cpp prog/tags.cpp
               prog/validate.cpp prog/validate_mon.cpp prog/wcn.cpp
               prog/xds2mtz.cpp
               $<TARGET_OBJECTS:options>)
add_executable(gemmi::prog ALIAS gemmi_prog)
target_link_libraries(gemmi_prog PRIVATE gemmi_cpp)
support_gz(gemmi_prog)
target_include_directories(gemmi_prog PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
target_compile_definitions(gemmi_prog PRIVATE GEMMI_ALL_IN_ONE=1)
set_target_properties(gemmi_prog PROPERTIES OUTPUT_NAME gemmi EXPORT_NAME prog)
if (WIN32 AND USE_WMAIN)
  # _UNICODE=1 is now set globally
  #target_compile_definitions(gemmi_prog PRIVATE _UNICODE=1)
  if(MINGW)
    target_link_options(gemmi_prog PRIVATE "-municode")
  endif()
endif()
if (NOT BUILD_GEMMI_PROGRAM)
  set_property(TARGET gemmi_prog PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# I'm not sure if this option should be here.
# The alternative is to use install/strip target, or cmake --install --strip.
if (STRIP_BINARY AND CMAKE_STRIP)
  add_custom_command(TARGET gemmi_prog POST_BUILD
                     COMMENT "Strip symbols from program with ${CMAKE_STRIP}"
                     COMMAND ${CMAKE_STRIP}
                     ARGS $<TARGET_FILE:gemmi_prog>)
endif()

### tests and examples ###

#add_executable(c_test EXCLUDE_FROM_ALL fortran/c_test.c)
#target_link_libraries(c_test PRIVATE cgemmi)

add_executable(cpptest EXCLUDE_FROM_ALL tests/main.cpp tests/cif.cpp
                                        src/mtz2cif.cpp src/eig3.cpp src/sprintf.cpp)
target_compile_definitions(cpptest PRIVATE USE_STD_SNPRINTF=1)
target_link_libraries(cpptest PRIVATE gemmi_headers)
target_include_directories(cpptest PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

add_executable(hello EXCLUDE_FROM_ALL examples/hello.cpp)
target_link_libraries(hello PRIVATE gemmi_headers)
add_executable(doc_example EXCLUDE_FROM_ALL
               docs/code/sym.cpp docs/code/elem.cpp docs/code/resinfo.cpp
               docs/code/cell.cpp src/resinfo.cpp)
target_link_libraries(doc_example PRIVATE gemmi_headers)
add_executable(doc_example2 EXCLUDE_FROM_ALL docs/code/cif_cc.cpp)
target_link_libraries(doc_example2 PRIVATE gemmi_headers)
add_executable(doc_maybegz EXCLUDE_FROM_ALL docs/code/maybegz.cpp docs/code/mutate.cpp)
target_link_libraries(doc_maybegz PRIVATE gemmi_cpp)
support_gz(doc_maybegz)
add_executable(doc_newmtz EXCLUDE_FROM_ALL docs/code/newmtz.cpp)
target_link_libraries(doc_newmtz PRIVATE gemmi_cpp)
support_gz(doc_newmtz)

# always compile these tests with assertions enabled
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
  target_compile_options(doc_example PRIVATE -UNDEBUG)
endif()

add_executable(test_disulf EXCLUDE_FROM_ALL tests/disulf.cpp)
target_link_libraries(test_disulf PRIVATE gemmi_cpp)
support_gz(test_disulf)

# auth_label requires <experimental/filesystem> and -lstdc++fs
add_executable(auth_label EXCLUDE_FROM_ALL examples/auth_label.cpp)
target_link_libraries(auth_label PRIVATE gemmi_headers)
if (NOT MSVC)
  target_link_libraries(auth_label PRIVATE stdc++fs)
endif()
support_gz(auth_label)

add_executable(check_conn EXCLUDE_FROM_ALL examples/check_conn.cpp)
target_link_libraries(check_conn PRIVATE gemmi_cpp)

enable_testing()

add_custom_target(print-version
  COMMAND gemmi_prog --version
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
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
  foreach(b stoi elem mod niggli pdb resinfo round sym writecif)
    add_executable(${b}-bm EXCLUDE_FROM_ALL benchmarks/${b}.cpp)
    if (b MATCHES "resinfo|pdb|writecif")
      target_link_libraries(${b}-bm PRIVATE gemmi_cpp)
    endif()
    target_link_libraries(${b}-bm PRIVATE gemmi_headers benchmark::benchmark)
    set_property(TARGET ${b}-bm PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                             "${CMAKE_BINARY_DIR}/benchmarks")
    add_dependencies(check ${b}-bm)
  endforeach()
endif()

### Python bindings ###
if (USE_PYTHON)
  message(STATUS "The python module will be built.")
  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  if (EXISTS "${CMAKE_HOME_DIRECTORY}/pybind11")
    message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/pybind11 (internal copy).")
    add_subdirectory(pybind11)
  else()
    # use pybind11-config (if available) to determine pybind11_DIR
    execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
    string(STRIP "${pybind11_DIR}" pybind11_DIR)
    find_package(pybind11 2.6 CONFIG REQUIRED)
    message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
  endif()
  pybind11_add_module(gemmi_py
          python/mol.cpp python/gemmi.cpp python/align.cpp
          python/ccp4.cpp python/chemcomp.cpp python/cif.cpp python/custom.cpp
          python/elem.cpp python/grid.cpp python/hkl.cpp
          python/meta.cpp python/monlib.cpp
          python/mtz.cpp python/read.cpp python/recgrid.cpp
          python/scaling.cpp python/search.cpp
          python/sf.cpp python/sym.cpp python/topo.cpp
          python/unitcell.cpp python/write.cpp
          $<TARGET_OBJECTS:gemmi_cpp>)
  target_link_libraries(gemmi_py PRIVATE gemmi_headers)
  set_property(TARGET gemmi_py PROPERTY OUTPUT_NAME gemmi)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(gemmi_py PRIVATE "-Wno-shadow")
  endif()
  support_gz(gemmi_py)
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()

# avoid GCC warning: the ABI of passing structure with 'complex float' member
# has changed in GCC 4.4
set_property(SOURCE python/recgrid.cpp python/hkl.cpp python/mtz.cpp prog/mtz.cpp
             PROPERTY COMPILE_FLAGS $<$<CXX_COMPILER_ID:GNU>:-Wno-psabi>)

include(CMakePackageConfigHelpers)
configure_package_config_file(tools/gemmi-config.cmake.in gemmi-config.cmake
                              INSTALL_DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
write_basic_package_version_file(gemmi-config-version.cmake
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

if (BUILD_GEMMI_PROGRAM)
  # Is it useful to have the program location in gemmi-targets.cmake?
  install(TARGETS gemmi_prog EXPORT GemmiTargets)
endif()

if (INSTALL_DEV_FILES)
  # In CMake 3.23 we'll be able to use FILE_SET to install the headers
  install(DIRECTORY include/gemmi DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(TARGETS gemmi_headers gemmi_cpp EXPORT GemmiTargets)
  install(EXPORT GemmiTargets FILE gemmi-targets.cmake NAMESPACE gemmi::
          DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
  install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/gemmi-config.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/gemmi-config-version.cmake"
          DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
elseif (BUILD_SHARED_LIBS AND BUILD_GEMMI_PROGRAM)
  install(TARGETS gemmi_cpp)
endif()

if (USE_PYTHON)
  if (DEFINED PYTHON_INSTALL_DIR)
    message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
    set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
  endif()
  file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
  # Using Python_SITEARCH/SITELIB is not good, because they are absolute
  # and don't respect CMAKE_INSTALL_PREFIX.
  # https://discourse.cmake.org/t/findpython3-how-to-specify-local-installation-directory-for-python-module/3580/5
  install(TARGETS gemmi_py DESTINATION "${Python_SITEARCH}")
endif()
