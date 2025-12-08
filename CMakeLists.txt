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
option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" ON)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
option(USE_WMAIN "(Windows only) take Unicode arguments in gemmi program" ON)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to libgemmi_cpp DLL" ON)
if (WIN32)
  set(GEMMI_INSTALL_CMAKEDIR "cmake" CACHE STRING
      "Install path for gemmi CMake files")
else()
  set(GEMMI_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/gemmi" CACHE STRING
      "Install path for gemmi CMake files")
endif()
option(STRIP_BINARY "Strip symbols from program" OFF)
option(ADD_SOVERSION "For shared library (on Linux):\
 set ABI version to 0, build versions to project version" OFF)
option(USE_NLOPT "Experimental, for testing only, don't use it" OFF)

# uncomment to show compilation times for each compilation unit
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")

if (ONLY_PROGRAM)
  set(USE_PYTHON OFF CACHE BOOL "" FORCE) # it's set to ON in pyproject.toml
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  set(INSTALL_DEV_FILES OFF CACHE BOOL "" FORCE)
else()
  if ((DEFINED ENV{FC} AND NOT "$ENV{FC}" STREQUAL "") OR CMAKE_Fortran_COMPILER)
    set(USE_FORTRAN ON CACHE BOOL "Build Fortran bindings" FORCE)
  endif()
  if (DEFINED Python_EXECUTABLE AND NOT DEFINED SKBUILD)
    set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
  endif()
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
# TODO: check that it's compatible with nanobind_add_module, which has
# its own logic of setting LTO options.
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
  message(STATUS "They can be built by running make in fortran/.")
  #enable_language(Fortran)
else()
  #message(STATUS "Skipping Fortran bindings. Add -D USE_FORTRAN=1 to build them.")
endif()

if (DEFINED ENV{EXTRA_WFLAGS})
  set(EXTRA_WARNINGS ON CACHE BOOL "Set extra warning flags" FORCE)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

include(CheckIncludeFile)  # for check_include_file

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (FETCH_ZLIB_NG)
  message(STATUS " =============== fetching and building zlib-ng ===============")
  include(FetchContent)
  # FetchContent_Declare can take EXCLUDE_FROM_ALL only since CMake 3.28;
  # for now use FetchContent_Populate instead of FetchContent_MakeAvailable.
  if (POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
  endif()
  FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
                               GIT_TAG        2.2.x
                               GIT_SHALLOW TRUE)
  FetchContent_GetProperties(zlib-ng)
  if (NOT zlib-ng_POPULATED)
    FetchContent_Populate(zlib-ng)
    set(ZLIB_ENABLE_TESTS OFF CACHE INTERNAL "")
    set(BUILD_SHARED_LIBS_here ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "build static zlib-ng")
    add_subdirectory(${zlib-ng_SOURCE_DIR} ${zlib-ng_BINARY_DIR} EXCLUDE_FROM_ALL)
    set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_here} CACHE BOOL "Build shared library" FORCE)
    set_property(TARGET zlib PROPERTY POSITION_INDEPENDENT_CODE ON)
    # zlib-ng, unless built WITH_NATIVE_INSTRUCTIONS=ON (-march=native),
    # builds multiple variants of some functions, using different CPU features.
    # IPO (LTO) has been reported to mess it up leading to illegal-instruction
    # crash. When zlib-ng is built as a standalone project, not subdirectory,
    # these flags are set in zlib-ng/CMakeLists.txt.
    if (NOT WITH_NATIVE_INSTRUCTIONS)
      set_target_properties(zlib PROPERTIES INTERPROCEDURAL_OPTIMIZATION OFF
                                            INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF
                                            INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL OFF)
    endif()
    message(STATUS " =============== end of zlib-ng configuration ===============\n")
  endif()
elseif (USE_ZLIB_NG)
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools")
  find_package(ZLIB_NG REQUIRED)
else()
  if (NOT INTERNAL_ZLIB)
    # CMake >= 3.24 respects variable ZLIB_USE_STATIC_LIBS
    find_package(ZLIB)
  endif()
  if (NOT ZLIB_FOUND)
    message(STATUS "The build will use zlib code from third_party/zlib.")
  endif()
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
target_compile_features(gemmi_headers INTERFACE cxx_std_14)
set_target_properties(gemmi_headers PROPERTIES EXPORT_NAME headers)

add_library(gemmi_cpp
            src/align.cpp src/assembly.cpp src/calculate.cpp src/ccp4.cpp
            src/crd.cpp src/ddl.cpp src/eig3.cpp src/flat.cpp src/fprime.cpp src/gz.cpp
            src/intensit.cpp src/json.cpp src/mmcif.cpp src/mmread_gz.cpp
            src/monlib.cpp src/mtz.cpp src/mtz2cif.cpp
            src/pdb.cpp src/polyheur.cpp src/read_cif.cpp
            src/resinfo.cpp src/riding_h.cpp
            src/select.cpp src/sprintf.cpp src/dssp.cpp src/symmetry.cpp
            src/to_json.cpp src/to_mmcif.cpp src/to_pdb.cpp src/topo.cpp
            src/xds_ascii.cpp)
add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
set_property(TARGET gemmi_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET gemmi_cpp PROPERTY CXX_VISIBILITY_PRESET hidden)
target_compile_definitions(gemmi_cpp PRIVATE GEMMI_BUILD)
target_include_directories(gemmi_cpp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

if (BUILD_SHARED_LIBS)
  target_compile_definitions(gemmi_cpp PUBLIC GEMMI_SHARED)
  if (ADD_SOVERSION)
    set_target_properties(gemmi_cpp PROPERTIES VERSION ${PROJECT_VERSION}
                                               SOVERSION "0")
  endif()
endif()
target_link_libraries(gemmi_cpp PUBLIC gemmi_headers)

if (FETCH_ZLIB_NG)
  target_compile_definitions(gemmi_cpp PRIVATE USE_ZLIB_NG=1)
  target_link_libraries(gemmi_cpp PRIVATE zlib)
elseif (ZLIB_NG_FOUND)
  target_compile_definitions(gemmi_cpp PRIVATE USE_ZLIB_NG=1)
  target_link_libraries(gemmi_cpp PRIVATE ZLIB_NG::ZLIB_NG)
elseif (ZLIB_FOUND)
  target_link_libraries(gemmi_cpp PRIVATE ZLIB::ZLIB)
else()
  target_sources(gemmi_cpp PRIVATE
      "third_party/zlib/adler32.c"
      "third_party/zlib/crc32.c"
      "third_party/zlib/gzlib.c"
      "third_party/zlib/gzread.c"
      "third_party/zlib/inflate.c"
      "third_party/zlib/inftrees.c"
      "third_party/zlib/inffast.c"
      "third_party/zlib/zutil.c")
  check_include_file(unistd.h has_unistd_h)
  target_compile_definitions(gemmi_cpp PRIVATE NO_GZCOMPRESS=1 DYNAMIC_CRC_TABLE=1)
  if (has_unistd_h)
    target_compile_definitions(gemmi_cpp PRIVATE Z_HAVE_UNISTD_H=1)
  endif()
  target_include_directories(gemmi_cpp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib")
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
    target_compile_options(gemmi_cpp PRIVATE /wd4267)
  endif()
endif()


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
  cmake_parse_arguments(arg "" "ADD_SRC" "" ${ARGN})
  add_executable(gemmi-${target} EXCLUDE_FROM_ALL prog/${target}.cpp ${arg_ADD_SRC}
                 $<TARGET_OBJECTS:options>)
  target_link_libraries(gemmi-${target} PRIVATE gemmi_cpp)
  target_include_directories(gemmi-${target} PRIVATE
                             "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
endfunction()

add_single_prog(align)
add_single_prog(blobs ADD_SRC prog/mapcoef.cpp)
add_single_prog(cif2json)
add_single_prog(cif2mtz)
add_single_prog(cifdiff)
add_single_prog(contact)
add_single_prog(contents)
add_single_prog(convert)
add_single_prog(crd ADD_SRC prog/monlib_opt.cpp)
add_single_prog(ecalc)
add_single_prog(fprime)
add_single_prog(grep)
add_single_prog(h)
add_single_prog(json2cif)
add_single_prog(map)
add_single_prog(map2sf)
add_single_prog(mask)
add_single_prog(mixmtz)
add_single_prog(merge)
add_single_prog(mondiff)
add_single_prog(mtz)
add_single_prog(mtz2cif)
add_single_prog(reindex)
add_single_prog(residues)
add_single_prog(rmsz)
add_single_prog(sf2map ADD_SRC prog/mapcoef.cpp)
add_single_prog(sfcalc)
add_single_prog(ss)
add_single_prog(sg)
add_single_prog(tags)
add_single_prog(validate ADD_SRC prog/validate_mon.cpp)
add_single_prog(wcn)
add_single_prog(xds2mtz)

# gemmi program with a number of subcommands
add_executable(gemmi_prog
               prog/align.cpp prog/blobs.cpp
               prog/cif2mtz.cpp prog/cif2json.cpp prog/cifdiff.cpp
               prog/contact.cpp prog/contents.cpp prog/convert.cpp prog/crd.cpp
               prog/ecalc.cpp prog/fprime.cpp
               prog/grep.cpp prog/h.cpp prog/json2cif.cpp
               prog/main.cpp prog/map.cpp prog/map2sf.cpp
               prog/mapcoef.cpp prog/mask.cpp prog/merge.cpp
               prog/mondiff.cpp prog/monlib_opt.cpp prog/mtz.cpp prog/mtz2cif.cpp
               prog/reindex.cpp prog/residues.cpp prog/rmsz.cpp prog/set.cpp
               prog/sf2map.cpp prog/sfcalc.cpp prog/sg.cpp prog/ss.cpp
               prog/tags.cpp prog/validate.cpp prog/validate_mon.cpp
               prog/wcn.cpp prog/xds2mtz.cpp
               $<TARGET_OBJECTS:options>)
add_executable(gemmi::prog ALIAS gemmi_prog)
target_link_libraries(gemmi_prog PRIVATE gemmi_cpp)
target_include_directories(gemmi_prog PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
target_compile_definitions(gemmi_prog PRIVATE GEMMI_ALL_IN_ONE=1)
set_target_properties(gemmi_prog PROPERTIES OUTPUT_NAME gemmi EXPORT_NAME prog)
if (BUILD_SHARED_LIBS)
  if (APPLE)
    set(origin "@loader_path")
  else ()
    set(origin "$ORIGIN")
  endif ()
  set_property(TARGET gemmi_prog
               PROPERTY INSTALL_RPATH "${origin}/../${CMAKE_INSTALL_LIBDIR}")
endif()
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

add_executable(cpptest EXCLUDE_FROM_ALL tests/main.cpp tests/cif.cpp tests/windowsh.cpp)
target_compile_definitions(cpptest PRIVATE USE_STD_SNPRINTF=1)
target_link_libraries(cpptest PRIVATE gemmi_cpp)
target_include_directories(cpptest PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

add_executable(doc_example EXCLUDE_FROM_ALL
               docs/code/sym.cpp docs/code/elem.cpp docs/code/resinfo.cpp
               docs/code/cell.cpp docs/code/mutate.cpp)
target_link_libraries(doc_example PRIVATE gemmi_cpp)
add_executable(doc_example2 EXCLUDE_FROM_ALL docs/code/cif_cc.cpp)
target_link_libraries(doc_example2 PRIVATE gemmi_headers)
add_executable(doc_newmtz EXCLUDE_FROM_ALL docs/code/newmtz.cpp)
target_link_libraries(doc_newmtz PRIVATE gemmi_cpp)

# always compile these tests with assertions enabled
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|Intel")
  target_compile_options(doc_example PRIVATE -UNDEBUG)
endif()

add_executable(test_disulf EXCLUDE_FROM_ALL tests/disulf.cpp)
target_link_libraries(test_disulf PRIVATE gemmi_cpp)

# auth_label requires C++17 for <filesystem>
add_executable(auth_label EXCLUDE_FROM_ALL examples/auth_label.cpp)
target_link_libraries(auth_label PRIVATE gemmi_cpp)

add_executable(check_conn EXCLUDE_FROM_ALL examples/check_conn.cpp)
target_link_libraries(check_conn PRIVATE gemmi_cpp)

enable_testing()

add_custom_target(print-version
  COMMAND gemmi_prog --version
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  COMMENT "gemmi --version"
)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -C \"$<CONFIG>\")
add_test(NAME cpptest COMMAND cpptest)

add_dependencies(check
    cpptest doc_example doc_example2 doc_newmtz
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
    if (b MATCHES "elem|resinfo|pdb|sym|writecif")
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

  # workaround for https://github.com/wjakob/nanobind/issues/1237
  set(NB_USE_SUBMODULE_DEPS ON)

  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  if (EXISTS "${CMAKE_HOME_DIRECTORY}/nanobind")
    message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/nanobind (internal copy).")
    add_subdirectory(nanobind)
  else()
    # Detect the installed nanobind package and import it into CMake
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nb_cmake_dir)
    list(APPEND CMAKE_PREFIX_PATH "${nb_cmake_dir}")
    find_package(nanobind 2.4.0 CONFIG REQUIRED)
    message(STATUS "Found nanobind ${nanobind_VERSION}: ${NB_DIR}")
  endif()
  nanobind_add_module(gemmi_py NOMINSIZE
          python/gemmi.cpp python/align.cpp
          python/ccp4.cpp python/chemcomp.cpp python/cif.cpp
          python/elem.cpp python/grid.cpp python/hkl.cpp
          python/meta.cpp python/mol.cpp python/monlib.cpp
          python/mtz.cpp python/read.cpp python/recgrid.cpp
          python/scaling.cpp python/search.cpp
          python/sf.cpp python/small.cpp python/sym.cpp python/topo.cpp
          python/unitcell.cpp python/write.cpp)

  if (STANDALONE_PYTHON_MODULE)
    target_sources(gemmi_py PRIVATE $<TARGET_OBJECTS:gemmi_cpp>)
    get_target_property(_gemmi_cpp_libs gemmi_cpp LINK_LIBRARIES)
    set_property(TARGET gemmi_py APPEND PROPERTY LINK_LIBRARIES ${_gemmi_cpp_libs})
  else()
    target_link_libraries(gemmi_py PRIVATE gemmi_cpp)
  endif()
  set_property(TARGET gemmi_py PROPERTY OUTPUT_NAME gemmi_ext)
  if (CMAKE_CONFIGURATION_TYPES)
    set(py_dir $<CONFIG>/py)
  else()
    set(py_dir py)
  endif()
  set_property(TARGET gemmi_py PROPERTY LIBRARY_OUTPUT_DIRECTORY ${py_dir}/gemmi)
  add_custom_command(TARGET gemmi_py PRE_BUILD
    COMMENT "Copying gemmi/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/python/gemmi/__init__.py"
            ${py_dir}/gemmi/__init__.py
    BYPRODUCTS ${py_dir}/gemmi/__init__.py)
  add_custom_command(TARGET gemmi_py POST_BUILD
    COMMENT "Copying gemmi/fetch.py"
    COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/python/gemmi/fetch.py"
            ${py_dir}/gemmi/fetch.py
    BYPRODUCTS ${py_dir}/gemmi/fetch.py)

  # nanobind gives warnings with -Wpedantic and -Wshadow
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    if (TARGET nanobind-static)
      target_compile_options(nanobind-static PRIVATE "-Wno-shadow")
    endif()
    target_compile_options(gemmi_py PRIVATE "-Wno-shadow")
  endif()
  if (CMAKE_CXX_FLAGS MATCHES "-Wpedantic")
    if (TARGET nanobind-static)
      target_compile_options(nanobind-static PRIVATE "-Wno-pedantic")
    endif()
    target_compile_options(gemmi_py PRIVATE "-Wno-pedantic")
  endif()
  if (CMAKE_CXX_FLAGS MATCHES "-Wredundant-decls")
    if (TARGET nanobind-static)
      target_compile_options(nanobind-static PRIVATE "-Wno-redundant-decls")
    endif()
  endif()

  if (USE_NLOPT)
    message(STATUS "Using NLopt.")
    find_package(NLopt REQUIRED)
    target_link_libraries(gemmi_py PRIVATE NLopt::nlopt)
    target_compile_definitions(gemmi_py PRIVATE WITH_NLOPT=1)
  endif()

  if (GENERATE_STUBS AND CMAKE_CROSSCOMPILING)
    message(WARNING "Stubs cannot be generated when cross-compiling - skipping.")
  elseif (GENERATE_STUBS)
    nanobind_add_stub(
      gemmi_stub
      MODULE gemmi
      OUTPUT ${py_dir}/gemmi/__init__.pyi
      PYTHON_PATH ${py_dir}
      DEPENDS gemmi_py
      PATTERN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/python/stubs.pat"
      MARKER_FILE ${py_dir}/gemmi/py.typed
    )
    nanobind_add_stub(
      gemmi_cif_stub
      MODULE gemmi.cif
      OUTPUT ${py_dir}/gemmi/cif.pyi
      PYTHON_PATH ${py_dir}
      DEPENDS gemmi_py
    )
  endif()
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
  if (BUILD_SHARED_LIBS AND NOT STANDALONE_PYTHON_MODULE)
    if (DEFINED SKBUILD)
      # When building via pip, Python_SITEARCH is /tmp.../wheel/platlib,
      # different from the actual install path. Best guess:
      set(rpath_rel "../../../../${CMAKE_INSTALL_LIBDIR}")
    else()
      file(RELATIVE_PATH rpath_rel
           "${Python_SITEARCH}/gemmi" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    endif()
    set_target_properties(gemmi_py PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
  endif()
  install(TARGETS gemmi_py DESTINATION "${Python_SITEARCH}/gemmi" COMPONENT py)
  install(DIRECTORY examples DESTINATION "${Python_SITEARCH}/gemmi" COMPONENT py
          FILES_MATCHING PATTERN "*.py"
          PATTERN "[._]*" EXCLUDE)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/gemmi/__init__.py"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/gemmi/fetch.py"
          DESTINATION "${Python_SITEARCH}/gemmi" COMPONENT py)
  if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILING)
    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/gemmi/py.typed"
            "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/gemmi/__init__.pyi"
            "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/gemmi/cif.pyi"
            DESTINATION "${Python_SITEARCH}/gemmi" COMPONENT py)
  endif()
endif()
