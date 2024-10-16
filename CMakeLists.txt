cmake_minimum_required(VERSION 3.30)
project(PyRecombine)

set(CMAKE_CXX_STANDARD 17)


set(Python_FIND_VIRTUALENV FIRST)
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GenerateExportHeader)


if (NOT CMAKE_SYSTEM_ARCHITECTURE MATCHES "x86")
    execute_process(COMMAND "${Python_EXECUTABLE}" "-c" [=[
import sysconfig
from pathlib import Path
p = Path(sysconfig.get_path('platlib'))
print(p.parent.parent)
    ]=]
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _out
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Adding ${_out}/cmake to search path")
    list(APPEND CMAKE_PREFIX_PATH "${_out}/cmake")
    find_package(MKL CONFIG REQUIRED)
endif()



Python_add_library(_PyRecombine MODULE WITH_SOABI
        src/py_recombine.c
        src/_recombine.cpp
        src/_recombine.h
)



# Distributing two shared libraries makes absolutely no sense here, so I'm going to manually
# compile recombine instead of using the build system that comes with it.

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm" OR (APPLE AND CMAKE_OSX_ARCHITECTURES MATCHES arm64))
    set(_no_mkl ON)
else()
    set(_no_mkl OFF)
endif()

if (WIN32)
    set(_bla_static OFF)
    set(_bla_sizeof_int 4)
    if (NOT _no_mkl)
        set(MKL_THREADING intel_threads)
        set(MKL_LINK static)
        set(MKL_INTERFACE ilp64)
        find_package(MKL CONFIG QUIET)
    endif()
elseif(APPLE)
    set(_bla_static ON)
    set(_bla_vendor Apple)
    set(_bla_sizeof_int 4) # Accelerate doesn't support 64-bit indexing
    if (NOT _no_mkl)
        # On Intel Mac, MKL is our best best, try to find that
        find_package(MKL CONFIG QUIET)
    endif()
else()
    set(_bla_sizeof_int 4)
    set(_bla_static ON)

    if (NOT _no_mkl)
        set(MKL_THREADING intel_threads)
        set(MKL_LINK static)
        set(MKL_INTERFACE ilp64)
        find_package(MKL CONFIG QUIET)
    endif()
endif()

if(NOT TARGET MKL::MKL)
    if (NOT DEFINED BLA_VENDOR)
        set(BLA_VENDOR ${_bla_vendor})
    endif()
    if (NOT DEFINED BLA_STATIC)
        set(BLA_STATIC ${_bla_static})
    endif()
    if (NOT DEFINED BLA_SIZEOF_INTEGER)
        set(BLA_SIZEOF_INTEGER ${_bla_sizeof_int})
    endif()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_link_libraries(_PyRecombine PRIVATE BLAS::BLAS LAPACK::LAPACK)
else()
    target_link_libraries(_PyRecombine PRIVATE MKL::MKL)
endif()

unset(_bla_vendor)
unset(_bla_static)
unset(_bla_sizeof_int)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

if (MSVC)
    set(OpenMP_CXX_FLAGS -openmp:experimental)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES INTERFACE_COMPILE_OPTIONS "-openmp:experimental")
    target_compile_options(_PyRecombine PRIVATE /fp:precise)
endif()


target_link_Libraries(_PyRecombine PRIVATE Python::NumPy OpenMP::OpenMP_CXX)

set(_rec_src ${CMAKE_CURRENT_LIST_DIR}/recombine/recombine)
target_sources(_PyRecombine PRIVATE
    ${_rec_src}/recombine/Compare.cpp
    ${_rec_src}/recombine/recombine.cpp
    ${_rec_src}/recombine/LinearAlgebraReductionTool.cpp
    ${_rec_src}/recombine/LinearAlgebraReductionTool.h
    ${_rec_src}/recombine/reweight.h
    ${_rec_src}/recombine/SafeInt3.hpp
    ${_rec_src}/recombine/TreeBufferHelper.h
    ${_rec_src}/recombine/TreeBufferHelper.cpp
    ${_rec_src}/recombine/RdToPowers.h
    ${_rec_src}/recombine/RdToPowers.cpp
)
target_include_directories(_PyRecombine PRIVATE "${_rec_src}" "${_rec_src}/recombine" "${CMAKE_CURRENT_BINARY_DIR}")

generate_export_header(_PyRecombine BASE_NAME RECOMBINE)


target_compile_definitions(_PyRecombine PRIVATE NOSIMPLEX REDUCTION_ALGO=svd NDEBUG)

SET_TARGET_PROPERTIES(_PyRecombine PROPERTIES
    LIBRARY_OUTPUT_NAME _recombine
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/src/pyrecombine
)



install(TARGETS
    _PyRecombine
    RUNTIME DESTINATION pyrecombine
    LIBRARY
    DESTINATION pyrecombine
    NAMELINK_SKIP
    ARCHIVE DESTINATION ${SKBUILD_NULL_DIR}
    COMPONENT Development
    EXCLUDE_FROM_ALL
    INCLUDES DESTINATION ${SKBUILD_NULL_DIR}
    COMPONENT Development
        FRAMEWORK DESTINATION roughpy
    EXCLUDE_FROM_ALL
)
install(DIRECTORY src/pyrecombine
    DESTINATION .
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "*.pyi")