cmake_minimum_required(VERSION 3.29)
project(PyRecombine)

set(CMAKE_CXX_STANDARD 17)


set(Python_FIND_VIRTUALENV FIRST)
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GenerateExportHeader)

# This function is borrowed from the RoughPy CMake helper scripts
function(get_brew_prefix _out_var _package)

    cmake_parse_arguments(
            ""
            "VERBOSE"
            "BREW_EXECUTABLE"
            ""
            ${ARGN}
    )

    set(_executable brew)
    if (_BREW_EXECUTABLE)
        set(_executable "${_BREW_EXECUTABLE}")
    elseif (DEFINED ROUGHPY_HOMEBREW_EXECUTABLE)
        set(_executable "${ROUGHPY_HOMEBREW_EXECUTABLE}")
    elseif (DEFINED ENV{ROUGHPY_HOMEBREW_EXECUTABLE})
        set(_executable "$ENV{ROUGHPY_HOMEBREW_EXECUTABLE}")
    endif ()

    if (NOT _VERBOSE AND ENV{VERBOSE})
        set(_VERBOSE ON)
    endif ()

    set(_verbosity_flag "")
    if (_VERBOSE)
        set(_verbosity_flag "ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE")

        execute_process(COMMAND ${_executable} info ${_package}
                RESULT_VARIABLE _result
                ERROR_VARIABLE _err
        )
        message(STATUS "${_package} info:\n${_result}\n${_err}")
        unset(_result)
        unset(_err)
    endif ()
    message(DEBUG "Locating ${_package} brew installation prefix")

    execute_process(COMMAND ${_executable}
            "--prefix" "${_package}"
            RESULT_VARIABLE _result
            OUTPUT_VARIABLE _out
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ${_verbosity_flat}
    )

    if (_result EQUAL 0)
        message(STATUS "Brew located ${_package} at ${_out}")
        set(${_out_var} "${_out}" PARENT_SCOPE)
    else ()
        message(DEBUG "Could not locate ${_package} using Brew")
        set(${_out_var} "${_out_var}-NOTFOUND")
    endif ()
endfunction()

if (APPLE)
    # Apple Clang disables OpenMP for some unknown reason, so we have
    # to make sure CMake can find the headers/library files so it can
    # be turned back on. Since this is usually Brew installed, this
    # is a sensible way to find it.
    get_brew_prefix(_omp_prefix libomp)
    if (_omp_prefix)
        message(STATUS "Setting OpenMP search dir to ${_omp_prefix}")
        list(PREPEND CMAKE_PREFIX_PATH "${_omp_prefix}")

        set(_omp_lib "${_omp_prefix}/lib/libomp.dylib")
        if (NOT EXISTS "${_omp_Lib}")
            execute_process(COMMAND brew reinstall libomp COMMAND_ECHO STDOUT)
            if (NOT EXISTS "${_omp_lib}")
                message(FATAL_ERROR "WTH Apple")
            endif()
        endif()

        link_directories(${_omp_prefix}/lib)
        include_directories(${_omp_prefix}/include)
    endif ()
endif ()

message(STATUS "Target architecture ${CMAKE_SYSTEM_PROCESSOR}")
if (NOT APPLE AND (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64"))
    execute_process(COMMAND "${Python_EXECUTABLE}" "-c" [=[
from importlib import metadata as ilm
dist = ilm.distribution("mkl-devel")
files = [dist.locate_file(path).resolve() for path in dist.files if path.name == "MKLConfig.cmake"]
assert len(files) > 0
print(files[0].parent)
    ]=]
            RESULT_VARIABLE _result
            OUTPUT_VARIABLE _out
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Adding ${_out} to search path")
    list(APPEND CMAKE_PREFIX_PATH "${_out}")

    set(MKL_ARCH intel64)
    set(MKL_THREADING intel_thread)
    if (WIN32)
        set(MKL_THREADING sequential)
    endif()
    set(MKL_INTERFACE ilp64)
    # If mkl-static is installed, we can use this.
    # It isn't clear that it is best to do so though.
     set(MKL_LINK static)
#    set(MKL_LINK dynamic)
    find_package(MKL CONFIG REQUIRED)
endif ()


Python_add_library(_PyRecombine MODULE WITH_SOABI
        src/py_recombine.c
        src/_recombine.cpp
        src/_recombine.h
        src/recombine_common.h
        src/aligned_vec.h
        src/recombine_internal.h
        src/Recombine.cpp
        src/TreeBufferHelper.h
        src/LinearAlgebraReductionTool.cpp
        src/LinearAlgebraReductionTool.h
        src/SafeInt3.hpp
        src/reweight.cpp
        src/reweight.h
        src/lapack_definitions.h
        src/RdToPowers.cpp
        src/RdToPowers.h
        src/checked_assign.h
        src/aligned_alloc.cpp
)


# Distributing two shared libraries makes absolutely no sense here, so I'm going to manually
# compile recombine instead of using the build system that comes with it.

if (NOT TARGET MKL::MKL)
    if (WIN32)
        set(_bla_static OFF)
        set(_bla_sizeof_int 8)
    elseif (APPLE)
        set(_bla_static ON)
        set(_bla_vendor Apple)
        set(_bla_sizeof_int 4) # Accelerate doesn't support 64-bit indexing
    else ()
        set(_bla_sizeof_int 8)
        set(_bla_static ON)
    endif ()

    if (NOT DEFINED BLA_VENDOR)
        set(BLA_VENDOR ${_bla_vendor})
    endif ()
    if (NOT DEFINED BLA_STATIC)
        set(BLA_STATIC ${_bla_static})
    endif ()
    if (NOT DEFINED BLA_SIZEOF_INTEGER)
        set(BLA_SIZEOF_INTEGER ${_bla_sizeof_int})
    endif ()

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    unset(_bla_vendor)
    unset(_bla_static)
    unset(_bla_sizeof_int)

    target_link_libraries(_PyRecombine PRIVATE BLAS::BLAS LAPACK::LAPACK)
    target_compile_definitions(_PyRecombine PRIVATE RECOMBINE_INT_SIZE=${BLA_SIZEOF_INTEGER})
else ()
    if (MKL_INTERFACE STREQUAL ilp64)
        set(_mkl_int 8)
    else()
        set(_mkl_int 4)
    endif()
    target_compile_definitions(_PyRecombine PRIVATE RECOMBINE_INT_SIZE=${_mkl_int})
    target_link_libraries(_PyRecombine PRIVATE MKL::MKL)
    if (MSVC)
#        target_link_options(_PyRecombine PRIVATE "/nodefaultlib:vcomp")
    endif()

    unset(_mkl_int)
endif ()


find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

if (MSVC)
    set(OpenMP_CXX_FLAGS -openmp:experimental)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES INTERFACE_COMPILE_OPTIONS "-openmp:experimental")
    target_compile_options(_PyRecombine PRIVATE /fp:precise)
endif ()


target_link_Libraries(_PyRecombine PRIVATE Python::NumPy OpenMP::OpenMP_CXX Threads::Threads)

target_include_directories(_PyRecombine PRIVATE "${_rec_src}" "${_rec_src}/recombine" "${CMAKE_CURRENT_BINARY_DIR}")

generate_export_header(_PyRecombine BASE_NAME RECOMBINE)


target_compile_definitions(_PyRecombine PRIVATE NOSIMPLEX REDUCTION_ALGO=svd NDEBUG)

SET_TARGET_PROPERTIES(_PyRecombine PROPERTIES
        LIBRARY_OUTPUT_NAME _recombine
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
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
        PATTERN "*.pyi"
)




add_custom_target(pytest
        COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_LIST_DIR}/src" "${Python_EXECUTABLE}" "-m" "pytest" "${CMAKE_CURRENT_LIST_DIR}/tests"
        SOURCES tests/test_recombine.py
        USES_TERMINAL

)
add_dependencies(pytest _PyRecombine)


enable_testing()
add_test(NAME PyTest
        COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_LIST_DIR}/src" "${Python_EXECUTABLE}" "-m" "pytest" "${CMAKE_CURRENT_LIST_DIR}/tests")