cmake_minimum_required(VERSION 3.30)
project(PyRecombine)

set(CMAKE_CXX_STANDARD 17)


set(Python_FIND_VIRTUALENV FIRST)
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)


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

add_subdirectory(recombine)


Python_add_library(_PyRecombine MODULE WITH_SOABI
        src/py_recombine.c
        src/_recombine.cpp
        src/_recombine.h
)


target_link_Libraries(_PyRecombine PRIVATE Python::NumPy recombine::recombine)


SET_TARGET_PROPERTIES(_PyRecombine PROPERTIES
    LIBRARY_OUTPUT_NAME _recombine
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/src/pyrecombine
)