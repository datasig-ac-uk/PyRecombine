
import os
import platform
import fnmatch
import sysconfig

from importlib import metadata as _ilm
from pathlib import Path as _Path



def _locate_file(module: str, name: str) -> _Path:
    """
    Locate a file within a module and return the correct absolute path to
    where the file is installed on the filesystem.

    This function raises an importlib.metadata.PackageNotFoundError upon
    failure.

    :param module: name of the module to locate a file within
    :param name: name of the file to find
    :return: Absolute path to the file as a pathlib.Path object
    """
    dist = _ilm.distribution(module)

    files = [dist.locate_file(path).resolve() for path in dist.files if fnmatch.fnmatch(path.name, name)]

    if not files:
        raise _ilm.PackageNotFoundError("Could not find specified file")

    file = files[0]
    assert file.exists()

    return file

if platform.system() == "Linux":
    from ctypes import CDLL as _CDLL

    try:
        _CDLL(str(_locate_file("intel_openmp", "libiomp5.so")))
        _CDLL(str(_locate_file("mkl", "libmkl_intel_ilp64.so.2")))
        _CDLL(str(_locate_file("mkl", "libmkl_core.so.2")))
        _CDLL(str(_locate_file("mkl", "libmkl_intel_thread.so.2")))


    except _ilm.PackageNotFoundError as e:
        raise ImportError("Could not find the MKL libraries") from e

elif platform.system() == "Windows":
    try:
        platdir_path = _locate_file("mkl", "libmkl_core.dll")
        os.add_dll_directory(str(platdir_path))
    except _ilm.PackageNotFoundError as e:
        raise ImportError("Could not find the MKL libraries") from e






from ._recombine import *