from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "coloring_utility",
        ["py_wrapper.cpp"],
        # MSVC specific flags: increase heap/stack if needed
        extra_compile_args=['/O2', '/std:c++14'],
    ),
]

setup(
    name="coloring_utility",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
