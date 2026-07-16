import sys
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "coloring_utility",
        [
            "src/py_wrapper.cpp",
            "src/compute.cpp",
            "src/algorithm.cpp",
            "src/algorithm_brelaz.cpp",
            "src/algorithm_greedyclique.cpp",
            "src/algorithm_igcol.cpp",
            "src/random_generator.cpp",
            "src/tabu_base.cpp",
            "src/tabucol.cpp",
        ],
        include_dirs=["3rdparty"],        # "togasat.hpp" should resolve
        cxx_std=14,                       # let pybind11 pick the right per-compiler flag
    ),
]

# MSVC and GCC/Clang need different flags for the same intent
# (optimization and pthreads on Unix - due to std::thread)
extra_compile_args = {
    "msvc": ["/O2"],
    "unix": ["-O3", "-pthread"],
}
extra_link_args = {
    "msvc": [],
    "unix": ["-pthread"],
}

class custom_build_ext(build_ext):
    def build_extensions(self):
        compiler_type = self.compiler.compiler_type
        # cxx_std=14 already set /std:c++14 or -std=c++14 via pybind11;
        # here we just layer on optimization + threading flags.
        ctype = "msvc" if compiler_type == "msvc" else "unix"
        for ext in self.extensions:
            ext.extra_compile_args = extra_compile_args[ctype]
            ext.extra_link_args = extra_link_args[ctype]
        super().build_extensions()

setup(
    name="coloring_utility",
    ext_modules=ext_modules,
    cmdclass={"build_ext": custom_build_ext},
)
