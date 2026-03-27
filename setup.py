from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "coloring_utility",
        [
            "py_wrapper.cpp",
            "compute.cpp",
            "algorithm.cpp",
            "algorithm_brelaz.cpp",
            "algorithm_greedyclique.cpp",
            "algorithm_igcol.cpp",
            "random_generator.cpp",
            "tabu_base.cpp",
            "tabucol.cpp",
        ],
        # MSVC specific flags: increase heap/stack if needed
        extra_compile_args=['/O3', '/std:c++14'],
    ),
]

setup(
    name="coloring_utility",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
