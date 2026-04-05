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
        # MSVC specific flags: increase heap/stack if needed
        extra_compile_args=['/O3', '/std:c++14'],
    ),
]

setup(
    name="coloring_utility",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
