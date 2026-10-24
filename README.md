
# ColoringUtility

ColoringUtility is a C++-based graph coloring solver with a Python library for comfortable usage of the solver.

## Installation

It is recommended to use ColoringUtility as a Python library. However, it is also available as a pure C++ command
line application.

### Python Library

To install ColoringUtility Python library, make sure you have `setuptools` installed. It is recommended to use
virtualenv to install ColoringUtility.

The following can be simply used for installing the library:

```bash
python setup.py build_ext --inplace
```

### Building the C++ Solver

#### Windows

There are two prerequisites that need to be met before building ColoringUtility on Windows.

1. Make sure that you have MSYS2 installed - we will need it for `mingw64`.
2. Make sure that you have CMake installed.

Afterwards, the following script should be sufficient to build the package:

```bash
build_cpp.bat
```

If you installed MSYS2 to a different location than the default, then the `.bat` file might need to be tweaked.

## Structure

- `src` - C++ source code of the solver itself.
- `example_coloring.py` - An example of how to use ColoringUtility from Python.
- `example_world_map_coloring.py` - An example of how to use ColoringUtility to create a 4-colored world map.
- `data` - Testing data for the solver in the DIMACS `col` format.
- `tests` - A rich test suite for the solver using `pytest`, testing the solver on both generated graphs and static testing data.

## Usage

### Python Wrapper

To start with a simple example of using ColoringUtility, just run the following:

```python
python example_coloring.py
```

Refer to the contents of the script

### C++ Solver

The C++ solver provides an optimized implementation of several specialized algorithm to perform graph coloring (see below for references).

## Algorithms

ColoringUtility uses several graph coloring algorithms, including:

- Brelaz's heuristic [1]
- Iterated greedy coloring algorithm [2]
- TabuCol algorithm [3, 4]
- RLS_b algorithm [5]

It also implements a clique-based lower bound, a spectral lower bound, and a SAT-based lower bound based on solver togasat:

- https://github.com/togatoga/togasat

## References

1. Brélaz, D. (1979). New methods to color the vertices of a graph. Communications of the ACM, 22(4), 251-256.
2. Culberson, J. C. (1996). Exploring the k-colorable landscape with iterated greedy. DIMACS Ser. Discrete Math., 26, 245-284.
3. Hertz, A., & de Werra, D. (1987). Using tabu search techniques for graph coloring. Computing, 39(4), 345-351.
4. Galinier, P., & Hao, J. K. (1999). Hybrid evolutionary algorithms for graph coloring. Journal of combinatorial optimization, 3(4), 379-397.
5. Chalupa, D., & Nielsen, P. (2021). Parameter-free and cooperative local search algorithms for graph colouring. Soft Computing, 25(24), 15035-15050.


