cmake_minimum_required(VERSION 3.16)

project(ColoringUtility VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Source Files ---
set(SOURCES
    src/main.cpp
    src/algorithm_greedyclique.cpp
    src/algorithm_igcol.cpp
    src/compute.cpp
    src/graphs.cpp
    src/common.cpp
    src/cli.cpp
    src/statistics.cpp
    src/random_generator.cpp
    src/tabu_base.cpp
    src/tabucol.cpp
    src/algorithm.cpp
    src/algorithm_brelaz.cpp
    src/h2col.cpp
)

set(HEADERS
    src/graphs.h
    src/algorithm_greedyclique.h
    src/algorithm_igcol.h
    src/common.h
    src/cli.h
    src/compute.h
    src/graphs_common.h
    src/statistics.h
    src/random_generator.h
    src/tabu_base.h
    src/tabucol.h
    src/algorithm.h
    src/algorithm_brelaz.h
    src/h2col.h
)

add_executable(ColoringUtility 
    ${SOURCES} 
    ${HEADERS} 
    ${FORMS}
)

if(MSVC)
    # MSVC equivalent to high optimization
    target_compile_options(ColoringUtility PRIVATE /O2)
else()
    target_compile_options(ColoringUtility PRIVATE -O3)
endif()

set(DIST_FILES
    build_py.bat
    coloring.py
    setup.py
    fractional_chrom_lb.py
    lovasz_theta_lb.py
)
source_group("DistFiles" FILES ${DIST_FILES})

set_target_properties(ColoringUtility PROPERTIES
    WIN32_EXECUTABLE FALSE
)
