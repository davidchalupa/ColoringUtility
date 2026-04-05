cmake_minimum_required(VERSION 3.16)

project(ColoringUtility VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Source Files ---
set(SOURCES
    main.cpp
    algorithm_greedyclique.cpp
    algorithm_igcol.cpp
    compute.cpp
    graphs.cpp
    common.cpp
    cli.cpp
    statistics.cpp
    random_generator.cpp
    tabu_base.cpp
    tabucol.cpp
    algorithm.cpp
    algorithm_brelaz.cpp
    h2col.cpp
)

set(HEADERS
    graphs.h
    algorithm_greedyclique.h
    algorithm_igcol.h
    common.h
    cli.h
    compute.h
    graphs_common.h
    statistics.h
    random_generator.h
    tabu_base.h
    tabucol.h
    algorithm.h
    algorithm_brelaz.h
    h2col.h
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
