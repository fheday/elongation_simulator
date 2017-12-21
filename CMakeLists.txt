cmake_minimum_required(VERSION 2.6)
project(gillespie)

#-fopt-info-vec-missed -g -fno-inline -flto
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native -ftree-vectorize -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -Wall -Wextra -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

set(common_sources concentrationsreader.cpp mrna_reader.cpp)

set(ribosome_simulator_sources ribosomesimulator.cpp)

set(enlongation_simulator_sources enlogationsimulator.cpp  enlongation_codon.cpp initiationterminationcodon.cpp mrnaelement.cpp translation.cpp ribosomesimulator.cpp)

#add_executable(gillespie main.cpp ${common_sources}) # Gillespie implementation more tighly tied to ribosome/mRNA simulator. disabling the canonical one.
add_executable(calculateCodonsTimes calculateCodonsTimes.cpp ${common_sources} ${ribosome_simulator_sources})
add_executable(run_translation run_translation.cpp ${common_sources} ${enlongation_simulator_sources})


set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/experiments)
##python interface

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)
add_library(translation_py EXCLUDE_FROM_ALL MODULE translation.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp ribosomesimulator.cpp mrnaelement.cpp enlongation_codon.cpp)
target_compile_definitions(translation_py PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(translation_py PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(translation_py PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")


add_library(ribosomesimulator_py EXCLUDE_FROM_ALL MODULE concentrationsreader.cpp ribosomesimulator.cpp)
target_compile_definitions(ribosomesimulator_py PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(ribosomesimulator_py PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(ribosomesimulator_py PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")
