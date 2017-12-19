cmake_minimum_required(VERSION 2.6)
project(gillespie)

#-fopt-info-vec-missed -g -fno-inline
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -flto -march=native -ftree-vectorize -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -g3 -Wall -Wextra -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined")

set(common_sources concentrationsreader.cpp enlogationsimulator.cpp  mrna_reader.cpp ribosomesimulator.cpp ribosomesimulator.cpp enlongation_codon.cpp initiationterminationcodon.cpp mrnaelement.cpp translation.h translation.cpp)

#add_executable(gillespie main.cpp ${common_sources}) # Gillespie implementation more tighly tied to ribosome/mRNA simulator. disabling the canonical one.
add_executable(calculateCodonsTimes calculateCodonsTimes.cpp ${common_sources})
add_executable(run_translation run_translation.cpp ${common_sources})


##python interface
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)
add_library(translation MODULE translation.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp ribosomesimulator.cpp mrnaelement.cpp
enlongation_codon.cpp)
target_compile_definitions(translation PRIVATE -DCOMIPLE_PYTHON_MODULE)
target_link_libraries(translation PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(translation PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")
