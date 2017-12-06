cmake_minimum_required(VERSION 2.6)
project(gillespie)

#-fopt-info-vec-missed -g -fno-inline
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -flto -march=native -ftree-vectorize -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -g3 -Wall -Wextra -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined")

set(common_sources concentrationsreader.cpp concentrationsreader.h enlogationsimulator.cpp enlogationsimulator.h gillespie.cpp gillespie.h mrna_reader.cpp mrna_reader.h ratecalculator.cpp ratecalculator.h reactionsset.cpp reactionsset.h ribosomesimulator.cpp ribosomesimulator.h ribosomesimulator.cpp enlongation_codon.h enlongation_codon.cpp initiationterminationcodon.h  initiationterminationcodon.cpp mrnaelement.h  mrnaelement.cpp translation.h translation.cpp)

add_executable(gillespie main.cpp ${common_sources})
add_executable(calculateCodonsTimes calculateCodonsTimes.cpp ${common_sources})
add_executable(run_translation run_translation.cpp ${common_sources})
install(TARGETS gillespie RUNTIME DESTINATION bin)


##python interface
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)
add_library(translation MODULE translation.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp ribosomesimulator.cpp mrnaelement.cpp gillespie.cpp
enlongation_codon.cpp reactionsset.cpp ratecalculator.cpp)
target_compile_definitions(translation PRIVATE -DCOMIPLE_PYTHON_MODULE)

target_link_libraries(translation PRIVATE pybind11::module ${PYTHON_LIBRARIES})
#target_link_libraries(translation PRIVATE ${PYTHON_LIBRARIES})
set_target_properties(translation PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")
