cmake_minimum_required(VERSION 2.6)
project(gillespie)

#-fopt-info-vec-missed -g -fno-inline
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native -ftree-vectorize -Wall -DNDEBUG -flto")
set(CMAKE_CXX_FLAGS_DEFAULT "${CMAKE_CXX_FLAGS_DEFAULT} -O0 -g3 -Wall -Wextra -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -Wall -Wextra")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")

set(common_sources concentrationsreader.cpp mrna_reader.cpp)
set(ribosome_simulator_sources ribosomesimulator.cpp ${common_sources})
set(enlongation_simulator_sources enlogationsimulator.cpp  enlongation_codon.cpp initiationterminationcodon.cpp mrnaelement.cpp translation.cpp ${ribosome_simulator_sources})

add_library(util_lib STATIC ${common_sources})
add_library(ribosome_simulator_lib STATIC ${ribosome_simulator_sources})
add_library(enlongation_simulator_lib STATIC  ${enlongation_simulator_sources})

add_executable(calculateCodonsTimes calculateCodonsTimes.cpp)
target_link_libraries(calculateCodonsTimes PRIVATE ribosome_simulator_lib)
add_executable(run_translation run_translation.cpp)
target_link_libraries(run_translation PRIVATE enlongation_simulator_lib)

#add_executable(calculateCodonsTimes calculateCodonsTimes.cpp ${common_sources} ${ribosome_simulator_sources})
#add_executable(run_translation run_translation.cpp ${common_sources} ${enlongation_simulator_sources})


set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/experiments)
##python interface

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)
add_library(translation EXCLUDE_FROM_ALL MODULE translation.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp ribosomesimulator.cpp mrnaelement.cpp enlongation_codon.cpp)
target_compile_definitions(translation PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(translation PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(translation PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(ribosomesimulator EXCLUDE_FROM_ALL MODULE concentrationsreader.cpp ribosomesimulator.cpp)
target_compile_definitions(ribosomesimulator PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(ribosomesimulator PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(ribosomesimulator PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")


add_subdirectory(Unit-Tests)
