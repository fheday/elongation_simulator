cmake_minimum_required(VERSION 3.9.2)
project(gillespie CXX)

#-fopt-info-vec-missed -g -fno-inline
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native -ftree-vectorize -Wall -DNDEBUG -flto -g")
set(CMAKE_CXX_FLAGS_DEFAULT "${CMAKE_CXX_FLAGS_DEFAULT} -O0 -g3 -Wall -Wextra -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined,address")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -Wall -Wextra")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")


set(common_sources concentrationsreader.cpp mrna_reader.cpp)
set(ribosome_simulator_sources ribosomesimulator.cpp ${common_sources})
set(elongation_simulator_sources   elongation_codon.cpp initiationterminationcodon.cpp mrnaelement.cpp translation.cpp ${ribosome_simulator_sources} elongation_simulation_manager.cpp elongation_simulation_processor.cpp)
add_library(util_lib STATIC ${common_sources})
add_library(ribosome_simulator_lib STATIC ${ribosome_simulator_sources})
add_library(elongation_simulator_lib STATIC  ${elongation_simulator_sources})

add_executable(calculateCodonsTimes calculateCodonsTimes.cpp)
target_link_libraries(calculateCodonsTimes PRIVATE ribosome_simulator_lib)
add_executable(run_translation run_translation.cpp)
target_link_libraries(run_translation PRIVATE elongation_simulator_lib)

##python interface
set (PYTHON_EXECUTABLE "python3.10")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories("eigen-3.3.7/")
include_directories("jsoncpp/")
add_subdirectory("pybind11/")

add_library(translation EXCLUDE_FROM_ALL MODULE translation.cpp elongation_simulation_manager.cpp elongation_simulation_processor.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp ribosomesimulator.cpp mrnaelement.cpp elongation_codon.cpp)
target_compile_definitions(translation PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(translation PRIVATE pybind11::module ${PYTHON_LIBRARIES} jsoncpp)
set_target_properties(translation PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(ribosomesimulator EXCLUDE_FROM_ALL MODULE concentrationsreader.cpp ribosomesimulator.cpp)
target_compile_definitions(ribosomesimulator PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(ribosomesimulator PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(ribosomesimulator PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")


file(COPY concentrations/Saccharomyces_cerevisiae.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/concentrations)

if(BUILD_TESTS)
    add_subdirectory(Unit-Tests)
endif(BUILD_TESTS)
