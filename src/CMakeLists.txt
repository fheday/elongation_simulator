cmake_minimum_required(VERSION 3.10)
project(gillespie CXX)

#-fopt-info-vec-missed -g -fno-inline
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -ffast-math -march=native -ftree-vectorize -Wall -DNDEBUG -flto=auto -g")
set(CMAKE_CXX_FLAGS_DEFAULT "${CMAKE_CXX_FLAGS_DEFAULT} -O0 -g3 -Wall -Wextra -fsanitize=leak -fno-omit-frame-pointer -fsanitize=undefined,address")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -Wall -Wextra")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")


set(common_sources concentrationsreader.cpp mrna_reader.cpp)
set(codon_simulator_sources codon_simulator.cpp ${common_sources})
set(elongation_simulator_sources   elongation_codon.cpp initiationterminationcodon.cpp mrnaelement.cpp sequence_simulator.cpp ${codon_simulator_sources} elongation_simulation_manager.cpp elongation_simulation_processor.cpp)
add_library(util_lib STATIC ${common_sources})
add_library(codon_simulator_lib STATIC ${codon_simulator_sources})
add_library(elongation_simulator_lib STATIC  ${elongation_simulator_sources})

add_executable(calculateCodonsTimes calculateCodonsTimes.cpp)
target_link_libraries(calculateCodonsTimes PRIVATE codon_simulator_lib)
add_executable(run_translation run_translation.cpp)
target_link_libraries(run_translation PRIVATE elongation_simulator_lib)

# #python interface
# set (PYTHON_EXECUTABLE "python3.10")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories("eigen-3.3.7/")
include_directories("jsoncpp/")

execute_process(COMMAND pybind11-config --cmakedir RESULT_VARIABLE pybind11_DIR)
# include(CMakePrintHelpers)
# cmake_print_variables(pybind11_DIR)

# add_subdirectory(execute_process(COMMAND pybind11-config --cmakedir))
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Julia interface
if(BUILD_JULIA)
    find_package(JlCxx)
    get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
    get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
endif(BUILD_JULIA)

add_library(sequence_simulator EXCLUDE_FROM_ALL MODULE sequence_simulator.cpp elongation_simulation_manager.cpp elongation_simulation_processor.cpp mrna_reader.cpp initiationterminationcodon.cpp
concentrationsreader.cpp codon_simulator.cpp mrnaelement.cpp elongation_codon.cpp)
target_compile_definitions(sequence_simulator PRIVATE -DCOMIPLE_PYTHON_MODULE) # enable python pathway in code.
target_link_libraries(sequence_simulator PRIVATE pybind11::module ${PYTHON_LIBRARIES} jsoncpp)
set_target_properties(sequence_simulator PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_library(codon_simulator EXCLUDE_FROM_ALL MODULE concentrationsreader.cpp codon_simulator.cpp)

target_compile_definitions(codon_simulator PRIVATE -DCOMIPLE_PYTHON_MODULE ) # enable python pathway in code.
# Julia interface
if(BUILD_JULIA)
    target_compile_definitions(codon_simulator PRIVATE -DCOMIPLE_JULIA_MODULE ) # enable Julia pathway in code.
    target_link_libraries(codon_simulator PRIVATE JlCxx::cxxwrap_julia)
endif(BUILD_JULIA)

target_link_libraries(codon_simulator PRIVATE pybind11::module ${PYTHON_LIBRARIES})
set_target_properties(codon_simulator PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                          SUFFIX "${PYTHON_MODULE_EXTENSION}")


file(COPY ../concentrations/Saccharomyces_cerevisiae.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/concentrations)

if(BUILD_TESTS)
  add_subdirectory(Unit-Tests)
endif(BUILD_TESTS)
