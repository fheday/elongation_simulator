cmake_minimum_required(VERSION 2.6)
project(gillespie)

add_definitions(-DCMAKE_BUILD)

add_executable(gillespie main.cpp)
install(TARGETS gillespie RUNTIME DESTINATION bin)


# cmake_minimum_required(VERSION 2.8.12)
# project(example)
# add_subdirectory(pybind11)
# pybind11_add_module(ribosomesimulator ribosomesimulator.cpp)
