cmake_minimum_required(VERSION 2.6)
project(gillespie)

add_executable(gillespie main.cpp)

install(TARGETS gillespie RUNTIME DESTINATION bin)
