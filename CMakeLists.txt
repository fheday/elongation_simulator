cmake_minimum_required(VERSION 2.6)
project(gillespie)

set(common_sources concentrationsreader.cpp concentrationsreader.h enlogationsimulator.cpp enlogationsimulator.h gillespie.cpp gillespie.h mrna_reader.cpp mrna_reader.h ratecalculator.cpp ratecalculator.h reactionsset.cpp reactionsset.h ribosomesimulator.cpp ribosomesimulator.h ribosomesimulator.cpp enlongation_codon.h enlongation_codon.cpp initiationterminationcodon.h  initiationterminationcodon.cpp mrnaelement.h  
translation.h translation.cpp)

add_definitions(-DCMAKE_BUILD)

add_executable(gillespie main.cpp ${common_sources})
add_executable(calculateCodonsTimes calculateCodonsTimes.cpp ${common_sources})
add_executable(run_translation run_translation.cpp ${common_sources})
install(TARGETS gillespie RUNTIME DESTINATION bin)


# cmake_minimum_required(VERSION 2.8.12)
# project(example)
# add_subdirectory(pybind11)
# pybind11_add_module(ribosomesimulator ribosomesimulator.cpp)
