
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(translatortester_test  translatortester_test.cpp)
target_link_libraries(translatortester_test ${GTEST_BOTH_LIBRARIES} util_lib ribosome_simulator_lib enlongation_simulator_lib)

add_test(AllTestsInFoo translatortester_test)
