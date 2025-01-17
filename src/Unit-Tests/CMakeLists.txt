
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(sequence_simulator_test  sequence_simulator_test.cpp)
add_executable(circularqueue_test circularqueue_test.cpp)
add_executable(sequence_simulator_manager_test sequence_simulator_manager_test.cpp)
add_executable(mrna_reader_test mrna_reader_test.cpp)
target_link_libraries(sequence_simulator_test ${GTEST_BOTH_LIBRARIES} util_lib codon_simulator_lib elongation_simulator_lib)
target_link_libraries(circularqueue_test ${GTEST_BOTH_LIBRARIES})
target_link_libraries(sequence_simulator_manager_test PRIVATE Threads::Threads ${GTEST_BOTH_LIBRARIES} jsoncpp elongation_simulator_lib)
target_link_libraries(mrna_reader_test ${GTEST_BOTH_LIBRARIES} util_lib)

add_test(AllTestsInFoo translatortester_test)
add_test(CircularQueue_test circularqueue_test)