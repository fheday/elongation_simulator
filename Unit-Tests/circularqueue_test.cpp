#include <gtest/gtest.h>
#include <stdlib.h>
#include "../circularbuffer.h"
TEST(CircularQueue_int_tester, basic_test) {
    utils::circular_buffer<int> queue(5);
    for (int i = 0; i < (int) queue.size(); i++) queue.put(i);
    for (int i = 0; i < (int) queue.size(); i++) ASSERT_EQ(queue.get(), i);
}

TEST(CircularQueue_int_tester, CircularQueue_int_to_vector) {
    utils::circular_buffer<int> queue(5);
    for (int i = 0; i < (int) queue.size(); i++) queue.put(i);
    auto queue_vector = queue.get_vector();
    for (int i = 0; i < (int) queue.size(); i++) ASSERT_EQ(queue_vector[i], i);
}

TEST(CircularQueue_vector_tester, basic_test) {
    utils::circular_buffer<std::vector<int>> queue(3);
    std::vector<int> x0{1, 2, 3};
    std::vector<int> x1{4, 5, 6};
    queue.put(x0);
    queue.put(x1);
    int count = 1;
    for (std::size_t i = 0; i < queue.size(); i++){
        std::vector<int> elem = queue.get();
        for (std::size_t j = 0; j < elem.size(); j++) ASSERT_EQ(elem[j], count++);
    }
}

TEST(CircularQueue_vector_tester, queue_vect_to_vect_of_vect_of_int) {
    utils::circular_buffer<std::vector<int>> queue(3);
    std::vector<int> x0{1, 2, 3};
    std::vector<int> x1{4, 5, 6};
    queue.put(x1);
    queue.put(x0);
    int count = 1;
    auto queue_vector = queue.get_vector();
    for (std::size_t i = 0; i < queue.size(); i++){
        std::vector<int> elem = queue_vector[i];
        for (std::size_t j = 0; j < elem.size(); j++) ASSERT_EQ(elem[j], count++);
    }
}
