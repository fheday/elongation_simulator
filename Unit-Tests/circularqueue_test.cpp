#include <gtest/gtest.h>
#include <stdlib.h>
#include "../circularbuffer.h"
TEST(CircularQueue_int_tester, basic_tests) {
    utils::circular_buffer<int> queue(5);
    queue.put(3);
    queue.put(4);
    ASSERT_EQ(queue.get(), 3);
    ASSERT_EQ(queue.get(), 4);
    // for (int i = 0; i < (int) queue.size(); i++) queue.push(i);
    // for (int i = 0; i < (int) queue.size(); i++) ASSERT_EQ(queue.pull(), i);

}

TEST(CircularQueue_int_tester, Initialize) { EXPECT_EQ(1, 1); }
