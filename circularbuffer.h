/*
 * @file Circularqueue.h
 *
 *  Created on: 30 Mar 2020
 *
 * @brief Implementation of custom circular buffer
 *
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#ifndef CIRCULARQUEUE_H_
#define CIRCULARQUEUE_H_

#include <cstddef>
#include <memory>
#include <vector>
namespace utils
{

  template <class T>
  class circular_buffer
  {
  public:
    explicit circular_buffer(size_t size) : buf_(std::unique_ptr<T[]>(new T[size])),
                                            max_size_(size)
    {
    }

    void put(T item) {
      buf_[head_] = item;
      if (full_)
      {
        tail_ = fast_mod((tail_ + 1) , max_size_);
      }

      head_ = fast_mod((head_ + 1) , max_size_);
      full_ = head_ == tail_;
    }

    T get() {
      if (empty()) return T();

      auto val = buf_[tail_];
      full_ = false;
      tail_ = fast_mod((tail_ + 1) , max_size_);

      return val;
    }

    T& peek_back() {
      if (head_ == 0) {
        return buf_[max_size_ - 1];
      } else {
        return buf_[head_ - 1];
      }
      
    }

    void replace(T old_item, T new_item) {
      std::size_t index = tail_;
      for (size_t i = 0; i < size(); i++) {
        if (buf_[index] == old_item) {
          buf_[index] = new_item;
          return;
        }
        index = fast_mod((index + 1) , max_size_);
      }
    }

    std::vector<T> get_vector(const bool inverted) {
      std::size_t size_ = size();
      std::vector<T> result(size_);
      for (std::size_t i = 0; i < result.size(); i++) {
        std::size_t ind = inverted ? size_ - 1 - i: i;
        result[ind] = buf_[fast_mod((tail_ + i) , max_size_)];
      }
      return result;
    }

    void reset() {
      head_ = tail_;
      full_ = false;
    }

    bool empty() const {
      return (!full_ && (head_ == tail_));
    }

    bool full() const {
      return full_;
    }

    size_t capacity() const {
      return max_size_;
    }

    size_t size() const {
      size_t size = max_size_;

      if (!full_) {
        if (head_ >= tail_) {
          size = head_ - tail_;
        } else {
          size = max_size_ + head_ - tail_;
        }
      }

      return size;
    }

  private:
    std::unique_ptr<T[]> buf_;
    size_t head_ = 0;
    size_t tail_ = 0;
    const size_t max_size_;
    bool full_ = 0;

    size_t fast_mod(const size_t input, const size_t ceil) {
      // get the modulus to be put into the buffer.
      // input has to be positive and to greater than ceil.
      return input < ceil ? input : input - ceil;
}
  };

} /* namespace utils */

#endif /* CIRCULARQUEUE_H_ */
