#ifndef UTILITIES_HPP_INCLUDED
#define UTILITIES_HPP_INCLUDED
#include <chrono>
#include <vector>

namespace boat {

template <int... Is>
struct seq {};

template <int N, int... Is>
struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

template <int... Is>
struct gen_seq<0, Is...> : seq<Is...> {};

class Timer {
 public:
  Timer() : start_(std::chrono::high_resolution_clock::now()) {}

  void stop() { end_ = std::chrono::high_resolution_clock::now(); }

  double time_interval_nanoseconds() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end_ - start_)
        .count();
  }
  double time_interval_microseconds() {
    return time_interval_nanoseconds() / 1000.0;
  }
  double time_interval_miliseconds() {
    return time_interval_nanoseconds() / 1000000.0;
  }
  double time_interval_seconds() {
    return time_interval_nanoseconds() / 1000000000.0;
  }

  double time_since_start_nanoseconds() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::high_resolution_clock::now() - start_).count();
  }
  double time_since_start_microseconds() {
    return time_since_start_nanoseconds() / 1000.0;
  }
  double time_since_start_miliseconds() {
    return time_since_start_nanoseconds() / 1000000.0;
  }
  double time_since_start_seconds() {
    return time_since_start_nanoseconds() / 1000000000.0;
  }

 private:
  std::chrono::high_resolution_clock::time_point start_;
  std::chrono::high_resolution_clock::time_point end_;
};
}

#endif  // UTILITIES_HPP_INCLUDED
