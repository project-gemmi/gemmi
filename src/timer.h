
#include <chrono>
#include <cstdio>

class Timer {
public:
  using Clock = std::chrono::steady_clock;
  Timer(bool enabled) : enabled_(enabled) {}
  void start() {
    if (enabled_)
      start_time_ = Clock::now();
  }
  bool enabled() const { return enabled_; }
  double count() const {
    std::chrono::duration<double> elapsed = Clock::now() - start_time_;
    return elapsed.count();
  };
  void print(const char* msg) const {
    if (enabled_)
      std::fprintf(stderr, "%s %g s.\n", msg, count());
  }
private:
  bool enabled_;
  std::chrono::time_point<Clock> start_time_;
};
