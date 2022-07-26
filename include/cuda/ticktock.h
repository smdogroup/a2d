#include <chrono>
#include <cstdio>
#include <iostream>

#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>


// define the ticktock function for timing
#define TICK(msg) startTimer(msg);
#define TOCK(msg) reportTimer(msg);
#define __TICK__ \
  Timer timer;   \
  timer.startTimer();
#define __TOCK__ timer.stopTimer();

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::time_point<Clock> TimePoint;
typedef std::chrono::duration<double> Duration;
typename Clock::time_point startTime;
typename Clock::time_point endTime;

void startTimer(const char *msg) {
  std::cout << msg << ": start" << std::endl;
  startTime = Clock::now();
}

void stopTimer() { endTime = Clock::now(); }

double getElapsedTime() {
  stopTimer();
  Duration d = endTime - startTime;
  return d.count();
}

void reportTimer(const char *msg) {
  double elapsed = getElapsedTime();
  if (elapsed < 1e-6)
    std::cout << msg << ": " << elapsed * 1e6 << " us" << std::endl;
  else if (elapsed < 1e-3)
    std::cout << msg << ": " << elapsed * 1e3 << " ms" << std::endl;
  else
    std::cout << msg << ": " << elapsed << " s" << std::endl;
}

// create a timer class
class Timer {
 public:
  typename std::chrono::high_resolution_clock::time_point start;
  typename std::chrono::high_resolution_clock::time_point end;

  void startTimer() {
    start = std::chrono::high_resolution_clock::now();
    std::cout << "Total: start" << std::endl;
  }

  void stopTimer() { end = std::chrono::high_resolution_clock::now(); }

  void reset() { start = std::chrono::high_resolution_clock::now(); }

  double getElapsedTime() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(end -
                                                                     start)
        .count();
  }

  ~Timer() {
    std::cout << "Total time: " << getElapsedTime() << " s" << std::endl;
  }
};
