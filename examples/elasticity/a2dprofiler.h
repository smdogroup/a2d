#ifndef A2D_PROFILER_H
#define A2D_PROFILER_H

#include <chrono>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

namespace A2D {

static int TIMER_RECURSION_COUNTER = 0;
static int TIMER_LOG_TOUCHED = false;
static std::string TIMER_OUTPUT_FILE = "profile.log";

/**
 * @brief A very simple scope-based timer.
 *
 * Usage:
 * {
 *   Timer t("Very brief info about what is timed");
 *
 *   // Code snippet to be timed
 *   // ...
 *
 * }
 */
class Timer {
 public:
  Timer(std::string _fun_name = "") {
    fun_name = _fun_name;

    if (!TIMER_LOG_TOUCHED) {
      fp = std::fopen(TIMER_OUTPUT_FILE.c_str(),
                      "w+");  // Clean the file, if exist
      auto now = std::chrono::system_clock::now();
      std::time_t now_time = std::chrono::system_clock::to_time_t(now);
      std::fprintf(fp, "========================\n%s========================\n",
                   std::ctime(&now_time));
      TIMER_LOG_TOUCHED = true;
    } else {
      fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
    }

    if (fun_name.empty()) {
      std::fprintf(
          fp, "%sentering an anonymous scope\n",
          std::string(NUM_SPACE_EACH_RECUR * TIMER_RECURSION_COUNTER, ' ')
              .c_str());
    } else {
      std::fprintf(
          fp, "%s%s executed\n",
          std::string(NUM_SPACE_EACH_RECUR * TIMER_RECURSION_COUNTER, ' ')
              .c_str(),
          fun_name.c_str());
    }
    std::fclose(fp);
    TIMER_RECURSION_COUNTER++;
    t_start = std::chrono::steady_clock::now();
  }

  ~Timer() {
    auto t_end = std::chrono::steady_clock::now();
    auto t_elapse =
        std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start)
            .count();

    TIMER_RECURSION_COUNTER--;
    fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
    if (fun_name.empty()) {
      std::fprintf(
          fp, "%-80s(%lld ms)\n",
          (std::string(NUM_SPACE_EACH_RECUR * TIMER_RECURSION_COUNTER, ' ') +
           "exiting an anonymous scope")
              .c_str(),
          t_elapse);
    } else {
      std::fprintf(
          fp, "%-80s(%lld ms)\n",
          (std::string(NUM_SPACE_EACH_RECUR * TIMER_RECURSION_COUNTER, ' ') +
           fun_name + " exits")
              .c_str(),
          t_elapse);
    }
    std::fclose(fp);
  }

  static void set_output_file(const std::string file_path) {
    TIMER_OUTPUT_FILE = file_path;
  };

 private:
  static const int NUM_SPACE_EACH_RECUR = 4;
  std::FILE *fp;
  std::chrono::time_point<std::chrono::steady_clock> t_start;
  std::string fun_name;
};

}  // namespace A2D

#endif  // A2D_PROFILER_H