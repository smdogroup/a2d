#ifndef A2D_PROFILER_H
#define A2D_PROFILER_H

#include <chrono>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

namespace A2D {

// Timer options
static bool TIMER_IS_ON = true;
static std::string TIMER_OUTPUT_FILE = "profile.log";
static int TIMER_NUM_SPACE = 4;

// Internally used variables for the timer
static int timer_recursion_counter = 0;
static int timer_log_touched = false;

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
    if (TIMER_IS_ON) {
      fun_name = _fun_name;

      if (!timer_log_touched) {
        fp = std::fopen(TIMER_OUTPUT_FILE.c_str(),
                        "w+");  // Clean the file, if exist
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::fprintf(fp,
                     "========================\n%s========================\n",
                     std::ctime(&now_time));
        timer_log_touched = true;
      } else {
        fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
      }

      if (fun_name.empty()) {
        std::fprintf(fp, "%sentering an anonymous scope\n",
                     std::string(TIMER_NUM_SPACE * timer_recursion_counter, ' ')
                         .c_str());
      } else {
        std::fprintf(
            fp, "%s%s executed\n",
            std::string(TIMER_NUM_SPACE * timer_recursion_counter, ' ').c_str(),
            fun_name.c_str());
      }
      std::fclose(fp);
      timer_recursion_counter++;
      t_start = std::chrono::steady_clock::now();
    }
  }

  ~Timer() {
    if (TIMER_IS_ON) {
      auto t_end = std::chrono::steady_clock::now();
      auto t_elapse =
          std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start)
              .count();

      timer_recursion_counter--;
      fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
      if (fun_name.empty()) {
        std::fprintf(
            fp, "%-80s(%lld ms)\n",
            (std::string(TIMER_NUM_SPACE * timer_recursion_counter, ' ') +
             "exiting an anonymous scope")
                .c_str(),
            t_elapse);
      } else {
        std::fprintf(
            fp, "%-80s(%lld ms)\n",
            (std::string(TIMER_NUM_SPACE * timer_recursion_counter, ' ') +
             fun_name + " exits")
                .c_str(),
            t_elapse);
      }
      std::fclose(fp);
    }
  }

 private:
  std::FILE *fp;
  std::chrono::time_point<std::chrono::steady_clock> t_start;
  std::string fun_name;
};

}  // namespace A2D

#endif  // A2D_PROFILER_H