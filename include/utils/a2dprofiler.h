#ifndef A2D_PROFILER_H
#define A2D_PROFILER_H

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace A2D {

struct TimerEntry {
  std::string msg;
  char type;  // '(' or ')'
  double t;   // in ms
};

/**
 * @brief A very simple scope-based timer that measures time between creation
 *        and destroy of the timer object.

 * Usage:
 * {
 *   Timer t("Very brief info about what is timed");
 *
 *   // Code snippet to be timed
 *   // ...
 * }
 */
class Timer {
 public:
  Timer(std::string _fun_name = "unknown") {
    if (TIMER_IS_ON) {
      fun_name = _fun_name;

      if (!timer_log_touched) {
        timer_log_touched = true;

        // Clean log files, if exist
        fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "w+");
        fp2 = std::fopen(TIMER_SHORT_OUTPUT_FILE.c_str(), "w+");

        // Get time stamp
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);

        // Write header
        std::fprintf(fp,
                     "========================\n%s========================\n",
                     std::ctime(&now_time));
        std::fprintf(fp2,
                     "========================\n%sshort log, threshold: %.1f "
                     "ms\n========================\n",
                     std::ctime(&now_time), TIMER_THRESHOLD_MS);
        std::fclose(fp2);

      } else {
        fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
      }

      // Assemble the line of message
      std::string message =
          std::string(TIMER_TAB * counter, ' ') + fun_name + " executed";

      // Flush the long log
      std::fprintf(fp, "%s\n", message.c_str());
      std::fclose(fp);

      // Save data for short log
      buffer.push_back(TimerEntry{message, '(', 0.0});

      counter++;
      t_start = std::chrono::steady_clock::now();
    }

    return;
  }

  ~Timer() {
    if (TIMER_IS_ON) {
      auto t_end = std::chrono::steady_clock::now();
      auto t_elapse =
          std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start)
              .count();

      counter--;

      // Assemble line of message
      char _message[256];
      snprintf(
          _message, 256, "%-80s(%.2f ms)",
          (std::string(TIMER_TAB * counter, ' ') + fun_name + " exits").c_str(),
          (double)t_elapse);
      std::string message(_message);

      // Flush the long log
      fp = std::fopen(TIMER_OUTPUT_FILE.c_str(), "a+");
      std::fprintf(fp, "%s\n", message.c_str());
      std::fclose(fp);

      // Save data for short log
      buffer.push_back(TimerEntry{message, ')', (double)t_elapse});

      // Once the most outer function returns, we filter the buffer such
      // that we only keep entry pairs whose elapse time is above threshold
      // We assume that lines in buffer are all paired
      if (counter == 0) {
        std::vector<int> istart, keep;
        int idx = 0;
        for (auto it = buffer.begin(); it != buffer.end(); it++) {
          if (it->type == '(') {
            istart.push_back(idx);
          }
          if (it->type == ')') {
            if (!istart.empty()) {
              int start_idx = istart.back();
              istart.pop_back();
              if (it->t > TIMER_THRESHOLD_MS) {
                keep.push_back(start_idx);
                keep.push_back(idx);  // End idx
              }
            }
          }
          idx++;
        }

        //  Now, we only keep the entries for expensive function calls
        std::sort(keep.begin(), keep.end());
        fp2 = std::fopen(TIMER_SHORT_OUTPUT_FILE.c_str(), "a+");
        for (auto ik = keep.begin(); ik != keep.end(); ik++) {
          std::fprintf(fp2, "%s\n", buffer[*ik].msg.c_str());
        }
        std::fclose(fp2);

        // Clean the buffer
        buffer.clear();
      }
    }

    return;
  }

  static void on() { TIMER_IS_ON = true; }
  static void off() { TIMER_IS_ON = false; }

  static void set_log_path(std::string log_path) {
    TIMER_OUTPUT_FILE = log_path;
  }

  static void set_short_log_path(std::string short_log_path) {
    TIMER_SHORT_OUTPUT_FILE = short_log_path;
  }

  static void set_threshold_ms(double t) { TIMER_THRESHOLD_MS = t; }

 private:
  std::FILE *fp, *fp2;
  std::chrono::time_point<std::chrono::steady_clock> t_start;
  std::string fun_name;
  inline static int counter = 0;
  inline static bool timer_log_touched = false;
  inline static int TIMER_TAB = 4;
  inline static std::vector<TimerEntry> buffer;

  // Timer options
  inline static std::string TIMER_OUTPUT_FILE = "profile.log";
  inline static std::string TIMER_SHORT_OUTPUT_FILE = "profile_short.log";
  inline static double TIMER_THRESHOLD_MS = 1.0;
  inline static bool TIMER_IS_ON = true;
};

}  // namespace A2D

#endif  // A2D_PROFILER_H