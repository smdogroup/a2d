#include <chrono>
#include <thread>

#include "utils/a2dprofiler.h"

int main() {
  A2D::Timer::set_threshold_ms(50.0);
  {
    A2D::Timer t("sleep(20)");
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    {
      A2D::Timer t;
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      {
        A2D::Timer t("sleep(40)");
        std::this_thread::sleep_for(std::chrono::milliseconds(40));
      }
    }
  }
  return 0;
}