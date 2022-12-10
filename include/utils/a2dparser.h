#ifndef A2D_PARSER_H
#define A2D_PARSER_H

#include <iostream>
#include <string>
#include <vector>

namespace A2D {

/**
 * @brief Parse command line arguments.
 *
 * Example usage:
 *
  ArgumentParser parser(argc, argv);
  int nx = parser.parse_option("--nx", 5);
  double lx = parser.parse_option("--lx", 3.4);
  bool plot = parser.parse_option("--grad_check_only");
  std::string prefix = parser.parse_option("--prefix", std::string("results"));
 *
 */
class ArgumentParser {
 public:
  ArgumentParser(int argc, char* argv[]) {
    // Populate arguments
    for (int i = 1; i < argc; i++) {
      args.push_back(argv[i]);
    }

    // Initialize help info
    help += "Usage: ./[executable] ";
  }

  /**
   * @brief parse and option with value: --option=val
   *
   * Returns default value if option is not found, otherwise returns the value
   * from command line arguments.
   *
   * @tparam ValType value type
   * @param option option string
   * @param default_val default value
   */
  template <typename ValType>
  ValType parse_option(const std::string option, ValType default_val) {
    // Add to help info
    char info[256];
    if constexpr (std::is_same<ValType, double>::value) {
      std::snprintf(info, sizeof(info), "%s [%.2e] ", option.c_str(),
                    default_val);
    } else if constexpr (std::is_same<ValType, int>::value) {
      std::snprintf(info, sizeof(info), "%s [%d] ", option.c_str(),
                    default_val);
    } else if constexpr (std::is_same<ValType, std::string>::value) {
      std::snprintf(info, sizeof(info), "%s [%s] ", option.c_str(),
                    default_val.c_str());
    }
    help += info;

    // find option in cmd arguments
    auto it = std::find(args.begin(), args.end(), option);

    // If found, convert to desired type
    if (it != args.end()) {
      if constexpr (std::is_same<ValType, double>::value) {
        return std::stod(*(it + 1));
      } else if constexpr (std::is_same<ValType, int>::value) {
        return std::stoi(*(it + 1));
      } else if constexpr (std::is_same<ValType, std::string>::value) {
        return *(it + 1);
      }
      throw std::runtime_error("invalid type");
    } else {
      return default_val;
    }
  }

  /**
   * @brief parse and option without value: --option
   *
   * Returns true if option is found, otherwise false
   *
   * @param option option string
   */
  bool parse_option(const std::string option) {
    // Add to help info
    help += option;
    help += " ";

    // find option in cmd arguments
    auto it = std::find(args.begin(), args.end(), option);

    // If found, return true, otherwise return false
    if (it != args.end()) {
      return true;
    } else {
      return false;
    }
  }

  void help_info() {
    // find option in cmd arguments
    auto it1 = std::find(args.begin(), args.end(), "--help");
    auto it2 = std::find(args.begin(), args.end(), "-h");

    // If -h or --help in cmd argument, print help information and exit
    if (it1 != args.end() or it2 != args.end()) {
      std::cout << help << std::endl;
      exit(0);
    }
  };

 private:
  std::vector<std::string> args;
  std::vector<std::string> options;
  std::string help;
};

/**
 * @brief Helper function: save the command line arguments to a txt file
 *
 * @param argc number of arguments
 * @param argv arguments
 * @param txt_path the txt path
 */
void save_cmd(int argc, char* argv[], const std::string txt_path) {
  std::FILE* cmd_fp = std::fopen(txt_path.c_str(), "w+");
  for (int i = 0; i < argc; i++) {
    std::fprintf(cmd_fp, "%s ", argv[i]);
  }
  std::fclose(cmd_fp);
}

/**
 * @brief Check if target equals one of valid_vals
 */
template <typename EntryType>
void assert_option_in(std::string target, std::vector<EntryType> valid_vals) {
  auto domain_it = std::find(valid_vals.begin(), valid_vals.end(), target);
  if (domain_it == valid_vals.end()) {
    std::printf("Agrument value %s is invalid! Valid options are: ",
                target.c_str());
    for (auto it = valid_vals.begin(); it != valid_vals.end(); it++) {
      std::printf("%s, ", it->c_str());
    }
    std::printf("\b\b.\n");
    exit(-1);
  }
}

}  // namespace A2D
#endif