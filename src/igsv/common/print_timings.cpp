#include <igsv/common/defs.h>
#include <igsv/common/print_timings.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace IGSV {
  template <typename StreamT>
  void print_timings(const std::map<std::string, double>& timings, //
                     StreamT& s,                                   //
                     bool print_total_hms) {
    s << std::endl << std::string(STRW + 20, '~') << std::endl;
    double total_time = 0.;
    for (const auto& t : timings) {
      total_time += t.second;
      s << "  " << (t.first + std::string(STRW - 5 - t.first.length(), ' ') + " : ") //
        << std::fixed << std::setprecision(4) << t.second << " sec" << std::endl;
    }
    s << "  " << std::string(STRW + 8, '-') << std::endl;
    s << "  " << ("::TOTAL" + std::string(STRW - 12, ' ') + " : ") //
      << std::fixed << std::setprecision(4) << total_time << " sec" << std::endl;
    if (print_total_hms) {
      const int total_hh  = std::floor(total_time / 3600.0);
      const int total_mm  = std::floor(total_time / 60.0);
      const int total_ss  = std::round(total_time);
      const int remain_mm = total_mm - total_hh * 60;
      const int remain_ss = total_ss - total_mm * 60;
      s << std::string(STRW - 2, ' ') << "[ " << total_hh << "h, " << remain_mm << "m, " << remain_ss << "s ]"
        << std::endl;
    }
    s << std::string(STRW + 20, '~') << std::endl << std::endl;
  }

  void log_section(std::string _msg) {
    std::cout << boost::format(TERMINAL_FMT_GREEN "\n#### %s" TERMINAL_FMT_RESET "\n") % _msg;
  }

  void log_error(std::string _msg) {
    std::cout << boost::format(TERMINAL_FMT_BOLDRED "error: %s" TERMINAL_FMT_RESET "\n") % _msg;
  }

  void log_warning(std::string _msg) {
    std::cout << boost::format(TERMINAL_FMT_BOLDYELLOW "warning: %s" TERMINAL_FMT_RESET "\n") % _msg;
  }

  void log_time(double time) {
    std::cout << boost::format(TERMINAL_FMT_BLUE "[%.2f sec]" TERMINAL_FMT_RESET "\n") % time;
  }

  void log_message(std::string _msg) { std::cout << boost::format("%s\n") % _msg; }

  void progress_bar(int i, // index of current element (0-based)
                    int N, // number of elements
                    int K, // number of progress bar parts
                    int& k // index of current progress bar part
  ) {
    if (i == -1) {
      std::cout << "[" << std::string(K, ' ') << "]";
      return;
    }

    if ((float)i >= (float)(k * (N - 1)) / (float)(K - 1)) {
      ++k;
      std::cout << std::string(K + 2, '\b') //
                << "["                      //
                << std::string(k, '*')      //
                << std::string(K - k, ' ')  //
                << "]"                      //
                << std::flush;
    }
    if (i == N - 1)
      std::cout << std::endl;
  }

} // namespace IGSV

// explicit template instantiation
#include <fstream>
#include <iostream>
#include <sstream>

template void IGSV::print_timings<std::basic_ostream<char, std::char_traits<char>>>(
    // std::basic_string<char, std::char_traits<char>, std::allocator<char>> const&,
    std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char>>, double,
             std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char>>>,
             std::allocator<std::pair<std::basic_string<char, std::char_traits<char>, std::allocator<char>> const,
                                      double>>> const&,
    std::basic_ostream<char, std::char_traits<char>>&, bool);

template void IGSV::print_timings<std::basic_ofstream<char, std::char_traits<char>>>(
    // std::basic_string<char, std::char_traits<char>, std::allocator<char>> const&,
    std::map<std::basic_string<char, std::char_traits<char>, std::allocator<char>>, double,
             std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char>>>,
             std::allocator<std::pair<std::basic_string<char, std::char_traits<char>, std::allocator<char>> const,
                                      double>>> const&,
    std::basic_ofstream<char, std::char_traits<char>>&, bool);
