#pragma once

#include <boost/format.hpp>
#include <map>
#include <string>

namespace IGSV {

  template <typename StreamT>
  void print_timings(const std::map<std::string, double>& timings, //
                     StreamT& s,                                   //
                     bool print_total_hms = false);

  void log_time(double);
  void log_section(std::string);
  void log_error(std::string);
  void log_warning(std::string);
  void log_message(std::string);

  // command line progress bar
  void progress_bar(int i, // index of current element (0-based)
                    int N, // number of elements
                    int K, // number of progress bar parts
                    int& k // index of current progress bar part
  );

} // namespace IGSV
