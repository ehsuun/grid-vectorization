#pragma once
#define STRW 40 // command line output : lhs string width

#undef PRINT_SETTINGS
#undef PRINT_TIMINGS
#undef PRINT_DEBUGINFO

#define PRINT_SETTINGS 1
#define PRINT_TIMINGS 1
#define PRINT_DEBUGINFO 0

// https://stackoverflow.com/a/9158263/1606707
// the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define TERMINAL_FMT_RESET "\033[0m"
#define TERMINAL_FMT_BLACK "\033[30m"              /* Black */
#define TERMINAL_FMT_RED "\033[31m"                /* Red */
#define TERMINAL_FMT_GREEN "\033[32m"              /* Green */
#define TERMINAL_FMT_YELLOW "\033[33m"             /* Yellow */
#define TERMINAL_FMT_BLUE "\033[34m"               /* Blue */
#define TERMINAL_FMT_MAGENTA "\033[35m"            /* Magenta */
#define TERMINAL_FMT_CYAN "\033[36m"               /* Cyan */
#define TERMINAL_FMT_WHITE "\033[37m"              /* White */
#define TERMINAL_FMT_BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define TERMINAL_FMT_BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define TERMINAL_FMT_BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define TERMINAL_FMT_BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define TERMINAL_FMT_BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define TERMINAL_FMT_BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define TERMINAL_FMT_BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define TERMINAL_FMT_BOLDWHITE "\033[1m\033[37m"   /* Bold White */

#define DEBUG_PRINT_INFO(fmt, ...)              \
  do {                                          \
    if (PRINT_DEBUGINFO)                        \
      fprintf(stdout, fmt "\n", ##__VA_ARGS__); \
  } while (0)

#define DEBUG_PRINT_WARNING(fmt, ...)                          \
  do {                                                         \
    if (PRINT_DEBUGINFO)                                       \
      fprintf(stderr, "warning: " fmt "\n", ##__VA_ARGS__); \
  } while (0)

#define DEBUG_PRINT_FN(fmt, ...)                                                                     \
  do {                                                                                               \
    if (PRINT_DEBUGINFO)                                                                             \
      fprintf(stdout, "\n" TERMINAL_FMT_MAGENTA "#### " fmt "\n" TERMINAL_FMT_RESET, ##__VA_ARGS__); \
  } while (0)

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#define BOOST_VARIANT_USE_RELAXED_GET_BY_DEFAULT
#define BOOST_PYTHON_STATIC_LIB
#define CGAL_INTERSECTION_VERSION 2
