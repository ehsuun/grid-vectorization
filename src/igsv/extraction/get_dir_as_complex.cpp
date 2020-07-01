#include <igsv/extraction/get_dir_as_complex.h>
namespace IGSV {
  std::complex<int> get_dir_as_complex(int r) {
    switch (r) {
    case 0:
      return std::complex<int>(1, 0);
    case 1:
      return std::complex<int>(0, 1);
    case 2:
      return std::complex<int>(-1, 0);
    default:
      return std::complex<int>(0, -1);
    }
  }
} // namespace IGSV
