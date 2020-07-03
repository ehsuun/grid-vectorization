#include <igsv/common/Options.h>
#include <igsv/common/defs.h>

#include <boost/filesystem.hpp>
#include <iostream>

namespace po = boost::program_options;

namespace IGSV {

  // ==========================================================================

  bool Options::parse_args(int argc, char* argv[]) {

    po::options_description all_options{ "options" };
    all_options.add_options()
        //
        ("help,h", "help screen")
        //
        ("input,i", po::value<std::string>(&in.filename_input), //
         "input filename (png)")
        //
        ("output,o", po::value<std::string>(&in.filename_output), //
         "output filename (svg)")
        //
        ("mask,m", po::value<std::string>(&in.filename_mask)->default_value(""), //
         "mask filename (png)")
        //
        ("threshold,b", po::value<int>(&in.threshold)->default_value(165), //
         "binarization threshold")
        //
        ("narrow-band-radius,n", po::value<float>(&tri.narrow_band_radius)->default_value(0.3),
         "maximal distance (relative to the average stroke width) of narrow band triangles from dark pixels")
        //
        ("scale-multiplier,s", po::value<float>(&uv.scale_multiplier)->default_value(1.0), //
         "scale multiplier")
        //
        ("mask-factor,f", po::value<float>(&in.mask_factor)->default_value(2.0), //
         "scale factor for user mask")
        //
        ("deserialize,d", po::value<std::string>(&in.filename_state), //
         "filename of a state used for deserialization")
        //
        ("no-gui", po::bool_switch(&misc.no_gui), //
         "If specified, the program will run in terminal-only mode. The parametrization will be serialized to the '${INPUT_FILENAME}.state' file, e.g. '../data/fish.png.state'.")
        //
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_options), vm);

    if (vm.count("help")) {
      std::cout << TERMINAL_FMT_YELLOW << R"IGSVTITLE(
  __    _______      _______.____    ____
 |  |  /  _____|    /       |\   \  /   /
 |  | |  |  __     |   (----` \   \/   /
 |  | |  | |_ |     \   \      \      /
 |  | |  |__| | .----)   |      \    /
 |__|  \______| |_______/        \__/

   Integer Grid Sketch Vectorization

)IGSVTITLE" << TERMINAL_FMT_RESET
                << all_options << std::endl;
      return false;
    }

    po::notify(vm);

    if (vm.count("input") == 0) {
      std::cerr << "--input [-i] not specified!\nabort." << std::endl;
      return false;
    }

    if (vm.count("output") == 0) {
      // if no output is specified, save to where the input is
      in.filename_output = boost::filesystem::path(in.filename_input).replace_extension("svg").string();
    }

    return true;
  }

  // ==========================================================================

  template <typename StreamT>
  void Options::print(StreamT& s) const {
    s << "options:"                                                       //
      << "\n  in.mask_factor            : " << in.mask_factor             //
      << "\n  in.threshold              : " << in.threshold               //
      << "\n  tri.adaptive              : " << tri.adaptive               //
      << "\n  tri.estimate_max_area     : " << tri.estimate_max_area      //
      << "\n  tri.max_area              : " << tri.max_area               //
      << "\n  tri.max_area_coeff        : " << tri.max_area_coeff         //
      << "\n  tri.narrow_band_radius    : " << tri.narrow_band_radius     //
      << "\n  tri.wdist_sigma           : " << tri.wdist_sigma            //
      << "\n  ff.use_weighted_reg       : " << ff.use_weighted_reg        //
      << "\n  ff.wsmooth                : " << ff.wsmooth                 //
      << "\n  ff.wtangent               : " << ff.wtangent                //
      << "\n  ff.wortho                 : " << ff.wortho                  //
      << "\n  ff.tangent_ratio_threshold: " << ff.tangent_ratio_threshold //
      << "\n  ff.tangent_ratio_streamlen: " << ff.tangent_ratio_streamlen //
      << "\n  uv.scale_multiplier       : " << uv.scale_multiplier        //
      << "\n  uv.wsnap                  : " << uv.wsnap << std::endl;
  }

  // ===========================================================================

  void Options::InitSerialization() {

    Add(in.mask_factor, "in__mask_factor");
    Add(in.threshold, "in__threshold");

    Add(tri.adaptive, "tri__adaptive");
    Add(tri.estimate_max_area, "tri__estimate_max_area");
    Add(tri.max_area, "tri__max_area");
    Add(tri.max_area_coeff, "tri__max_area_coeff");
    Add(tri.narrow_band_radius, "tri__narrow_band_radius");
    Add(tri.wdist_sigma, "tri__wdist_sigma");

    Add(ff.use_weighted_reg, "ff__use_weighted_reg");
    Add(ff.wsmooth, "ff__wsmooth");
    Add(ff.wtangent, "ff__wtangent");
    Add(ff.wortho, "ff__wortho");
    Add(ff.tangent_ratio_threshold, "ff__tangent_ratio_threshold");
    Add(ff.tangent_ratio_streamlen, "ff__tangent_ratio_streamlen");

    Add(uv.scale_multiplier, "uv__scale_multiplier");
    Add(uv.wsnap, "uv__wsnap");
  }

  // ===========================================================================

  void Options::set_default() {
    in.threshold               = 165;
    in.mask_factor             = 2.0f;
    tri.adaptive               = true;
    tri.estimate_max_area      = true;
    tri.max_area_coeff         = 0.25f;
    tri.narrow_band_radius     = 0.3f;
    tri.wdist_sigma            = 0.1f;
    ff.use_weighted_reg        = true;
    ff.wsmooth                 = 50.0f;
    ff.wtangent                = 1.0f;
    ff.wortho                  = 0.02f;
    ff.tangent_ratio_threshold = 0.48;
    ff.tangent_ratio_streamlen = 4.0;
    uv.scale_multiplier        = 1.0f;
    uv.wsnap                   = 25.0f;
  }

  // ===========================================================================

} // namespace IGSV

#include <fstream>

template void IGSV::Options::print<std::basic_ostream<char, std::char_traits<char>>>(
    std::basic_ostream<char, std::char_traits<char>>&) const;

template void IGSV::Options::print<std::basic_ofstream<char, std::char_traits<char>>>(
    std::basic_ofstream<char, std::char_traits<char>>&) const;
