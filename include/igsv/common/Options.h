//=============================================================================
//
//  CLASS : Options
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <boost/program_options.hpp>
#include <igl/serialize.h>
#include <string>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ======================================================

  class Options : public igl::Serializable {

  public:
    Options() { set_default(); }

    // parse command line arguments using boost::program_options
    bool parse_args(int argc, char* argv[]);

    // set default values
    void set_default();

    // serialization
    void InitSerialization() override;

    // output values to a stream
    template <typename StreamT>
    void print(StreamT&) const;

    // ===========================================================================

    struct {

      std::string filename_input  = "";
      std::string filename_output = "";
      std::string filename_mask   = "";
      std::string filename_state  = "";
      int threshold;
      float mask_factor;
    } in{};

    // ===========================================================================

    struct {
      bool adaptive;
      bool estimate_max_area;
      float max_area = 1.0f;
      float max_area_coeff;
      float narrow_band_radius;
      float wdist_sigma;
    } tri{};

    // ===========================================================================

    struct {
      bool use_weighted_reg;
      float wsmooth;
      float wtangent;
      float wortho;
      float tangent_ratio_threshold;
      float tangent_ratio_streamlen;
    } ff{};

    // ===========================================================================

    struct {
      float scale_multiplier;
      float wsnap;
    } uv{};

    // ===========================================================================

    struct {
      bool compute_from_scratch = true;
      bool no_gui               = false;
    } misc;

    // ===========================================================================
  };

} // namespace IGSV
