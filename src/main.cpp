#include <igsv/common/VectorizationData.h>
#include <igsv/common/print_timings.h>

#ifdef IGSV_WITH_GUI
#include <igsv/gui2/Gui.h>
#endif

int main(int argc, char* argv[]) {

  IGSV::VectorizationData data;

  // parse options
  if (!data.opts.parse_args(argc, argv))
    return 0;

  // command line mode
  if (data.opts.misc.no_gui) {

    // compute from scratch and serialize
    if (data.compute_no_gui()) {
      IGSV::log_section("serialize");
      igl::serialize(data, "state", data.opts.in.filename_input + ".state");
      IGSV::log_message("state saved to " + data.opts.in.filename_input + ".state");
      return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
  }

  // check if gui mode available
#ifndef IGSV_WITH_GUI

  // gui was *not* compiled
  IGSV::log_error(
      "requested gui mode, but gui was *not* compiled\n"
      "       use the flag --no-gui or compile with IGSV_COMPILE_GUI on\n"
      "abort.");
  return EXIT_FAILURE;

#else

  // gui is available

  if (data.opts.in.filename_state.empty()) {
    // compute from scratch
    data.compute_no_gui();

  } else {
    // deserialize from a saved state
    igl::deserialize(data, "state", data.opts.in.filename_state);

    if (!data.init_sketch())
      return EXIT_FAILURE;

    if (!data.init_mesh(false))
      return EXIT_FAILURE;

    if (!data.init_frame_field(false, false))
      return EXIT_FAILURE;

    if (!data.init_parametrization(false))
      return EXIT_FAILURE;

    if (!data.init_extraction())
      return EXIT_FAILURE;
  }

  IGSV::Gui2 gui2(data.input,                       //
                  data.detail_mask,                 //
                  data.opts.tri.narrow_band_radius, //
                  data.opts.uv.scale_multiplier,    //
                  data.opts.in.mask_factor,         //
                  data.V,                           //
                  data.UV,                          //
                  data.BC,                          //
                  data.X0_nonunit_comb,             //
                  data.X1_nonunit_comb,             //
                  data.F,                           //
                  data.FUV,                         //
                  data.F_isNarrow,                  //
                  data.QV,                          //
                  data.QE,                          //
                  data.CH                           //
  );

  gui2._callback = [&data]() -> bool {        //
    return data.init_parametrization(true) && //
           data.init_extraction();
  };

  gui2.update_structures();
  gui2.render();

  return EXIT_SUCCESS;

#endif
}
