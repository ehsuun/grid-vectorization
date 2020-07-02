#include <CoMISo/Utils/StopWatch.hh>
#include <boost/format.hpp>
#include <igsv/igsv.h>

namespace IGSV {

  // ===========================================================================

  void VectorizationData::InitSerialization() {
    Add(opts, "opts");
    Add(mask_read_from_file, "mask_read_from_file");
    Add(V, "V");
    Add(F, "F");
    Add(UV, "UV");
    Add(FUV, "FUV");
    Add(X0_arg, "X0_arg");
    Add(X1_arg, "X1_arg");
    Add(X0_mag, "X0_mag");
    Add(X1_mag, "X1_mag");
    Add(E_periodJumps, "E_periodJumps");
    Add(F_labels, "F_labels");
  }

  // ===========================================================================

  bool VectorizationData::PreSerialization() const {
    // export_info();
    DEBUG_PRINT_FN("Serialize");
    return true;
  }

  // ===========================================================================

  void VectorizationData::PostSerialization() const { export_svg(); }

  // ===========================================================================

  bool VectorizationData::compute_no_gui() {

    if (!init_sketch())
      return false;

    if (!init_mesh(true))
      return false;

    if (!init_frame_field(true, true))
      return false;

    if (!init_parametrization(true))
      return false;

    if (!init_extraction())
      return false;

    export_svg();

    double total_time = 0.0;
    std::for_each(timings.begin(), timings.end(),
                  [&total_time](const auto& label_time_pair) { total_time += label_time_pair.second; });

    log_message("\nall done!");
    log_time(total_time);

    return true;
  }

  // ===========================================================================

  bool VectorizationData::init_sketch() {

    log_section("read sketch");

    COMISO::StopWatch timer;
    timer.start();

    if (!read_sketch())
      return false;

    if (!process_sketch())
      return false;

    if (!compute_constraints())
      return false;

    timings["Sketch::init"] = timer.stop() / 1000.0;

    log_message((boost::format("# stroke pixels: %d") % n_dark).str());
    log_message((boost::format("avg stroke width: %.2f pixels") % sw_avg).str());
    log_time(timings["Sketch::init"]);

    return true;
  }

  // ===========================================================================

  bool VectorizationData::init_mesh(bool _compute_from_scratch) {

    log_section("generate mesh");

    COMISO::StopWatch timer;
    timer.start();

    if (opts.tri.estimate_max_area) {
      if (!estimate_max_area())
        return false;
    }

    log_message((boost::format("max tri area: %.2f") % opts.tri.max_area).str());

    if (_compute_from_scratch) {
      if (!triangulate())
        return false;

    } else {
      if (!process_mesh())
        return false;
    }

    if (!find_nearest())
      return false;

    if (!compute_narrow_band())
      return false;

    if (!compute_snapping_weights())
      return false;

    timings["Mesh::init"] = timer.stop() / 1000.0;

    log_message((boost::format("# verts: %d") % V_isNarrow.count()).str());
    log_message((boost::format("# faces: %d") % F_isNarrow.count()).str());
    log_time(timings["Mesh::init"]);

    return true;
  }

  // ===========================================================================

  bool VectorizationData::init_frame_field(bool _compute_from_scratch, bool _trace_streamlines) {

    log_section("compute frame field");

    COMISO::StopWatch timer_ff_all, timer;

    timer_ff_all.start();

    if (_compute_from_scratch) {
      timer.start();
      if (!compute_frame_field_narrow_band())
        return false;
      timings["FF::solve"] = timer.stop() / 1000.0;
    }

    log_message((boost::format("# iterations : %d") % n_ff_iters).str());
    log_message((boost::format("# small angle: %d of %d verts (%.1f%%)") //
                 % n_ff_degenerate                                       //
                 % V_isNarrow.count()                                    //
                 % (100. * n_ff_degenerate / (double)V_isNarrow.count()))
                    .str());

    timer.start();
    if (!comb_frame_field())
      return false;
    timings["FF::comb"] = timer.stop() / 1000.0;

    log_message("field combing: ok");

    if (_trace_streamlines) {
      timer.start();
      std::cout << "streamline tracing: ";
      if (!trace_streamlines())
        return false;
      timings["FF::trace"] = timer.stop() / 1000.0;
    }

    timer.start();
    if (!simplify_cut_graph_iterate())
      return false;
    timings["FF::simplify"] = timer.stop() / 1000.0;

    log_message("cut graph simplification: ok");

    log_time(timer_ff_all.stop() / 1000.0);

    return true;
  }

  // ===========================================================================

  bool VectorizationData::init_parametrization(bool _compute_from_scratch) {

    log_section("parametrize");

    log_message((boost::format("scale : %.1e") % opts.uv.scale_multiplier).str());
    log_message((boost::format("maskf : %.1e") % opts.in.mask_factor).str());

    // TODO : replace this, only need to re-compute the mask factor
    if (!process_sketch())
      return false;

    COMISO::StopWatch timer, timer_uv_all;

    timer_uv_all.start();

    timer.start();
    if (!compute_scale())
      return false;
    timings["UV:compute_scale"] = timer.stop() / 1000.0;

    if (_compute_from_scratch) {
      double _solve_time(0);
      timer.start();
      if (!parametrize(_solve_time))
        return false;
      timings["UV::mi_solve"] = timer.stop() / 1000.0;
    }

    log_time(timer_uv_all.stop() / 1000.0);

    return store_orientation();
  }

  // ===========================================================================

  bool VectorizationData::init_extraction() {

    log_section("extract");

    COMISO::StopWatch timer;
    timer.start();

    clear_extraction();

    if (!sanitize_parametrization())
      return false;

    if (!precompute_transition_fn())
      return false;

    // TODO : disable in command line mode
    if (!extract_line_map())
      return false;

    if (!extract_qvertices())
      return false;

    if (!extract_qports())
      return false;

    if (!extract_qedges())
      return false;

    if (!assign_pixels_to_qedges())
      return true;

    if (!init_base_labels(false))
      return true;

    // for (int iter = 0; iter < 1; ++iter) { ... }
    if (!init_base_labels(true))
      return true;

    if (!init_chains())
      return true;

    if (!parametrize_chains())
      return true;

    if (!assign_samples_to_chains())
      return true;

    if (!remove_short_chains())
      return true;

    // need to re-assign samples in case some chains were removed
    if (!assign_samples_to_chains())
      return true;

    if (!determine_chain_adjacency())
      return true;

    if (!fit_chains())
      return true;

    timings["EX::everything"] = timer.stop() / 1000.0;

    log_time(timings["EX::everything"]);

    return true;
  }

  // ===========================================================================

  bool VectorizationData::read_sketch() {
    DEBUG_PRINT_FN("Read sketch");
    return IGSV::read_sketch(opts.in.filename_input, opts.in.filename_mask, opts.in.mask_factor, input, detail_mask,
                             mask_read_from_file);
  }

  // ===========================================================================

  bool VectorizationData::process_sketch() {
    DEBUG_PRINT_FN("Process sketch");
    return IGSV::process_sketch(input, detail_mask, opts.in.threshold, grey, bw, dt, idt_01, sw, sw_01, sw_max, sw_avg,
                                sw_stddev);
  }

  // ===========================================================================

  bool VectorizationData::compute_constraints() {
    DEBUG_PRINT_FN("Compute constraints");
    return IGSV::compute_constraints(grey, bw, sw_avg, grey_blurred, n_dark, pix_weight, pix2band, PX_I, PX_J, PX_XY,
                                     PX_SobelTangent, PX_SobelWeight);
  }

  // ===========================================================================

  bool VectorizationData::estimate_max_area() {
    DEBUG_PRINT_FN("Estimate max area");
    /*
    -- T : equilateral triangle with side a
    -- goal: sw_avg = diameter of circumcircle of T (= d) with max_area
    -- r = a / √3        ==> a = r * √3
    -- d = 2 * a / √3    ==> a = d * √3 / 2
    -- d = sw_avg
    -- Area(T) = a² * √3 / 4
               = (d * √3 / 2)² * √3 / 4
               = d² * 3√3 / 16
              = sw_avg² * 3√3 / 16
    */
    opts.tri.max_area = std::pow(sw_avg - sw_stddev, 2) * 3.0 * std::sqrt(3.0) / 16.0;
    // scale down as needed to increase the density
    opts.tri.max_area *= opts.tri.max_area_coeff;
    opts.tri.estimate_max_area = false;
    DEBUG_PRINT_INFO("  max_area = %0.4f", opts.tri.max_area);
    return true;
  }

  // ===========================================================================

  bool VectorizationData::triangulate() {
    // GENERATE the mesh, then process
    DEBUG_PRINT_FN("Triangulate");
    return IGSV::triangulate(bw, dt, opts.tri.max_area, opts.tri.adaptive, opts.tri.narrow_band_radius, sw_avg, V, F,
                             UV_generated, BC, TT, TTi, VT, VTi, F_isBlack);
  }

  // ===========================================================================

  bool VectorizationData::process_mesh() {
    // DO NOT GENERATE the mesh, only process
    DEBUG_PRINT_FN("Process mesh");
    return IGSV::process_mesh(bw, V, F, UV_generated, BC, TT, TTi, VT, VTi, F_isBlack);
  }

  // ===========================================================================

  bool VectorizationData::find_nearest() {
    DEBUG_PRINT_FN("Find nearest");
    return IGSV::find_nearest(V, F, VT, PX_XY, tri, PX_nearestID, PX_nearestXY, PX_fid, PX_bc);
  }

  // ===========================================================================

  bool VectorizationData::compute_narrow_band() {
    DEBUG_PRINT_FN("Compute narrow band");
    return IGSV::compute_narrow_band(V, F, TT, TTi, VT, VTi, BC, bw, dt, detail_mask, opts.in.mask_factor,
                                     opts.tri.narrow_band_radius, sw_avg, V_isBlack, V_isNarrow, V_isNarrowBoundary,
                                     F_isNarrow, E_isNarrow, E_isNarrowBoundary);
  }

  // ===========================================================================

  bool VectorizationData::compute_snapping_weights(bool print_message) {
    if (print_message)
      DEBUG_PRINT_FN("Compute snapping weights");
    return IGSV::compute_snapping_weights(V, idt_01, V_snappingWeights);
  }

  // ===========================================================================

  bool VectorizationData::compute_frame_field_narrow_band() {

    DEBUG_PRINT_FN("Compute frame field (narrow band only)");

    // temps vars
    int vs, fs;

    //// index maps between base mesh and narrow and sub-mesh
    std::map<int, int> V_base_to_sub, V_sub_to_base;
    std::map<int, int> F_base_to_sub, F_sub_to_base;
    std::set<int> used_base_indices;

    // -- temp faces
    Eigen::MatrixXi F0_temp;
    F0_temp.resize(F.rows(), 3);
    fs = 0;
    for (int fb = 0; fb < F.rows(); ++fb)
      if (F_isNarrow(fb)) {

        F0_temp.row(fs) << F.row(fb);

        used_base_indices.insert(F(fb, 0));
        used_base_indices.insert(F(fb, 1));
        used_base_indices.insert(F(fb, 2));

        F_base_to_sub.insert(std::make_pair(fb, fs));
        F_sub_to_base.insert(std::make_pair(fs, fb));

        fs++;
      }
    F0_temp.conservativeResize(fs, 3);

    //// narrow band sub mesh
    Eigen::MatrixXi F0;
    Eigen::MatrixXd V0;

    // -- vertices
    V0.resize(V.rows(), 3);
    vs = 0;
    for (auto vb : used_base_indices) {
      V_base_to_sub.insert(std::make_pair(vb, vs));
      V_sub_to_base.insert(std::make_pair(vs, vb));
      V0.row(vs) = V.row(vb);
      vs++;
    }
    V0.conservativeResize(vs, 3);

    // -- faces
    F0.setZero(F0_temp.rows(), 3);
    for (int f = 0; f < F0_temp.rows(); ++f)
      for (int k = 0; k < 3; ++k)
        F0(f, k) = V_base_to_sub[F0_temp(f, k)];

    //// convert nearest vertex indices
    Eigen::VectorXi nearest;
    nearest.setZero(PX_nearestID.rows(), 1);
    for (int i = 0; i < PX_nearestID.rows(); ++i)
      nearest(i) = V_base_to_sub[PX_nearestID(i)];

    //// compute frame field
    Eigen::VectorXd V_arg0, V_mag0, V_arg1, V_mag1;
    Eigen::VectorXd F_arg0, F_mag0, F_arg1, F_mag1;
    IGSV::FrameFieldGenerator ffgenerator(V0,
                                          F0,              //
                                          PX_SobelTangent, //
                                          PX_SobelWeight,  //
                                          nearest,         //
                                          V_arg0,          //
                                          V_mag0,          //
                                          V_arg1,          //
                                          V_mag1,          //
                                          F_arg0,          //
                                          F_mag0,          //
                                          F_arg1,          //
                                          F_mag1           //
    );

    log_message((boost::format("w_smooth  : %.2f") % opts.ff.wsmooth).str());
    log_message((boost::format("w_tangent : %.2f") % opts.ff.wtangent).str());
    log_message((boost::format("w_regular : %.2f") % opts.ff.wortho).str());
    log_message(std::string("adaptive? : ") + std::string(opts.ff.use_weighted_reg ? "yes" : "no"));

    ffgenerator.set_wsmooth(opts.ff.wsmooth);
    ffgenerator.set_wtangent(opts.ff.wtangent);
    ffgenerator.set_wortho(opts.ff.wortho);
    ffgenerator.set_weighted_reg(opts.ff.use_weighted_reg);

    n_ff_iters      = 0;
    n_ff_degenerate = 0;

    log_message("solve ...");

    //// solve
    if (!ffgenerator.compute())
      return false;

    n_ff_iters      = ffgenerator._n_iters;
    n_ff_degenerate = ffgenerator._n_degenerate;

    //// copy
    X0_arg_vtx.setZero(V.rows(), 1);
    X0_arg_vtx.setZero(V.rows(), 1);
    X0_mag_vtx.setZero(V.rows(), 1);
    X1_arg_vtx.setZero(V.rows(), 1);
    X1_mag_vtx.setZero(V.rows(), 1);
    X0_arg.setZero(F.rows(), 1);
    X0_mag.setZero(F.rows(), 1);
    X1_arg.setZero(F.rows(), 1);
    X1_mag.setZero(F.rows(), 1);

    for (vs = 0; vs < V0.rows(); ++vs) {
      const int vb   = V_sub_to_base[vs];
      X0_arg_vtx(vb) = V_arg0(vs);
      X0_mag_vtx(vb) = V_mag0(vs);
      X1_arg_vtx(vb) = V_arg1(vs);
      X1_mag_vtx(vb) = V_mag1(vs);
    }

    for (fs = 0; fs < F0.rows(); ++fs) {
      const int fb = F_sub_to_base[fs];
      X0_arg(fb)   = F_arg0(fs);
      X0_mag(fb)   = F_mag0(fs);
      X1_arg(fb)   = F_arg1(fs);
      X1_mag(fb)   = F_mag1(fs);
    }

    return true;
  }

  // ===========================================================================

  bool VectorizationData::comb_frame_field() {
    DEBUG_PRINT_FN("Comb frame field");
    return IGSV::custom_comb_frame_field(
        V, F, TT, TTi, VT, VTi, V_isNarrow, E_isNarrow, F_isNarrow, X0_arg, X1_arg, X0_mag, X1_mag, V_singularityIdx,
        E_periodJumps, E_customCut, TT_dualConnections, BIS0_unit, BIS1_unit, BIS0_unit_comb, BIS1_unit_comb, X0_unit,
        X1_unit, X0_nonunit, X1_nonunit, X0_unit_comb, X1_unit_comb, X0_nonunit_comb, X1_nonunit_comb, F_frameField);
  }

  // ===========================================================================

  bool VectorizationData::trace_streamlines() {
    DEBUG_PRINT_FN("Trace streamlines");
    return IGSV::trace_streamlines(TT, F_isNarrow, F_isBlack, BC, F_frameField, E_periodJumps, tri, bw, sw_avg,
                                   opts.uv.scale_multiplier, opts.ff.tangent_ratio_threshold,
                                   opts.ff.tangent_ratio_streamlen, F_labels, FaceDataList);
  }

  // ===========================================================================

  bool VectorizationData::simplify_cut_graph_iterate() {
    DEBUG_PRINT_FN("Simplify cut graph (all)");
    return IGSV::simplify_cut_graph_iterate(V, F, TT, TTi, V_singularityIdx, V_isNarrowBoundary, E_customCut);
  }

  // ===========================================================================

  bool VectorizationData::compute_scale() {
    DEBUG_PRINT_FN("Compute scale");
    return IGSV::compute_scale(BC, detail_mask, opts.in.mask_factor, sw_avg, opts.uv.scale_multiplier, F_uvScale);
  }

  // ===========================================================================

  bool VectorizationData::parametrize(double& time) {
    DEBUG_PRINT_FN("Parametrize");
    IGSV::IGMSolver igmsolver(V, X0_unit_comb, X1_unit_comb, F, TT, TTi, VT, VTi, V_isBlack, V_isNarrowBoundary,
                              V_snappingWeights, V_singularityIdx, E_periodJumps, E_isNarrow, E_isNarrowBoundary,
                              E_customCut, F_isNarrow, F_labels, F_uvScale, UV, FUV, FK_snapping_constraints,
                              I_snapping_connections);

    igmsolver.set_snapping_weight((double)opts.uv.wsnap);

    return igmsolver.solve(time);
  }

  // ===========================================================================

  bool VectorizationData::store_orientation() {
    DEBUG_PRINT_FN("Store orientation");
    return IGSV::store_orientation(UV, FUV, F_isNarrow, F_uvOrientation);
  }

  // ===========================================================================

  bool VectorizationData::sanitize_parametrization() {
    DEBUG_PRINT_FN("Sanitize parametrization");
    return IGSV::sanitize_parametrization(UV);
  }

  // ===========================================================================

  bool VectorizationData::precompute_transition_fn() {
    DEBUG_PRINT_FN("Precompute trans fn");
    return IGSV::precompute_transition_fn(UV, FUV, TT, TTi, E_periodJumps, R_transfn, T_transfn, T_integer_transfn,
                                          T_decimal_transfn, flip_transfn);
  }

  // ===========================================================================

  bool VectorizationData::extract_line_map() {
    DEBUG_PRINT_FN("Extract line map");
    return igm.compute(V, F, UV, FUV, F_isNarrow);
  }

  // ===========================================================================

  bool VectorizationData::extract_qvertices() {
    DEBUG_PRINT_FN("Extract Q-Vertices");
    return IGSV::extract_qvertices(V, F, UV, FUV, F_isNarrow, QV);
  }

  // ===========================================================================

  bool VectorizationData::extract_qports() {
    DEBUG_PRINT_FN("Extract Q-Ports");
    return IGSV::extract_qports(V, F, UV, FUV, VT, VTi, F_isNarrow, QV, QP, F_qports, QP_QP_connections);
  }

  // ===========================================================================

  bool VectorizationData::extract_qedges() {
    DEBUG_PRINT_FN("Extract Q-Edges");
    return IGSV::extract_qedges(UV, FUV, TT, TTi, E_periodJumps, R_transfn, T_integer_transfn, E_isNarrow, F_isNarrow,
                                F_labels, F_qports, QV, QP, QE, QP_QP_connections, QVPair_to_QE);
  }

  // ===========================================================================

  bool VectorizationData::assign_pixels_to_qedges() {
    DEBUG_PRINT_FN("Assign pixels to Q-Edges");
    return IGSV::assign_pixels_to_qedges(grey, sw_01, PX_XY, PX_I, PX_J, PX_fid, X0_unit_comb, X1_unit_comb, QV, QP, QE,
                                         sw_avg, PX);
  }

  // ===========================================================================

  bool VectorizationData::init_base_labels(bool split_at_nodes) {
    DEBUG_PRINT_FN("Init base labels");
    return IGSV::init_base_labels(QP_QP_connections, split_at_nodes, QV, QP, QE);
  }

  // ===========================================================================

  bool VectorizationData::init_chains() {
    DEBUG_PRINT_FN("Init chains");
    return IGSV::init_chains(QV, QE, CH);
  }

  // ===========================================================================

  bool VectorizationData::parametrize_chains() {
    DEBUG_PRINT_FN("Parametrize chains");
    return IGSV::parametrize_chains(QV, QE, CH);
  }

  // ===========================================================================

  bool VectorizationData::assign_samples_to_chains() {
    DEBUG_PRINT_FN("Assign samples to chains");
    return IGSV::assign_samples_to_chains(QE, PX, CH);
  }

  // ===========================================================================

  bool VectorizationData::remove_short_chains() {
    DEBUG_PRINT_FN("Remove short and split long chains");
    return IGSV::remove_short_chains(QV, QE, CH);
  }

  // ===========================================================================

  bool VectorizationData::determine_chain_adjacency() {
    DEBUG_PRINT_FN("Determine chain adjacency");
    return IGSV::determine_chain_adjacency(QE, CH, CH_CH_connections, CH_CH_connections_enabled);
  }

  // ===========================================================================

  bool VectorizationData::fit_chains() {
    DEBUG_PRINT_FN("Fit chains");
    return IGSV::fit_chains(QV, QVPair_to_QE, QE, CH, sw_avg, opts.tri.narrow_band_radius);
  }

  // ===========================================================================

  bool VectorizationData::export_svg() const {
    // DEBUG_PRINT_FN("Export SVG");

    log_section("export svg");
    return IGSV::export_svg(opts.in.filename_output, CH, input.cols, input.rows);
  }

  // ===========================================================================

  // bool VectorizationData::export_info() const {
  //
  //   DEBUG_PRINT_FN("Export info");
  //
  //   ////// write a human-readable state file
  //   std::ofstream s;
  //   s.open(get_filename(FILE_FLAG::INFO));
  //   s << get_filename(FILE_FLAG::STATE) << std::endl << std::endl;
  //   opts.print(s);
  //   print_timings(s);
  //   s.close();
  //
  //   ////// print a message
  //   std::cout << "    saved " << get_filename(FILE_FLAG::INFO) << std::endl;
  //   opts.print(std::cout);
  //   print_timings(std::cout);
  //
  //   return true;
  // }

  // ===========================================================================

  template <typename StreamT>
  void VectorizationData::print_timings(StreamT& s) const {
    IGSV::print_timings(timings, s);
  }

  // ===========================================================================

  void VectorizationData::clear() {

    timings.clear();

    mask_read_from_file = false;

    sw_avg    = 0;
    sw_max    = 0;
    sw_stddev = 0;
    n_dark    = 0;

    pix_weight.resize(0, 1);
    pix2band.resize(0, 1);
    PX_I.resize(0, 1);
    PX_J.resize(0, 1);
    PX_XY.resize(0, 1);
    PX_fid.resize(0, 1);
    PX_bc.resize(0, 3);
    PX_SobelTangent.resize(0, 1);
    PX_SobelWeight.resize(0, 1);

    V.resize(0, 3);
    F.resize(0, 3);
    BC.resize(0, 3);
    UV_generated.resize(0, 2);
    UV.resize(0, 2);
    FUV.resize(0, 3);
    TT.resize(0, 3);
    TTi.resize(0, 3);
    VT.clear();
    VTi.clear();
    tri.clear();

    PX_nearestXY.resize(0, 1);
    PX_nearestID.resize(0, 1);

    V_snappingWeights.resize(0, 1);
    V_singularityIdx.resize(0, 1);
    V_isBlack.resize(0, 1);
    V_isNarrow.resize(0, 1);
    V_isNarrowBoundary.resize(0, 1);
    V_isNonManifold.resize(0, 1);

    E_isNarrow.resize(0, 3);
    E_isNarrowBoundary.resize(0, 3);
    E_customCut.resize(0, 3);
    E_periodJumps.resize(0, 3);

    F_isBlack.resize(0, 1);
    F_isNarrow.resize(0, 1);
    F_labels.resize(0, 1);

    F_uvScale.resize(0, 1);

    F_frameField[0].resize(0, 1);
    F_frameField[1].resize(0, 1);
    F_frameField[2].resize(0, 1);
    F_frameField[3].resize(0, 1);

    FaceDataList.clear();
    TT_dualConnections.clear();

    FK_snapping_constraints.resize(0, 3);
    I_snapping_connections.resize(0, 4);

    X0_arg_vtx.resize(0, 1);
    X0_mag_vtx.resize(0, 1);
    X1_arg_vtx.resize(0, 1);
    X1_mag_vtx.resize(0, 1);

    X0_arg.resize(0, 1);
    X0_mag.resize(0, 1);
    X1_arg.resize(0, 1);
    X1_mag.resize(0, 1);

    X0_unit.resize(0, 3);
    X1_unit.resize(0, 3);
    X0_nonunit.resize(0, 3);
    X1_nonunit.resize(0, 3);
    X0_unit_comb.resize(0, 3);
    X1_unit_comb.resize(0, 3);
    X0_nonunit_comb.resize(0, 3);
    X1_nonunit_comb.resize(0, 3);

    F_uvOrientation.resize(0, 1);

    clear_extraction();
  }

  // ===========================================================================

  void VectorizationData::clear_extraction() {

    R_transfn.resize(0, 3);
    T_transfn.resize(0, 3);
    T_integer_transfn.resize(0, 3);
    T_decimal_transfn.resize(0, 3);
    flip_transfn.resize(0, 3);

    igm.clear();

    QV.clear();
    QP.clear();
    QE.clear();
    PX.clear();
    CH.clear();

    F_qports.clear();
    QP_QP_connections.clear();
    QVPair_to_QE.clear();
    CH_CH_connections.clear();
    CH_CH_connections_enabled.clear();
  }

  // ===========================================================================

} // namespace IGSV