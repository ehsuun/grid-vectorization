#include <igl/colormap.h>
#include <igsv/common/CustomColors.h>
#include <igsv/common/VectorizationData.h>
#include <igsv/extraction_entities/Bezier.h>
#include <igsv/gui/Gui.h>
#include <igsv/parametrization/pix_coords_from_xy.h>
#include <imgui/imgui.h>

void IGSV::IGSVGui::store_entities() {

  // ===========================================================================

  if (data.F_uvOrientation.rows() > 0) {

    F__uvFlipped.clr.resize(data.F_uvOrientation.rows(), 3);
    F__uvFlipped.clr.fill(1.0);

    for (int f = 0; f < data.F_uvOrientation.rows(); ++f)
      if (data.F_isNarrow(f))
        switch (data.F_uvOrientation(f)) {

        case -1: // flipped
          F__uvFlipped.clr.row(f) << CustomColors::Coral;
          break;

        case 0: // degenerate
          F__uvFlipped.clr.row(f) << CustomColors::Yellow;
          break;

        default: // regular
          F__uvFlipped.clr.row(f) << CustomColors::PaleGreen;
        }

  } else {

    F__uvFlipped.clr = CustomColors::Black;
  }

  // ===========================================================================

  if (data.X0_arg.rows() > 0) {
    F__ffAngles.clr.resize(data.F.rows(), 3);
    F__ffAngles.clr.fill(1.0);
    for (int f = 0; f < data.F.rows(); ++f)
      if (data.F_isNarrow(f)) {
        // compute the difference of complex arguments : be careful about periodicity!
        double ua   = data.X0_arg(f);
        double va   = data.X1_arg(f);
        double diff = std::abs(ua - va);
        while (diff > +M_PI)
          diff -= 2 * M_PI;
        while (diff < -M_PI)
          diff += 2 * M_PI;
        if (diff < 0.35 /* ~ 20 degrees */)
          F__ffAngles.clr.row(f) << CustomColors::Orange;
      }
  }

  // ===========================================================================

  {
    L__faces.pos.resize(data.F.rows(), 3);
    L__faces.clr = CustomColors::Black;
    L__faces.str.clear();
    L__faces.str.resize(data.F.rows());

    std::stringstream ss;
    int i, j;
    for (int f = 0; f < data.F.rows(); ++f) {
      pix_coords_from_xy(data.input.cols, data.input.rows, data.BC(f, 0), data.BC(f, 1), i, j);
      ss.str(std::string());
      ss << "d=" << (int)data.dt.at<unsigned char>(i, j);
      L__faces.str[f] = ss.str();
      L__faces.pos.row(f) << data.BC.row(f);
    }
  }

  // ===========================================================================

  if (data.X0_arg_vtx.rows() > 0) {
    L__ffAngles.pos.resize(data.V.rows(), 3);
    L__ffAngles.clr.resize(data.V.rows(), 3);
    L__ffAngles.clr.fill(0.0);
    L__ffAngles.str.clear();
    L__ffAngles.str.resize(data.V.rows());
    for (int v = 0; v < data.V.rows(); ++v) {
      double ua   = data.X0_arg_vtx(v);
      double va   = data.X1_arg_vtx(v);
      double diff = std::abs(ua - va);
      while (diff > +M_PI)
        diff -= 2 * M_PI;
      while (diff < -M_PI)
        diff += 2 * M_PI;
      std::stringstream ss;
      ss << "a=" << std::round(diff / M_PI * 180.0);
      L__ffAngles.pos.row(v) << data.V.row(v);
      L__ffAngles.str[v] = ss.str();
      if (diff < 0.35 /* ~ 20 degrees */) {
        L__ffAngles.clr.row(v) << CustomColors::Red;
      }
    }
  }

  // ===========================================================================

  {
    F__narrowBand.clr.resize(data.F.rows(), 3);
    F__narrowBand.clr.fill(0.5);
    for (int f = 0; f < data.F.rows(); ++f)
      if (data.F_isNarrow(f))
        F__narrowBand.clr.row(f).fill(1.0);
  }

  // ===========================================================================

  {
    // V__SobelTangents.clr = CustomColors::MediumVioletRed;
    V__SobelTangents.clr = CustomColors::SeaGreen;
    V__SobelTangents.pos.setZero(data.n_dark, 3);
    V__SobelTangents.pos.col(0) << data.PX_XY.real();
    V__SobelTangents.pos.col(1) << data.PX_XY.imag();
    V__SobelTangents.vec.setZero(data.n_dark, 3);
    V__SobelTangents.vec.col(0) << data.PX_SobelWeight.cwiseMax(0.1).array() * data.PX_SobelTangent.real().array();
    V__SobelTangents.vec.col(1) << data.PX_SobelWeight.cwiseMax(0.1).array() * data.PX_SobelTangent.imag().array();

    V__SobelNormals.clr = CustomColors::Pink;
    V__SobelNormals.pos = V__SobelTangents.pos;
    V__SobelNormals.vec.setZero(data.n_dark, 3);
    V__SobelNormals.vec.col(0) = -V__SobelTangents.vec.col(1);
    V__SobelNormals.vec.col(1) = V__SobelTangents.vec.col(0);
  }

  // ===========================================================================

  {
    E__blackPixNearestVertex.pos0.setZero(data.n_dark, 3);
    E__blackPixNearestVertex.pos0.col(0) = data.PX_XY.real();
    E__blackPixNearestVertex.pos0.col(1) = data.PX_XY.imag();
    E__blackPixNearestVertex.pos1.setZero(data.n_dark, 3);
    E__blackPixNearestVertex.pos1.col(0) = data.PX_nearestXY.real();
    E__blackPixNearestVertex.pos1.col(1) = data.PX_nearestXY.imag();
    E__blackPixNearestVertex.clr         = CustomColors::Black;
  }

  // ===========================================================================

  // {
  //   V__ffVector0.clr = CustomColors::U();
  //   V__ffVector1.clr = CustomColors::V();
  //
  //   V__ffVector0.pos.resize(data.F.rows() * 2, 3);
  //   V__ffVector1.pos.resize(data.F.rows() * 2, 3);
  //   V__ffVector0.vec.resize(data.F.rows() * 2, 3);
  //   V__ffVector1.vec.resize(data.F.rows() * 2, 3);
  //
  //   int row = 0;
  //   for (int f = 0; f < data.F.rows(); ++f)
  //     if (data.F_isNarrow(f)) {
  //
  //       V__ffVector0.pos.row(row) << data.BC(f, 0), data.BC(f, 1), 0.0;
  //       V__ffVector1.pos.row(row) << data.BC(f, 0), data.BC(f, 1), 0.0;
  //
  //       if (show_bisector_field) {
  //         V__ffVector0.vec.row(row) << data.BIS0_unit_comb(f, 0), data.BIS0_unit_comb(f, 1), 0.0;
  //         V__ffVector1.vec.row(row) << data.BIS1_unit_comb(f, 0), data.BIS1_unit_comb(f, 1), 0.0;
  //
  //       } else {
  //         V__ffVector0.vec.row(row) << data.X0_nonunit_comb(f, 0), data.X0_nonunit_comb(f, 1), 0.0;
  //         V__ffVector1.vec.row(row) << data.X1_nonunit_comb(f, 0), data.X1_nonunit_comb(f, 1), 0.0;
  //       }
  //       row++;
  //     }
  //   V__ffVector0.pos.conservativeResize(row, 3);
  //   V__ffVector1.pos.conservativeResize(row, 3);
  //   V__ffVector0.vec.conservativeResize(row, 3);
  //   V__ffVector1.vec.conservativeResize(row, 3);
  // }

  {
    // V__ffVector0.clr = CustomColors::Orange;
    // V__ffVector1.clr = CustomColors::Orange;
    V__ffVector0.clr = CustomColors::U();
    V__ffVector1.clr = CustomColors::V();
    V__ffVector0.pos.resize(data.F.rows() * 4, 3);
    V__ffVector1.pos.resize(data.F.rows() * 4, 3);
    V__ffVector0.vec.resize(data.F.rows() * 4, 3);
    V__ffVector1.vec.resize(data.F.rows() * 4, 3);
    int row = 0;
    for (int f = 0; f < data.F.rows(); ++f)
      if (data.F_isNarrow(f)) {
        V__ffVector0.pos.row(row) << data.BC(f, 0), data.BC(f, 1), 0.0;
        V__ffVector1.pos.row(row) << data.BC(f, 0), data.BC(f, 1), 0.0;
        V__ffVector0.pos.row(row + 1) << data.BC(f, 0), data.BC(f, 1), 0.0;
        V__ffVector1.pos.row(row + 1) << data.BC(f, 0), data.BC(f, 1), 0.0;
        V__ffVector0.vec.row(row) << -data.X0_nonunit_comb(f, 0), -data.X0_nonunit_comb(f, 1), 0.0;
        V__ffVector1.vec.row(row) << -data.X1_nonunit_comb(f, 0), -data.X1_nonunit_comb(f, 1), 0.0;
        V__ffVector0.vec.row(row + 1) << +data.X0_nonunit_comb(f, 0), +data.X0_nonunit_comb(f, 1), 0.0;
        V__ffVector1.vec.row(row + 1) << +data.X1_nonunit_comb(f, 0), +data.X1_nonunit_comb(f, 1), 0.0;
        row += 2;
      }
    V__ffVector0.pos.conservativeResize(row, 3);
    V__ffVector1.pos.conservativeResize(row, 3);
    V__ffVector0.vec.conservativeResize(row, 3);
    V__ffVector1.vec.conservativeResize(row, 3);
  }

  // ===========================================================================

  if (data.FaceDataList.size() > 0) {
    int ne = 0;
    E__streamlines.pos0.resize(1e7, 3);
    E__streamlines.pos1.resize(1e7, 3);
    E__streamlines.clr.resize(1e7, 3);
    // E__streamlines.clr = CustomColors::Orange;
    P__streamlines.pos.resize(1e7, 3);
    P__streamlines.clr.resize(1e7, 3);

    int nl = 0;
    L__nTracedBlackPix.str.clear();
    L__nTracedBlackPix.pos.resize(data.F.rows(), 3);
    L__nTracedBlackPix.clr.resize(data.F.rows(), 3);

    for (const auto& fdata : data.FaceDataList) {
      // show labels for all faces
      std::stringstream ss;
      ss << "f: " << fdata.f << std::endl;
      ss << "  u: " << fdata.c[0] << std::endl;
      ss << "  v: " << fdata.c[1];
      L__nTracedBlackPix.pos.row(nl) << fdata.bx, fdata.by, 0.0;
      L__nTracedBlackPix.clr.row(nl) << CustomColors::Black;
      L__nTracedBlackPix.str.push_back(ss.str());
      nl++;

      if (selected_face > -1 & selected_face != fdata.f)
        continue; // only show streamlines for selected face

      for (int mm = 0; mm < 4; mm++) {
        for (int i = 0; i < fdata.samples[mm].rows() - 1; i++) {

          // E__streamlines.pos0.row(ne) << fdata.samples[mm].row(i);
          // E__streamlines.pos1.row(ne) << fdata.samples[mm].row(i + 1);
          // P__streamlines.pos.row(ne) << fdata.samples[mm].row(i + 1);
          // if (streamline_random_colors) {
          //   E__streamlines.clr.row(ne) << CustomColors::Palette(2 * fdata.f + (mm % 2));
          // } else {
          // if (mm % 2 == 0) {
          //   E__streamlines.clr.row(ne) << streamline_color1[0], streamline_color1[1], streamline_color1[2];
          //   P__streamlines.clr.row(ne) << streamline_color1[0], streamline_color1[1], streamline_color1[2];
          // } else {
          //   E__streamlines.clr.row(ne) << streamline_color2[0], streamline_color2[1], streamline_color2[2];
          //   P__streamlines.clr.row(ne) << streamline_color2[0], streamline_color2[1], streamline_color2[2];
          // }
          // }

          if ((mm % 2 == 0 && data.F_labels[fdata.f] == CORNER_TYPE_INTEGER_U) ||
              (mm % 2 == 1 && data.F_labels[fdata.f] == CORNER_TYPE_INTEGER_V)) {
            E__streamlines.pos0.row(ne) << fdata.samples[mm].row(i);
            E__streamlines.pos1.row(ne) << fdata.samples[mm].row(i + 1);
            E__streamlines.clr.row(ne) << CustomColors::Palette(2 * fdata.f + (mm % 2));
            P__streamlines.pos.row(ne) << fdata.samples[mm].row(i + 1);
            ne++;
          }
        }
      }
    }
    L__nTracedBlackPix.pos.conservativeResize(nl, 3);
    L__nTracedBlackPix.clr.conservativeResize(nl, 3);
    P__streamlines.pos.conservativeResize(ne, 3);
    P__streamlines.clr.conservativeResize(ne, 3);
    E__streamlines.pos0.conservativeResize(ne, 3);
    E__streamlines.pos1.conservativeResize(ne, 3);
    E__streamlines.clr.conservativeResize(ne, 3);
  }

  // ===========================================================================

  {
    P__tangentDir.pos.resize(data.F.rows(), 3);
    P__tangentDir.clr.resize(data.F.rows(), 3);
    int r = 0;
    for (int f = 0; f < data.F.rows(); ++f)
      if (data.F_isBlack(f)) {
        P__tangentDir.pos.row(r) << data.BC.row(f);

        switch (data.F_labels(f)) {

        case 1:
          P__tangentDir.clr.row(r) = CustomColors::U();
          break;

        case 2:
          P__tangentDir.clr.row(r) = CustomColors::V();
          break;

        default:
          P__tangentDir.clr.row(r) = CustomColors::Black;
        }

        r++;
      }
    P__tangentDir.pos.conservativeResize(r, 3);
    P__tangentDir.clr.conservativeResize(r, 3);
  }

  // ===========================================================================

  {
    E__periodJumps.pos0.resize(3 * data.F.rows(), 3);
    E__periodJumps.pos1.resize(3 * data.F.rows(), 3);
    E__periodJumps.clr.resize(3 * data.F.rows(), 3);
    int row = 0;
    for (int f = 0; f < data.F.rows(); f++)
      for (int k = 0; k < 3; k++)
        if (data.E_periodJumps(f, k) > 0) {
          E__periodJumps.pos0.row(row) << data.V.row(data.F(f, k));
          E__periodJumps.pos1.row(row) << data.V.row(data.F(f, (k + 1) % 3));
          E__periodJumps.clr.row(row++) << CustomColors::Palette(data.E_periodJumps(f, k));
        }
    E__periodJumps.pos0.conservativeResize(row, 3);
    E__periodJumps.pos1.conservativeResize(row, 3);
    E__periodJumps.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    P__singularities.pos.resize(data.V.rows(), 3);
    P__singularities.clr.resize(data.V.rows(), 3);
    int row = 0;
    for (unsigned i = 0; i < data.V.rows(); ++i)
      if (data.V_singularityIdx(i) > 0) {
        P__singularities.pos.row(row) << data.V.row(i);
        P__singularities.clr.row(row++) << CustomColors::MediumOrchid;
      }
    P__singularities.pos.conservativeResize(row, 3);
    P__singularities.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    E__customCuts.clr.resize(3 * data.F.rows(), 3);
    E__customCuts.pos0.resize(3 * data.F.rows(), 3);
    E__customCuts.pos1.resize(3 * data.F.rows(), 3);
    int row = 0;
    if (data.E_customCut.rows() == data.F.rows() && data.E_customCut.cols() == 3)
      for (int f = 0; f < data.F.rows(); ++f)
        for (int k = 0; k < 3; ++k)
          if (data.E_customCut(f, k) & data.E_isNarrow(f, k)) {
            E__customCuts.clr.row(row) << CustomColors::Magenta;
            E__customCuts.pos0.row(row) << data.V.row(data.F(f, k));
            E__customCuts.pos1.row(row++) << data.V.row(data.F(f, (k + 1) % 3));
          }
    E__customCuts.pos0.conservativeResize(row, 3);
    E__customCuts.pos1.conservativeResize(row, 3);

    // E__customCuts.pos0.conservativeResize(row + data.TT_dualConnections.size(), 3);
    // E__customCuts.pos1.conservativeResize(row + data.TT_dualConnections.size(), 3);
    // for (const auto& d_conn : data.TT_dualConnections) {
    //   E__customCuts.clr.row(row) << CustomColors::Aquamarine;
    //   E__customCuts.pos0.row(row)   = data.BC.row(d_conn.first);
    //   E__customCuts.pos1.row(row++) = data.BC.row(d_conn.second);
    // }
  }

  // ===========================================================================

  {
    if (data.FUV.rows() == data.F.rows()) {
      L__corner.clr.resize(3 * data.F.rows(), 3);
      L__corner.pos.resize(3 * data.F.rows(), 3);
      L__corner.str.clear();
      L__corner.str.reserve(3 * data.F.rows());
      int row        = 0;
      const double t = 0.6; // in [0,1] : controls label placement between vertex and triangle barycenter
      for (int f = 0; f < data.F.rows(); ++f)
        if (data.F_isNarrow(f))
          for (int k = 0; k < 3; ++k) {
            std::stringstream ss;
            // ss << "[" << f << "," << k << "] : "                                                   // face, edge
            // index
            //    << "(" << std::fixed << std::setprecision(3) << data.UV(data.FUV(f, k), 0)         // u coord
            //    << "," << std::fixed << std::setprecision(3) << data.UV(data.FUV(f, k), 1) << ")"; // v coord
            ss << data.FUV(f, k) << " (" << k << ")";
            L__corner.pos.row(row) << (1.0 - t) * data.BC.row(f) + t * data.V.row(data.F(f, k));
            L__corner.clr.row(row) << CustomColors::Palette(data.FUV(f, k));
            L__corner.str.push_back(ss.str());
            row++;
          }
      L__corner.pos.conservativeResize(row, 3);
      L__corner.clr.conservativeResize(row, 3);
    }
  }

  // ===========================================================================

  {
    E__narrowBandBoundary.clr = CustomColors::Black;
    E__narrowBandBoundary.pos0.resize(data.F.rows() * 3, 3);
    E__narrowBandBoundary.pos1.resize(data.F.rows() * 3, 3);
    int row = 0;
    for (int f = 0; f < data.F.rows(); ++f)
      for (int k = 0; k < 3; ++k)
        if (data.E_isNarrowBoundary(f, k)) {
          E__narrowBandBoundary.pos0.row(row) << data.V.row(data.F(f, k));
          E__narrowBandBoundary.pos1.row(row) << data.V.row(data.F(f, (k + 1) % 3));
          row++;
        }
    E__narrowBandBoundary.pos0.conservativeResize(row, 3);
    E__narrowBandBoundary.pos1.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    E__lineMapSegments.pos0.resize(1e6, 3);
    E__lineMapSegments.pos1.resize(1e6, 3);
    E__lineMapSegments.clr.resize(1e6, 3);
    int row = 0;
    for (int k = 0; k < 2; ++k) {
      std::vector<LineSegment>* list;
      Color3d color;
      if (k == 0) {
        list  = &(data.igm.segmentsU);
        color = CustomColors::V();
      } else {
        list  = &(data.igm.segmentsV);
        color = CustomColors::U();
      }
      for (const auto& s : (*list))
        if (s.numpts == 2) {
          E__lineMapSegments.pos0.row(row) << s.xyz[0].x(), s.xyz[0].y(), 0.0;
          E__lineMapSegments.pos1.row(row) << s.xyz[1].x(), s.xyz[1].y(), 0.0;
          E__lineMapSegments.clr.row(row++) << color;
        }
    }
    E__lineMapSegments.pos0.conservativeResize(row, 3);
    E__lineMapSegments.pos1.conservativeResize(row, 3);
    E__lineMapSegments.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    // P__qvertices.pos = data.QV_xy;
    // P__qvertices.clr = CustomColors::Black;
    P__qvertices.pos.resize(data.QV.size(), 3);
    P__qvertices.clr.resize(data.QV.size(), 3);
    for (int row = 0; row < data.QV.size(); ++row) {
      P__qvertices.pos.row(row) << data.QV[row].x, data.QV[row].y, 0.0;
      P__qvertices.clr.row(row) << CustomColors::Palette((int)data.QV[row].type);
    }
  }

  // ===========================================================================

  {
    E__qports.pos0.resize(data.QP.size(), 3);
    E__qports.pos1.resize(data.QP.size(), 3);
    E__qports.clr.resize(data.QP.size(), 3);
    int qvi;
    for (int qpi = 0; qpi < data.QP.size(); ++qpi) {
      qvi = data.QP[qpi].qvi;
      E__qports.pos0.row(qpi) << data.QV[qvi].x, data.QV[qvi].y, 0.0;
      E__qports.pos1.row(qpi) << data.QP[qpi].x_end, data.QP[qpi].y_end, 0.0;
      E__qports.clr.row(qpi) << CustomColors::Palette(data.QP[qpi].pj); // color by period jump
    }
  }

  // ===========================================================================

  {
    // std::cout << "    E__qedges" << std::endl;

    E__qedges.pos0.resize(1e6, 3);
    E__qedges.pos1.resize(1e6, 3);
    E__qedges.clr.resize(1e6, 3);

    int row = 0;
    int qv0, qv1, i0, i1;

    // TODO
    // TODO
    // TODO

    int chid = 0;
    double t, r, g, b;
    for (auto& chain : data.CH)
      if (!chain.is_removed || show_unused_qedges) { // not removed *or* we want to show everything
        for (auto qei : chain.qedges) {

          qv0 = data.QE[qei].qvi0;
          qv1 = data.QE[qei].qvi1;
          i0  = chain.vindex_map[qv0];
          i1  = chain.vindex_map[qv1];

          if (chain.is_removed) {
            E__qedges.clr.row(row) << CustomColors::Red;

          } else {
            switch (selected_qedge_color) {

            case QEDGE_COLOR_UNIFORM:
              E__qedges.clr.row(row) << CustomColors::MediumSeaGreen;
              break;

            case QEDGE_COLOR_CHAIN_ID:
              E__qedges.clr.row(row) << CustomColors::Palette(chid);
              break;

            case QEDGE_COLOR_CHAIN_T:
              t = 0.5 * (chain.t[i0] + chain.t[i1]);
              igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, t, r, g, b);
              E__qedges.clr.row(row) << r, g, b;
              break;

            default:
              std::cerr << "    WARNING: invalid qedge color scheme" << std::endl;
              continue;
            }
          }

          E__qedges.pos0.row(row) << data.QV[qv0].x, data.QV[qv0].y, 0.0;
          E__qedges.pos1.row(row) << data.QV[qv1].x, data.QV[qv1].y, 0.0;

          row++;
        }

        if (!chain.is_removed)
          chid++;
      }

    //// NON-TANGENT EDGES
    if (show_unused_qedges)
      for (int qei = 0; qei < data.QE.size(); ++qei) {
        if (!data.QE[qei].is_used) {
          qv0 = data.QE[qei].qvi0;
          qv1 = data.QE[qei].qvi1;
          E__qedges.pos0.row(row) << data.QV[qv0].x, data.QV[qv0].y, 0.0;
          E__qedges.pos1.row(row) << data.QV[qv1].x, data.QV[qv1].y, 0.0;
          E__qedges.clr.row(row) << CustomColors::Orange;
          row++;
        } else if (data.QE[qei].chain == -1) {
          qv0 = data.QE[qei].qvi0;
          qv1 = data.QE[qei].qvi1;
          E__qedges.pos0.row(row) << data.QV[qv0].x, data.QV[qv0].y, 0.0;
          E__qedges.pos1.row(row) << data.QV[qv1].x, data.QV[qv1].y, 0.0;
          E__qedges.clr.row(row) << CustomColors::Magenta;
          row++;
        }
      }

    //// COLOR BY CHAIN PARAMETRIZATION

    E__qedges.pos0.conservativeResize(row, 3);
    E__qedges.pos1.conservativeResize(row, 3);
    E__qedges.clr.conservativeResize(row, 3);

    std::cout << row << " q-edges" << std::endl;
  }

  // ===========================================================================

  {
    // std::cout << "    E__pixSampleConnections" << std::endl;

    E__pixSampleConnections.pos0.resize(data.PX.size(), 3);
    E__pixSampleConnections.pos1.setZero(data.PX.size(), 3);
    E__pixSampleConnections.clr = CustomColors::Magenta;
    Sample* sample;
    int row = 0;
    for (int pid = 0; pid < data.PX.size(); ++pid) {
      if (data.PX[pid].samples.empty())
        continue;
      sample = &(data.PX[pid].samples.front()); // get the closest sample
      E__pixSampleConnections.pos0.row(row) << sample->ex, sample->ey, 0.0;
      E__pixSampleConnections.pos1.row(row) << sample->px, sample->py, 0.0;
      row++;
    }
    E__pixSampleConnections.pos0.conservativeResize(row, 3);
    E__pixSampleConnections.pos1.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    // std::cout << "    L__qexData" << std::endl;

    L__qexData.pos.resize(1e6, 3);
    L__qexData.clr.resize(1e6, 3);
    L__qexData.str.clear();

    int row = 0;
    int qvi, qpi, qei, qv0, qv1;
    std::stringstream ss;

    //// qverts
    for (qvi = 0; qvi < data.QV.size(); ++qvi) {
      L__qexData.pos.row(row) << data.QV[qvi].x, data.QV[qvi].y, 0.0;
      L__qexData.clr.row(row) << CustomColors::Blue;

      ss.str(std::string());
      ss << "qv #" << qvi;
      L__qexData.str.push_back(ss.str());

      row++;
    }

    //// qports
    static const double _t = 0.2;
    for (qpi = 0; qpi < data.QP.size(); ++qpi) {
      qvi = data.QP[qpi].qvi;
      L__qexData.pos.row(row) <<                                //
          (1. - _t) * data.QV[qvi].x + _t * data.QP[qpi].x_end, //
          (1. - _t) * data.QV[qvi].y + _t * data.QP[qpi].y_end, //
          0.0;
      L__qexData.clr.row(row) << CustomColors::Black;

      ss.str(std::string());
      ss << "qp #" << qpi;
      L__qexData.str.push_back(ss.str());

      row++;
    }

    //// q-edges
    for (qei = 0; qei < data.QE.size(); ++qei) {
      qv0 = data.QE[qei].qvi0;
      qv1 = data.QE[qei].qvi1;
      L__qexData.pos.row(row) <<                   //
          0.5 * (data.QV[qv0].x + data.QV[qv1].x), // midpoint
          0.5 * (data.QV[qv0].y + data.QV[qv1].y), //
          0.0;
      L__qexData.clr.row(row) << CustomColors::Green;
      ss.str(std::string());

      ss << "qe #" << qei;
      L__qexData.str.push_back(ss.str());

      row++;
    }
    L__qexData.pos.conservativeResize(row, 3);
    L__qexData.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    // std::cout << "    E__bezierCurves" << std::endl;

    // #define N_BEZIER_SAMPLES 50
    //     E__bezierCurves.pos0.resize(data.CU.size() * N_BEZIER_SAMPLES, 3);
    //     E__bezierCurves.pos1.resize(data.CU.size() * N_BEZIER_SAMPLES, 3);
    //     E__bezierCurves.clr.resize(data.CU.size() * N_BEZIER_SAMPLES, 3);
    //     int row = 0;
    //     double t;
    //     const double dt = 1.0 / (N_BEZIER_SAMPLES - 1);
    //     std::complex<double> p;
    //     int cui = 0;
    //     for (const auto& curve : data.CU) {
    //       t = 0;
    //       p = Bezier::BezierPoint(curve.CP, t);
    //       for (int e = 0; e < N_BEZIER_SAMPLES; ++e) {
    //         E__bezierCurves.pos0.row(row) << p.real(), p.imag(), 0.0;
    //         t += dt;
    //         p = Bezier::BezierPoint(curve.CP, t);
    //         E__bezierCurves.pos1.row(row) << p.real(), p.imag(), 0.0;
    //         E__bezierCurves.clr.row(row) << CustomColors::Palette(cui);
    //         row++;
    //       }
    //       cui++;
    //     }
  }

  // ===========================================================================

  {
    // std::cout << "    E__bezierPolygons" << std::endl;

    E__bezierPolygons.pos0.resize(1e6, 3);
    E__bezierPolygons.pos1.resize(1e6, 3);
    E__bezierPolygons.clr.resize(1e6, 3);
    int row = 0;
    int cui = 0;
    for (const auto& _chain : data.CH) {
      // for (int k = 0; k < curve.CP.rows() - 1; ++k) {
      //   E__bezierPolygons.pos0.row(row) << curve.CP(k).real(), curve.CP(k).imag(), 0.0;
      //   E__bezierPolygons.pos1.row(row) << curve.CP(k + 1).real(), curve.CP(k + 1).imag(), 0.0;
      //   E__bezierPolygons.clr.row(row) << CustomColors::Palette(cui);
      //   row++;
      // }
      for (int i = 0; i < _chain.CP_spline.cols(); ++i)
        for (int k = 0; k < 3; ++k) {
          E__bezierPolygons.pos0.row(row) << _chain.CP_spline(k, i).real(), _chain.CP_spline(k, i).imag(), 0.0;
          E__bezierPolygons.pos1.row(row) << _chain.CP_spline(k + 1, i).real(), _chain.CP_spline(k + 1, i).imag(), 0.0;
          E__bezierPolygons.clr.row(row) << CustomColors::Palette(cui);
          row++;
        }
      ++cui;
    }
    E__bezierPolygons.pos0.conservativeResize(row, 3);
    E__bezierPolygons.pos1.conservativeResize(row, 3);
    E__bezierPolygons.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    // std::cout << "    E__splineCurves" << std::endl;

    const int n_points_per_segment = 20;

    E__splineCurves.pos0.resize(1e6, 3);
    E__splineCurves.pos1.resize(1e6, 3);
    E__splineCurves.clr.resize(1e6, 3);

    double t0, t1, t, dt, c0, c1, c2, c3, d0, d1, d2, d3, x, y;

    int row = 0;
    int cui = 0;

    for (const auto& _chain : data.CH) {
      for (int i = 0; i < _chain.CubicSpline.rows(); ++i) {

        t0 = _chain.CubicSpline(i, 0);
        t1 = _chain.CubicSpline(i, 1);

        c0 = _chain.CubicSpline(i, 2);
        c1 = _chain.CubicSpline(i, 3);
        c2 = _chain.CubicSpline(i, 4);
        c3 = _chain.CubicSpline(i, 5);

        d0 = _chain.CubicSpline(i, 6);
        d1 = _chain.CubicSpline(i, 7);
        d2 = _chain.CubicSpline(i, 8);
        d3 = _chain.CubicSpline(i, 9);

        dt = (t1 - t0) / ((double)n_points_per_segment - 1.);

        // std::cout << "segment : t0=" << t0 << ", t1=" << t1 << std::endl;

        t = 0;
        // std::cout << "  t=" << t << std::endl;
        x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
        y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;

        for (int k = 0; k < n_points_per_segment - 1; ++k) {
          E__splineCurves.pos0.row(row) << x, y, 0.0;

          t += dt;
          // std::cout << "  t=" << t << std::endl;
          x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
          y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;

          E__splineCurves.pos1.row(row) << x, y, 0.0;
          E__splineCurves.clr.row(row) << CustomColors::Palette(cui);
          row++;
        }
      }
      cui++;
      // break;
    }
    E__splineCurves.pos0.conservativeResize(row, 3);
    E__splineCurves.pos1.conservativeResize(row, 3);
    E__splineCurves.clr.conservativeResize(row, 3);
  }

  // ===========================================================================

  {
    //
    E__snappingConnections.pos0.resize(data.I_snapping_connections.rows(), 3);
    E__snappingConnections.pos1.resize(data.I_snapping_connections.rows(), 3);
    E__snappingConnections.clr.resize(data.I_snapping_connections.rows(), 3);

    int f, ci, cj, vi, vj, type;
    const double t = 0.75;
    for (int row = 0; row < data.I_snapping_connections.rows(); ++row) {
      f    = data.I_snapping_connections(row, 0);
      ci   = data.I_snapping_connections(row, 1);
      cj   = data.I_snapping_connections(row, 2);
      type = data.I_snapping_connections(row, 3);
      vi   = data.F(f, ci);
      vj   = data.F(f, cj);

      E__snappingConnections.pos0.row(row) << t * data.V.row(vi) + (1. - t) * data.BC.row(f);
      E__snappingConnections.pos1.row(row) << t * data.V.row(vj) + (1. - t) * data.BC.row(f);
      E__snappingConnections.clr.row(row) << (type == 0 ? CustomColors::U() : CustomColors::V());
    }
  }

  // ===========================================================================

  {
    P__snappingConstraints.pos.resize(data.FK_snapping_constraints.rows() * 3, 3); //
    P__snappingConstraints.clr.resize(data.FK_snapping_constraints.rows() * 3, 3); //

    const double t = 0.75;

    int row = 0;

    for (int f = 0; f < data.FK_snapping_constraints.rows(); ++f)
      for (int k = 0; k < 3; ++k) {

        switch (data.FK_snapping_constraints(f, k)) {

        case 0:
          continue;

        case 1:
          P__snappingConstraints.clr.row(row) << CustomColors::U();
          break;

        case 2:
          P__snappingConstraints.clr.row(row) << CustomColors::V();
          break;

        case 3:
          P__snappingConstraints.clr.row(row) << CustomColors::Magenta;
          break;

        default:
          std::cerr << "WARNING: unexpected value in FK_snapping_constraints" << std::endl;
          continue;
        }

        P__snappingConstraints.pos.row(row) << t * data.V.row(data.F(f, k)) + (1. - t) * data.BC.row(f);
        row++;
      }

    P__snappingConstraints.pos.conservativeResize(row, 3); //
    P__snappingConstraints.clr.conservativeResize(row, 3); //
  }

  // ===========================================================================

  {
    L__pixels.pos.resize(data.PX_fid.rows(), 3);
    L__pixels.clr = CustomColors::Red;
    L__pixels.str.clear();
    L__pixels.str.resize(data.PX_fid.rows());
    int f = 0;
    std::stringstream ss;
    for (int k = 0; k < data.PX_fid.rows(); ++k) {

      f = data.PX_fid(k);

      ss.str(std::string());
      ss << "i=" << data.PX_I(k) << std::endl //
         << "j=" << data.PX_J(k)
         << std::endl
         // << "idt= " << (int)(data.idt.at<unsigned char>(data.PX_I(k), data.PX_J(k))) << std::endl
         << "sw= " << data.sw.at<double>(data.PX_I(k), data.PX_J(k)) << std::endl
          // << "swb= " << data.sw_blur.at<float>(data.PX_I(k), data.PX_J(k)) << std::endl
          ;

      if (data.FUV.rows() > 0) {
        double u = 0;
        double v = 0;
        for (int l = 0; l < 3; ++l) {
          u += data.PX_bc(k, l) * data.UV(data.FUV(f, l), 0);
          v += data.PX_bc(k, l) * data.UV(data.FUV(f, l), 1);
        }
        ss << "u=" << u << std::endl << "v=" << v << std::endl;
      }

      L__pixels.str[k] = ss.str();
      L__pixels.pos.row(k) << data.PX_XY(k).real(), data.PX_XY(k).imag(), 0.0;
    }
  }

  // ===========================================================================

  L__vertex.clr.resize(1, 3);
  L__vertex.clr << CustomColors::Black; // all black
  L__vertex.pos = data.V;
  L__vertex.str.resize(data.V.rows());
  for (int v = 0; v < data.V.rows(); ++v) {
    std::stringstream ss;
    ss << "[" << std::fixed << std::setprecision(2) << data.V(v, 0) //
       << "," << std::fixed << std::setprecision(2) << data.V(v, 1) << "]";
    L__vertex.str[v] = ss.str();
  }

  // ===========================================================================

  {
    L__curves.pos.resize(data.CH.size(), 3);
    L__curves.clr.resize(data.CH.size(), 3);
    L__curves.str.resize(data.CH.size());
    int cui = 0;
    int qvi;
    std::stringstream ss;
    for (const auto& _chain : data.CH) {
      qvi = _chain.qverts[_chain.qverts.size() / 2]; // middle q-vertex

      ss.str(std::string()); // clear
      ss << "cu#" << cui
         << std::endl
         // << " cost = " << _chain.fitting_cost << std::endl
         // << "_cost = "
         << (_chain.fitting_cost / (double)_chain.samples.size());
      L__curves.str[cui] = ss.str();
      L__curves.pos.row(cui) << data.QV[qvi].x, data.QV[qvi].y, 0.0;
      L__curves.clr.row(cui) << CustomColors::Palette(cui);
      cui++;
    }
  }

  // ===========================================================================

  if (data.V_isNonManifold.rows() > 0) {
    P__isNonManifold.clr = CustomColors::Red;
    P__isNonManifold.pos.resize(data.V_isNonManifold.count(), 3);
    int row = 0;
    for (int v = 0; v < data.V_isNonManifold.rows(); ++v)
      if (data.V_isNonManifold(v))
        P__isNonManifold.pos.row(row++) << data.V.row(v);
  }

  // ===========================================================================
}
