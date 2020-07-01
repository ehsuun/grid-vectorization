//
// Created by Tibor Stanko on 08/06/2020.
//

// #include "igsv/Gui.h"
#include "igsv/State.h"

int main(int argc, char* argv[]) {

  if (argc != 3) {
    std::cout << "format: ./IGSSV_bin input_file output_file" << std::endl;
    return EXIT_FAILURE;
  }

  igsv::State _s;

  if (!_s.read_sketch(argv[1])) {
    std::cout << "error: read_sketch() failed, filename=" << argv[1] << std::endl;
    return EXIT_FAILURE;
  }

  if (!_s.triangulate_sketch()) {
    std::cout << "error: triangulate_sketch() failed" << std::endl;
    return EXIT_FAILURE;
  }

  if (!_s.compute_frame_field()) {
    std::cout << "error: compute_frame_field() failed" << std::endl;
    return EXIT_FAILURE;
  }

  if (!_s.compute_parametrization()) {
    std::cout << "error: compute_parametrization() failed" << std::endl;
    return EXIT_FAILURE;
  }

  if (!_s.extract_curves()) {
    std::cout << "error: extract_curves() failed" << std::endl;
    return EXIT_FAILURE;
  }

  if (!_s.export_svg(argv[2])) {
    std::cout << "error: export_svg() failed, filename=" << argv[2] << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

  // // ---- Render ---- dark pixels
  // std::cout << "gui: dark pixels" << std::endl;
  // igsv::drawable::points dark_pixels(_s.band._n_dark);
  // for (int i = 0; i < _s.band._n_dark; ++i) {
  //   dark_pixels.xy[i][0] = _s.band._PX_XY(i).real();
  //   dark_pixels.xy[i][1] = _s.band._PX_XY(i).imag();
  //   if (_s.band._PX_fid(i) == -1) {
  //     dark_pixels.rgb[i][0] = 1.0;
  //     dark_pixels.rgb[i][1] = 0.0;
  //     dark_pixels.rgb[i][2] = 0.0;
  //   } else {
  //     dark_pixels.rgb[i][0] = 0.0;
  //     dark_pixels.rgb[i][1] = 1.0;
  //     dark_pixels.rgb[i][2] = 0.0;
  //   }
  // }
  //
  // // ---- Render ---- nearest verts
  // std::cout << "gui: nearest verts" << std::endl;
  // igsv::drawable::edges nearest_verts(_s.band._n_dark * 2, _s.band._n_dark);
  // for (int i = 0; i < _s.band._n_dark; ++i) {
  //   nearest_verts.xy[2 * i + 0][0] = _s.band._PX_XY(i).real();
  //   nearest_verts.xy[2 * i + 0][1] = _s.band._PX_XY(i).imag();
  //   nearest_verts.xy[2 * i + 1][0] = _s.band._PX_nearestXY(i).real();
  //   nearest_verts.xy[2 * i + 1][1] = _s.band._PX_nearestXY(i).imag();
  //   nearest_verts.idx[i][0]        = 2 * i + 0;
  //   nearest_verts.idx[i][1]        = 2 * i + 1;
  //   if (_s.band._PX_fid(i) == -1) {
  //     dark_pixels.rgb[i][0] = 1.0;
  //     dark_pixels.rgb[i][1] = 0.0;
  //     dark_pixels.rgb[i][2] = 0.0;
  //   } else {
  //     dark_pixels.rgb[i][0] = 0.0;
  //     dark_pixels.rgb[i][1] = 1.0;
  //     dark_pixels.rgb[i][2] = 0.0;
  //   }
  // }
  //
  // // ----- Render : Sobel tangents -----
  // std::cout << "gui: sobel tangents" << std::endl;
  // igsv::drawable::edges sobel_tangents(_s.band._n_dark * 2, _s.band._n_dark);
  // for (int k = 0; k < _s.band._n_dark; ++k) {
  //   sobel_tangents.xy[2 * k][0]      = _s.band._PX_XY(k).real();
  //   sobel_tangents.xy[2 * k][1]      = _s.band._PX_XY(k).imag();
  //   sobel_tangents.xy[2 * k + 1][0]  = _s.band._PX_XY(k).real() + _s.band._PX_SobelTangent(k).real();
  //   sobel_tangents.xy[2 * k + 1][1]  = _s.band._PX_XY(k).imag() + _s.band._PX_SobelTangent(k).imag();
  //   sobel_tangents.rgb[2 * k][0]     = 1.0f;
  //   sobel_tangents.rgb[2 * k][1]     = 0.0f;
  //   sobel_tangents.rgb[2 * k][2]     = 0.0f;
  //   sobel_tangents.rgb[2 * k + 1][0] = 1.0f;
  //   sobel_tangents.rgb[2 * k + 1][1] = 0.0f;
  //   sobel_tangents.rgb[2 * k + 1][2] = 0.0f;
  //   sobel_tangents.idx[k][0]         = 2 * k;
  //   sobel_tangents.idx[k][1]         = 2 * k + 1;
  // }
  //
  // // ----- Render : Narrow band mesh -----
  // std::cout << "gui: nb mesh" << std::endl;
  // igsv::drawable::triangles mesh_nb(_s.mesh._V_nb.rows(), _s.mesh._F_nb.rows());
  // igsv::drawable::edges mesh_nb_edges(_s.mesh._V_nb.rows(), _s.mesh._F_nb.rows() * 3);
  // for (int v = 0; v < _s.mesh._V_nb.rows(); ++v) {
  //   for (int k = 0; k < 2; ++k) {
  //     mesh_nb.xy[v][k]       = _s.mesh._V_nb(v, k);
  //     mesh_nb_edges.xy[v][k] = _s.mesh._V_nb(v, k);
  //   }
  //   mesh_nb.rgb[v]       = { 0.9f, 0.6f, 0.3f };
  //   mesh_nb_edges.rgb[v] = { 0.8f, 0.4f, 0.2f };
  // }
  //
  // for (int f = 0; f < _s.mesh._F_nb.rows(); ++f)
  //   for (int k = 0; k < 3; ++k) {
  //     mesh_nb.idx[f][k] = _s.mesh._F_nb(f, k);
  //
  //     mesh_nb_edges.idx[3 * f + k][0] = _s.mesh._F_nb(f, k);
  //     mesh_nb_edges.idx[3 * f + k][1] = _s.mesh._F_nb(f, (k + 1) % 3);
  //   }

  // // ----- Render : Narrow band mesh with UVs -----
  // std::cout << "gui: nb mesh with uvs" << std::endl;
  // igsv::drawable::triangles mesh_nb_uvs(_s.mesh._F_nb.rows() * 3, _s.mesh._F_nb.rows());
  // for (int f = 0; f < _s.mesh._F_nb.rows(); ++f) {
  //   // int v0 = _s.mesh._F_nb(f, 0);
  //   // int v1 = _s.mesh._F_nb(f, 1);
  //   // int v2 = _s.mesh._F_nb(f, 2);
  //   //
  //   // int u0 = _s.mesh._FUV(f, 0);
  //   // int u1 = _s.mesh._FUV(f, 1);
  //   // int u2 = _s.mesh._FUV(f, 2);
  //
  //   for (int k = 0; k < 3; ++k) {
  //     int v = _s.mesh._F_nb(f, k);
  //     int u = _s.param._FUV(f, k);
  //
  //     mesh_nb_uvs.xy[3 * f + k][0] = _s.mesh._V_nb(v, 0);
  //     mesh_nb_uvs.xy[3 * f + k][1] = _s.mesh._V_nb(v, 1);
  //
  //     mesh_nb_uvs.uv[3 * f + k][0] = _s.param._UV(u, 0);
  //     mesh_nb_uvs.uv[3 * f + k][1] = _s.param._UV(u, 1);
  //
  //     mesh_nb_uvs.idx[f][k] = 3 * f + k;
  //   }
  // }
  // for (int v = 0; v < _s.mesh._V_nb.rows(); ++v) {
  //   for (int k = 0; k < 2; ++k) {
  //     mesh_nb.xy[v][k]       = _s.mesh._V_nb(v, k);
  //     mesh_nb.uv[v][k]       = _s.param._UV(v, k);
  //   }
  // }

  // for (int f = 0; f < _s.mesh._F_nb.rows(); ++f)
  //   for (int k = 0; k < 3; ++k) {
  //     mesh_nb.idx[f][k] = _s.mesh._F_nb(f, k);
  //   }

  // // ----- Render : Frame field vectors -----
  // std::cout << "gui: frame field vectors" << std::endl;
  // igsv::drawable::edges frame_field_vectors(_s.mesh._V_nb.rows() * 4, _s.mesh._V_nb.rows() * 2);
  // for (int v = 0; v < _s.mesh._V_nb.rows(); ++v) {
  //
  //   frame_field_vectors.xy[4 * v + 0][0] = _s.mesh._V_nb(v, 0);
  //   frame_field_vectors.xy[4 * v + 0][1] = _s.mesh._V_nb(v, 1);
  //   frame_field_vectors.xy[4 * v + 2][0] = _s.mesh._V_nb(v, 0);
  //   frame_field_vectors.xy[4 * v + 2][1] = _s.mesh._V_nb(v, 1);
  //   frame_field_vectors.xy[4 * v + 1][0] = _s.mesh._V_nb(v, 0) + std::cos(_s.field._X0_arg_vtx(v));
  //   frame_field_vectors.xy[4 * v + 1][1] = _s.mesh._V_nb(v, 1) + std::sin(_s.field._X0_arg_vtx(v));
  //   frame_field_vectors.xy[4 * v + 3][0] = _s.mesh._V_nb(v, 0) + std::cos(_s.field._X1_arg_vtx(v));
  //   frame_field_vectors.xy[4 * v + 3][1] = _s.mesh._V_nb(v, 1) + std::sin(_s.field._X1_arg_vtx(v));
  //
  //   frame_field_vectors.rgb[4 * v + 0][0] = frame_field_vectors.rgb[4 * v + 1][0] = 1.0;
  //   frame_field_vectors.rgb[4 * v + 0][1] = frame_field_vectors.rgb[4 * v + 1][1] = 0.0;
  //   frame_field_vectors.rgb[4 * v + 0][2] = frame_field_vectors.rgb[4 * v + 1][2] = 0.0;
  //
  //   frame_field_vectors.rgb[4 * v + 2][0] = frame_field_vectors.rgb[4 * v + 3][0] = 0.0;
  //   frame_field_vectors.rgb[4 * v + 2][1] = frame_field_vectors.rgb[4 * v + 3][1] = 0.0;
  //   frame_field_vectors.rgb[4 * v + 2][2] = frame_field_vectors.rgb[4 * v + 3][2] = 1.0;
  //
  //   frame_field_vectors.idx[2 * v + 0][0] = 4 * v + 0;
  //   frame_field_vectors.idx[2 * v + 0][1] = 4 * v + 1;
  //   frame_field_vectors.idx[2 * v + 1][0] = 4 * v + 2;
  //   frame_field_vectors.idx[2 * v + 1][1] = 4 * v + 3;
  // }
  //
  // // ----- Render : Frame field vectors -----
  // std::cout << "gui: frame field vectors (face)" << std::endl;
  // igsv::drawable::edges frame_field_vectors_face(_s.mesh._BC_nb.rows() * 4, _s.mesh._BC_nb.rows() * 2);
  // for (int f = 0; f < _s.mesh._BC_nb.rows(); ++f) {
  //
  //   frame_field_vectors_face.xy[4 * f + 0][0] = _s.mesh._BC_nb(f, 0);
  //   frame_field_vectors_face.xy[4 * f + 0][1] = _s.mesh._BC_nb(f, 1);
  //   frame_field_vectors_face.xy[4 * f + 2][0] = _s.mesh._BC_nb(f, 0);
  //   frame_field_vectors_face.xy[4 * f + 2][1] = _s.mesh._BC_nb(f, 1);
  //   frame_field_vectors_face.xy[4 * f + 1][0] = _s.mesh._BC_nb(f, 0) + _s.field._X0_nonunit_comb(f, 0);
  //   frame_field_vectors_face.xy[4 * f + 1][1] = _s.mesh._BC_nb(f, 1) + _s.field._X0_nonunit_comb(f, 1);
  //   frame_field_vectors_face.xy[4 * f + 3][0] = _s.mesh._BC_nb(f, 0) + _s.field._X1_nonunit_comb(f, 0);
  //   frame_field_vectors_face.xy[4 * f + 3][1] = _s.mesh._BC_nb(f, 1) + _s.field._X1_nonunit_comb(f, 1);
  //
  //   frame_field_vectors_face.rgb[4 * f + 0][0] = frame_field_vectors_face.rgb[4 * f + 1][0] = 1.0;
  //   frame_field_vectors_face.rgb[4 * f + 0][1] = frame_field_vectors_face.rgb[4 * f + 1][1] = 0.0;
  //   frame_field_vectors_face.rgb[4 * f + 0][2] = frame_field_vectors_face.rgb[4 * f + 1][2] = 0.0;
  //
  //   frame_field_vectors_face.rgb[4 * f + 2][0] = frame_field_vectors_face.rgb[4 * f + 3][0] = 0.0;
  //   frame_field_vectors_face.rgb[4 * f + 2][1] = frame_field_vectors_face.rgb[4 * f + 3][1] = 0.0;
  //   frame_field_vectors_face.rgb[4 * f + 2][2] = frame_field_vectors_face.rgb[4 * f + 3][2] = 1.0;
  //
  //   frame_field_vectors_face.idx[2 * f + 0][0] = 4 * f + 0;
  //   frame_field_vectors_face.idx[2 * f + 0][1] = 4 * f + 1;
  //   frame_field_vectors_face.idx[2 * f + 1][0] = 4 * f + 2;
  //   frame_field_vectors_face.idx[2 * f + 1][1] = 4 * f + 3;
  // }
  //
  // // ----- Render : Cuts -----
  // igsv::drawable::edges custom_cuts;
  // for (int f = 0; f < _s.mesh._F_nb.rows(); ++f)
  //   for (int k = 0; k < 3; ++k)
  //     if (_s.field._E_customCut(f, k)) {
  //       int v0 = _s.mesh._F_nb(f, k);
  //       int v1 = _s.mesh._F_nb(f, (k + 1) % 3);
  //       custom_cuts.xy.push_back({ static_cast<float>(_s.mesh._V_nb(v0, 0)), //
  //                                  static_cast<float>(_s.mesh._V_nb(v0, 1)) });
  //       custom_cuts.xy.push_back({ static_cast<float>(_s.mesh._V_nb(v1, 0)), //
  //                                  static_cast<float>(_s.mesh._V_nb(v1, 1)) });
  //       if (_s.mesh._TT(f, k) == -1) {
  //         custom_cuts.rgb.push_back({ 0., 1., 0. });
  //         custom_cuts.rgb.push_back({ 0., 1., 0. });
  //       } else {
  //         custom_cuts.rgb.push_back({ 1., 0., 1. });
  //         custom_cuts.rgb.push_back({ 1., 0., 1. });
  //       }
  //       custom_cuts.idx.push_back({ static_cast<unsigned int>(custom_cuts.xy.size() - 2), //
  //                                   static_cast<unsigned int>(custom_cuts.xy.size() - 1) });
  //     }
  //
  // // ---- Render ---- tangent direction labels
  // std::cout << "gui: tandir labels" << std::endl;
  // igsv::drawable::points tandir_labels(_s.mesh._F_nb.rows());
  // float r, g, b;
  // for (int i = 0; i < _s.mesh._F_nb.rows(); ++i) {
  //   tandir_labels.xy[i][0] = _s.mesh._BC_nb(i, 0);
  //   tandir_labels.xy[i][1] = _s.mesh._BC_nb(i, 1);
  //   if (_s.field._F_labels(i) == igsv::CORNER_TYPE_INTEGER_U) {
  //     tandir_labels.rgb[i][0] = 1.0;
  //     tandir_labels.rgb[i][1] = 0.0;
  //     tandir_labels.rgb[i][2] = 0.0;
  //   } else if (_s.field._F_labels(i) == igsv::CORNER_TYPE_INTEGER_V) {
  //     tandir_labels.rgb[i][0] = 0.0;
  //     tandir_labels.rgb[i][1] = 0.0;
  //     tandir_labels.rgb[i][2] = 1.0;
  //   } else {
  //     tandir_labels.rgb[i][0] = 0.3;
  //     tandir_labels.rgb[i][1] = 0.3;
  //     tandir_labels.rgb[i][2] = 0.3;
  //   }
  // }

  // // ----- Render : base mesh (two triangles) -----
  // igsv::drawable::triangles _base(4, 2);
  // _base.xy[0]  = { 0.0f, 0.0f };
  // _base.xy[1]  = { (float)_s.sketch._input.cols, 0.0f };
  // _base.xy[2]  = { (float)_s.sketch._input.cols, (float)_s.sketch._input.rows };
  // _base.xy[3]  = { 0.0f, (float)_s.sketch._input.rows };
  // _base.uv[0]  = { 0.0f, 0.0f };
  // _base.uv[1]  = { 1.0f, 0.0f };
  // _base.uv[2]  = { 1.0f, 1.0f };
  // _base.uv[3]  = { 0.0f, 1.0f };
  // _base.rgb[0] = { 1.0f, 1.0f, 1.0f };
  // _base.rgb[1] = { 1.0f, 1.0f, 1.0f };
  // _base.rgb[2] = { 1.0f, 1.0f, 1.0f };
  // _base.rgb[3] = { 1.0f, 1.0f, 1.0f };
  // _base.idx[0] = { 0, 1, 2 };
  // _base.idx[1] = { 2, 3, 0 };
  //
  // // ----- Render : LineMap -----
  // std::cout << "gui: line map" << std::endl;
  // igsv::drawable::edges line_map;
  // for (const auto& s : _s.extr._lineMap.segmentsU) {
  //   line_map.xy.push_back({ (float)s.xyz[0].x(), (float)s.xyz[0].y() });
  //   line_map.xy.push_back({ (float)s.xyz[1].x(), (float)s.xyz[1].y() });
  //   line_map.rgb.push_back({ 1.0f, 0.0f, 0.0f });
  //   line_map.rgb.push_back({ 1.0f, 0.0f, 0.0f });
  //   line_map.idx.push_back({ (unsigned)line_map.xy.size() - 2, (unsigned)line_map.xy.size() - 1 });
  // }
  // for (const auto& s : _s.extr._lineMap.segmentsV) {
  //   line_map.xy.push_back({ (float)s.xyz[0].x(), (float)s.xyz[0].y() });
  //   line_map.xy.push_back({ (float)s.xyz[1].x(), (float)s.xyz[1].y() });
  //   line_map.rgb.push_back({ 0.0f, 0.0f, 1.0f });
  //   line_map.rgb.push_back({ 0.0f, 0.0f, 1.0f });
  //   line_map.idx.push_back({ (unsigned)line_map.xy.size() - 2, (unsigned)line_map.xy.size() - 1 });
  // }
  //
  // // ---- Render ---- q-vertices
  // std::cout << "gui: q-vertices" << std::endl;
  // igsv::drawable::points qvertices(_s.extr._QV.size());
  // {
  //   int qvid = 0;
  //   for (const auto& qv : _s.extr._QV) {
  //     qvertices.xy[qvid][0]  = (float)qv.x;
  //     qvertices.xy[qvid][1]  = (float)qv.y;
  //     qvertices.rgb[qvid][0] = 0.2f;
  //     qvertices.rgb[qvid][1] = 0.8f;
  //     qvertices.rgb[qvid][2] = 0.4f;
  //     qvid++;
  //   }
  // }
  //
  // // ----- Render : qedges -----
  // std::cout << "gui: qedges" << std::endl;
  // igsv::drawable::edges qedges(_s.extr._QV.size(), _s.extr._QE.size());
  // {
  //   int qvid = 0;
  //   for (const auto& qv : _s.extr._QV) {
  //     qedges.xy[qvid][0]  = (float)qv.x;
  //     qedges.xy[qvid][1]  = (float)qv.y;
  //     qedges.rgb[qvid][0] = 0.2f;
  //     qedges.rgb[qvid][1] = 0.6f;
  //     qedges.rgb[qvid][2] = 0.3f;
  //     qvid++;
  //   }
  // }
  // {
  //   int qeid = 0;
  //   for (const auto& qe : _s.extr._QE) {
  //     qedges.idx[qeid][0] = qe.qvi0;
  //     qedges.idx[qeid][1] = qe.qvi1;
  //     qeid++;
  //   }
  // }
  //
  // // ----- Render : Bezier curves -----
  // std::cout << "gui: Bezier curves" << std::endl;
  // igsv::drawable::edges bezier_curves;
  // {
  //   const int n_points_per_segment = 100;
  //
  //   double t0, t1, t, dt, c0, c1, c2, c3, d0, d1, d2, d3, x, y;
  //
  //   bool start_new_chain = true;
  //
  //   int row = 0;
  //   int cui = 0;
  //
  //   for (const auto& _chain : _s.extr._CH) {
  //
  //     start_new_chain = true;
  //
  //     for (int i = 0; i < _chain.CubicSpline.rows(); ++i) {
  //
  //       t0 = _chain.CubicSpline(i, 0);
  //       t1 = _chain.CubicSpline(i, 1);
  //
  //       c0 = _chain.CubicSpline(i, 2);
  //       c1 = _chain.CubicSpline(i, 3);
  //       c2 = _chain.CubicSpline(i, 4);
  //       c3 = _chain.CubicSpline(i, 5);
  //
  //       d0 = _chain.CubicSpline(i, 6);
  //       d1 = _chain.CubicSpline(i, 7);
  //       d2 = _chain.CubicSpline(i, 8);
  //       d3 = _chain.CubicSpline(i, 9);
  //
  //       dt = (t1 - t0) / ((double)n_points_per_segment - 1.);
  //
  //       // std::cout << "segment : t0=" << t0 << ", t1=" << t1 << std::endl;
  //
  //       t = 0;
  //       // std::cout << "  t=" << t << std::endl;
  //       x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
  //       y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;
  //
  //       for (int k = 0; k < n_points_per_segment - 1; ++k) {
  //
  //         if (start_new_chain) {
  //           bezier_curves.xy.push_back({ (float)x, (float)y });
  //           bezier_curves.rgb.push_back({ 1.0f, 0.0f, 0.0f });
  //           start_new_chain = false;
  //         }
  //
  //         t += dt;
  //         // std::cout << "  t=" << t << std::endl;
  //         x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
  //         y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;
  //
  //         bezier_curves.xy.push_back({ (float)x, (float)y });
  //         if (k == n_points_per_segment - 2)
  //           bezier_curves.rgb.push_back({ 0.0f, 1.0f, 0.0f });
  //         else
  //           bezier_curves.rgb.push_back({ 0.0f, 0.0f, 0.0f });
  //         bezier_curves.idx.push_back({ (unsigned)bezier_curves.xy.size() - 2, //
  //                                       (unsigned)bezier_curves.xy.size() - 1 });
  //       }
  //     }
  //   }
  // }
  //
  // //###########################################################################
  // //###########################################################################
  // //###########################################################################
  //
  // std::cout << "gui: init" << std::endl;
  //
  // cv::Mat grid_texture(100, 100, CV_8UC3);
  // grid_texture.setTo(cv::Vec3b(255, 255, 255));
  //
  // cv::copyMakeBorder(grid_texture, grid_texture, 1, 1, 0, 0, CV_HAL_BORDER_CONSTANT, cv::Vec3b(255, 0, 0));
  // cv::copyMakeBorder(grid_texture, grid_texture, 0, 0, 1, 1, CV_HAL_BORDER_CONSTANT, cv::Vec3b(0, 0, 255));
  //
  // igsv::Gui gui(_s.sketch._input, _s.sketch._mask, grid_texture);
  // gui.add_texture("grid", &grid_texture);
  // gui.add_texture("input", &_s.sketch._input);
  // // gui.add_texture("bw", &_s.sketch._bw);
  // // gui.add_texture("dt_01", &_s.sketch._dt, true);
  // // gui.add_texture("idt_01", &_s.sketch._idt_01, true);
  // // gui.add_texture("sw_01", &_s.sketch._sw_01, true);
  //
  // // gui.add_points("dark pixels", dark_pixels);
  // // gui.add_points("tangent direction labels", tandir_labels);
  // gui.add_points("q-vertices", qvertices);
  //
  // // gui.add_edges("sobel tangents", sobel_tangents);
  // // gui.add_edges("frame field vectors (vertex)", frame_field_vectors);
  // // gui.add_edges("combed frame field", frame_field_vectors_face);
  // // gui.add_edges("custom cuts", custom_cuts);
  // // gui.add_edges("nearest verts", nearest_verts);
  // gui.add_edges("line map", line_map);
  // gui.add_edges("q-edges", qedges);
  // gui.add_edges("bezier curves", bezier_curves);
  //
  // // gui.add_edges("mesh nb edges", mesh_nb_edges);
  //
  // // gui.add_triangles("mesh nb", mesh_nb);
  // gui.add_triangles("base triangle mesh", _base);
  // gui.add_triangles("mesh nb uv", mesh_nb_uvs);
  //
  // std::cout << "gui: render" << std::endl;
  // gui.render();

  return EXIT_SUCCESS;
}
