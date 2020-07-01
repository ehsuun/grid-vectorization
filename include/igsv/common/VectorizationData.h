//=============================================================================
//
//  CLASS : VectorizationData
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <map>
#include <string>

#include <Eigen/Core>
#include <igl/serialize.h>
#include <opencv2/opencv.hpp>

#include "../common/CornerType.h"
#include "../common/FaceData.h"
#include "../common/Options.h"
#include "../common/defs.h"

#include "../extraction_entities/Chain.h"
#include "../extraction_entities/LineSegment.h"
#include "../extraction_entities/Pixel.h"
#include "../extraction_entities/QEdge.h"
#include "../extraction_entities/QPort.h"
#include "../extraction_entities/QVertex.h"

#include "../parametrization/_cgal_delaunay_2.h" // using IndexedDT

//== NAMESPACES ===============================================================

namespace IGSV {

  //== TYPEDEFS ===============================================================

  using VectorXb  = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
  using MatrixX3b = Eigen::Matrix<bool, Eigen::Dynamic, 3>;
  using MatrixXb  = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

  //== CLASS DEFINITION =======================================================

  class VectorizationData : public igl::Serializable {
  public:
    VectorizationData() { clear(); }

    void clear();
    void clear_extraction();

    void InitSerialization() override;
    bool PreSerialization() const override;
    void PostSerialization() const override;

    // Wrappers
    bool compute_no_gui();

    bool init_sketch();
    bool init_mesh(bool _compute_from_scratch = true);
    bool init_frame_field(bool _compute_from_scratch = true, bool _trace_streamlines = true);
    bool init_parametrization(bool _compute_from_scratch = true);
    bool init_extraction();

    // Methods : sketch
    bool read_sketch();
    bool process_sketch();
    bool compute_constraints();
    bool estimate_max_area();

    // Methods : mesh
    bool triangulate();
    bool process_mesh();
    bool find_nearest();
    bool compute_narrow_band();
    bool compute_snapping_weights(bool print_message = true);

    // Methods : frame field
    bool compute_frame_field_narrow_band();
    bool comb_frame_field();
    bool trace_streamlines();
    bool simplify_cut_graph_iterate();

    // Methods : parametrization
    bool compute_scale();
    bool parametrize(double& time);

    // Methods : parametrization post-process
    bool store_orientation();
    bool sanitize_parametrization();
    bool precompute_transition_fn();

    // Methods : extraction
    bool extract_line_map();

    // Methods : QEx-like extraction
    bool extract_qvertices();
    bool extract_qports();
    bool extract_qedges();

    // Methods : extraction, get samples and base labels
    bool assign_pixels_to_qedges();
    bool init_base_labels(bool split_at_nodes);
    bool init_chains();
    bool parametrize_chains();
    bool assign_samples_to_chains();
    bool remove_short_chains();
    bool determine_chain_adjacency();
    bool fit_chains();

    // Methods : export, serialization, timings, stats
    // bool export_info() const;
    bool export_svg() const;
    // bool check_if_state_exists() const;

    template <typename StreamT>
    void print_timings(StreamT&) const;

    // std::string get_filename(FILE_FLAG, bool create_dir = false) const;

    // ===========================================================================
    // TIMINGS
    // ===========================================================================
    std::map<std::string, double> timings;

    // ===========================================================================
    // OPTIONS
    // ===========================================================================
    Options opts;

    // ===========================================================================
    // SKETCH
    // ===========================================================================
    cv::Mat input;              // input image
    cv::Mat detail_mask;        // detail mask
    bool mask_read_from_file{}; // has detail mask been read?
    cv::Mat grey;               // greyscale
    cv::Mat grey_blurred;       // blurred greyscale, used to compute Sobel tangents
    cv::Mat bw;                 // binarized
    cv::Mat dt;                 // "skeleton" of white regions; dist: white pixel -> nearest dark
    cv::Mat idt_01;             // -||- normalized to 01
    cv::Mat sw;                 // per-pixel stroke width
    cv::Mat sw_01;              //  -||- normalized to 01
    double sw_max{};            // stroke width : max
    double sw_avg{};            // stroke width : avg
    double sw_stddev{};         // stroke width : stddev

    // ===========================================================================
    // NARROW BAND (dark pixels)
    // ===========================================================================
    unsigned n_dark;                  // number of dark pixels
    Eigen::MatrixXd pix_weight;       // weights (all pixels)
    Eigen::MatrixXi pix2band;         // conversion: pixel indices (i,j) -> band indices (k), pix2band(i,j) = k
    Eigen::VectorXi PX_I;             // dark pixel indices (ij)
    Eigen::VectorXi PX_J;             // dark pixel indices (ij)
    Eigen::VectorXi PX_fid;           // dark pixel triangle indices
    Eigen::MatrixXd PX_bc;            // dark pixel barycentric coords w.r.t. its face
    Eigen::VectorXcd PX_XY;           // dark pixel position (complex)
    Eigen::VectorXcd PX_SobelTangent; // Sobel tangent (complex)
    Eigen::VectorXd PX_SobelWeight;   // Sobel tangent confidence (scalar)

    // ===========================================================================
    // MESH, FRAME FIELD, PARAMETRIZATION
    // ===========================================================================
    IndexedDT tri;                                       // CGAL triangulation
    Eigen::MatrixXd V;                                   // vertices
    Eigen::MatrixXi F;                                   // triangles
    Eigen::MatrixXd BC;                                  // triangle barycenters
    Eigen::MatrixXd UV_generated;                        // generated UVs (for image texturing)
    Eigen::MatrixXd UV;                                  // integer-grid map, parametric coordinates
    Eigen::MatrixXi FUV;                                 // integer-grid map, parametric indices
    Eigen::MatrixXi TT;                                  // triangle-triangle adjacency
    Eigen::MatrixXi TTi;                                 // triangle-triangle adjacency
    std::vector<std::vector<int>> VT;                    // vertex-triangle adjacency
    std::vector<std::vector<int>> VTi;                   // vertex-triangle adjacency
    Eigen::VectorXcd PX_nearestXY;                       // dark pixels: nearest vertex in the mesh (xy, complex)
    Eigen::VectorXi PX_nearestID;                        // dark pixels: nearest vertex in the mesh (index)
    Eigen::VectorXd V_snappingWeights;                   // per-vertex weights used for snapping
    Eigen::VectorXi V_singularityIdx;                    // indices of singularities
    VectorXb V_isBlack;                                  // which vertices are dark?
    VectorXb V_isNarrow;                                 // which vertices are in the narrow band?
    VectorXb V_isNarrowBoundary;                         // which vertices are on the bd of the narrow band?
    VectorXb V_isNonManifold;                            // which vertices are non-manifold in the initial narrow band?
    MatrixX3b E_isNarrow;                                // which edges are in the narrow band?
    MatrixX3b E_isNarrowBoundary;                        // which edges are on the bd the narrow band?
    MatrixX3b E_customCut;                               // which edges are cuts?
    Eigen::MatrixXi E_periodJumps;                       // period jumps
    VectorXb F_isBlack;                                  // which faces are dark?
    VectorXb F_isNarrow;                                 // which faces are in the narrow band?
    Eigen::VectorXi F_uvOrientation;                     // orientation of uv triangles
    Eigen::VectorXi F_labels;                            // corner types
    Eigen::VectorXd F_uvScale;                           // per-face scale, isotropic
    Eigen::VectorXcd F_frameField[4];                    // per-face unit vectors as complex numbers
    std::vector<IGSV::FaceData> FaceDataList;            // per-face data
    std::vector<std::pair<int, int>> TT_dualConnections; // dual spanning tree (combing)
    Eigen::MatrixXi FK_snapping_constraints;             //
    Eigen::MatrixXi I_snapping_connections;              //
    Eigen::VectorXd X0_arg_vtx;                          // per-vertex frame field
    Eigen::VectorXd X0_mag_vtx;                          // per-vertex frame field
    Eigen::VectorXd X1_arg_vtx;                          // per-vertex frame field
    Eigen::VectorXd X1_mag_vtx;                          // per-vertex frame field
    Eigen::VectorXd X0_arg;                              // per-face frame field
    Eigen::VectorXd X0_mag;                              // per-face frame field
    Eigen::VectorXd X1_arg;                              // per-face frame field
    Eigen::VectorXd X1_mag;                              // per-face frame field
    Eigen::MatrixXd BIS0_unit;                           // per-face unit bisectors
    Eigen::MatrixXd BIS1_unit;                           // per-face unit bisectors
    Eigen::MatrixXd BIS0_unit_comb;                      // per-face unit bisectors, combed
    Eigen::MatrixXd BIS1_unit_comb;                      // per-face unit bisectors, combed
    Eigen::MatrixXd X0_unit;                             // per-face unit vectors
    Eigen::MatrixXd X1_unit;                             // per-face unit vectors
    Eigen::MatrixXd X0_nonunit;                          // per-face non-unit vectors
    Eigen::MatrixXd X1_nonunit;                          // per-face non-unit vectors
    Eigen::MatrixXd X0_unit_comb;                        // per-face unit vectors, combed
    Eigen::MatrixXd X1_unit_comb;                        // per-face unit vectors, combed
    Eigen::MatrixXd X0_nonunit_comb;                     // per-face non-unit vectors, combed
    Eigen::MatrixXd X1_nonunit_comb;                     // per-face non-unit vectors, combed
    Eigen::MatrixXcd R_transfn;                          // trans fn : rotation
    Eigen::MatrixXcd T_transfn;                          // trans fn : non-integer translation
    Eigen::MatrixXcd T_integer_transfn;                  // trans fn : integer part of the translation (rounded T)
    Eigen::MatrixXcd T_decimal_transfn;                  // trans fn : decimal part of the translation
    MatrixX3b flip_transfn;                              // trans fn : is the direction flipped?

    // frame field stats to show
    int n_ff_iters      = 0;
    int n_ff_degenerate = 0;

    // ===========================================================================
    // EXTRACTION
    // ===========================================================================
    LineMap igm;                                        // traced parametrization
    std::vector<QVertex> QV;                            // list of qvertices
    std::vector<QPort> QP;                              // list of qports
    std::vector<QEdge> QE;                              // list of qedges
    std::vector<Pixel> PX;                              // list of dark pixels
    std::vector<Chain> CH;                              // list of chains
    std::vector<std::vector<int>> F_qports;             // list of qports per face
    std::vector<std::pair<int, int>> QP_QP_connections; // connections between qports
    std::map<std::pair<int, int>, int> QVPair_to_QE;    // qvertex pair --> qedge
    std::vector<std::pair<int, int>> CH_CH_connections; // connections between chains
    std::vector<bool> CH_CH_connections_enabled;        // is a specific chain connection enabled?

    // ===========================================================================
  };
} // namespace IGSV
