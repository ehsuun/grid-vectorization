#pragma once
#include "BaseGui.h"
namespace IGSV {
  class VectorizationData;
  class IGSVGui : public BaseGui {
  private:
    GuiEdges E__bezierCurves;          // [  ex ]  Bezier/B-spline curves
    GuiEdges E__bezierPolygons;        // [  ex ]  Bezier/B-spline curves: control polygons
    GuiEdges E__blackPixNearestVertex; // [ tri ]  nearest neighbor in the mesh for each black pixel
    GuiEdges E__customCuts;            // [  ff ]  custom cuts in the combed frame field
    GuiEdges E__lineMapSegments;       // [  ex ]  integer-grid map as line segments
    GuiEdges E__narrowBandBoundary;    // [ tri ]  edges at the narrow band boundary
    GuiEdges E__periodJumps;           // [  ff ]  period jumps
    GuiEdges E__pixSampleConnections;  // [  ex ]  pixel-sample connections
    GuiEdges E__qedges;                // [  ex ]  q-edges
    GuiEdges E__qports;                // [  ex ]  q-ports
    GuiEdges E__snappingConnections;   // [  uv ]  snapping connections between corners
    GuiEdges E__splineCurves;          // [  ex ]  Cubic spline curves
    GuiEdges E__splinePolygons;        // [  ex ]  Cubic spline curves: polygons
    GuiEdges E__streamlines;           // [  ff ]  traced streamlines
                                       //
    GuiFColors F__ffAngles;            // [  ff ]  per-face frame field angles
    GuiFColors F__isCombingOk;         // [  ff ]  combing, problematic faces
    GuiFColors F__narrowBand;          // [ tri ]  narrow band
    GuiFColors F__uvFlipped;           // [  uv ]  flipped triangles
                                       //
    GuiLabels L__corner;               // [ tri ]  labels for corners
    GuiLabels L__curves;               // [  ex ]  labels for curves
    GuiLabels L__faces;                // [ tri ]  labels for faces
    GuiLabels L__ffAngles;             // [  ff ]  labels for frame angles
    GuiLabels L__nTracedBlackPix;      // [  ff ]  labels for black-pixels-along-streamlines counts
    GuiLabels L__pixels;               // [  ex ]  labels for pixels
    GuiLabels L__qexData;              // [  ex ]  labels for qex data
    GuiLabels L__vertex;               // [ tri ]  labels for vertices
                                       //
    GuiPoints P__BRUSH;                // [ DEBUG ]  points traced with a brush
    GuiPoints P__isNonManifold;        // [ DEBUG ]  non-manifold vertices in the initial narrow band mesh
    GuiPoints P__qvertices;            // [  ex ]  q-vertices
    GuiPoints P__singularities;        // [  ff ]  singularities
    GuiPoints P__snappingConstraints;  // [  uv ]  snapping constraints
    GuiPoints P__streamlines;          // [  ff ]  traced streamlines, points
    GuiPoints P__tangentDir;           // [  ff ]  selected tangent direction
                                       //
    GuiVectors V__ffVector0;           // [  ff ]  first frame vector
    GuiVectors V__ffVector1;           // [  ff ]  second frame vector
    GuiVectors V__SobelTangents;       // [  in ]  tangent constraints (Sobel)
    GuiVectors V__SobelNormals;        // [  in ]  tangent constraints (Sobel)

  public:
    explicit IGSVGui(VectorizationData& data);

    void draw_custom_window() override;
    void draw_viewer_menu() override;

    void init(igl::opengl::glfw::Viewer* _viewer) override;
    void store_entities() override;

    bool key_down(int key, int modifiers) override;

    bool mouse_down(int button, int modifier) override;
    bool mouse_up(int button, int modifier) override;
    bool mouse_move(int mouse_x, int mouse_y) override;

  private:
    void update_mesh();
    void update_narrow_band();
    void update_frame_field();
    void update_parametrization();
    void update_extraction();
    void serialize_data();

    bool show_unused_qedges  = false;
    bool show_bisector_field = false;

    bool streamline_random_colors = true;
    float streamline_color1[3];
    float streamline_color2[3];

    enum BrushMode { Inactive, White, Black };
    BrushMode brush_mode;

  public:
    enum QEdgeColorScheme {
      QEDGE_COLOR_UNIFORM = 0,
      QEDGE_COLOR_CHAIN_ID,
      QEDGE_COLOR_CHAIN_T,
      NUMBER_OF_QEDGE_COLORS
    };

  private:
    const char* QEdgeColorLabels[NUMBER_OF_QEDGE_COLORS] = {
      "uniform",  //
      "chain id", //
      "chain t",  //
    };

    int selected_qedge_color = QEDGE_COLOR_UNIFORM;
  };
}; // namespace IGSV
