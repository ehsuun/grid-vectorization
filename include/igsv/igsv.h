//
// Created by Tibor Stanko on 16/06/2020.
//

#ifndef IGSV_BASE_IGSV_H
#define IGSV_BASE_IGSV_H

#include "common/CornerType.h"
#include "common/CustomColors.h"
#include "common/defs.h"
#include "common/FaceData.h"
#include "common/Options.h"
#include "common/print_timings.h"
#include "common/VectorizationData.h"

#include "parametrization/FrameFieldGenerator.h"
#include "parametrization/IGMSolver.h"
#include "parametrization/compute_constraints.h"
#include "parametrization/compute_narrow_band.h"
#include "parametrization/compute_scale.h"
#include "parametrization/compute_snapping_weights.h"
#include "parametrization/custom_comb_frame_field.h"
#include "parametrization/find_nearest.h"
#include "parametrization/process_uv.h"
#include "parametrization/read_sketch.h"
#include "parametrization/simplify_cut_graph.h"
#include "parametrization/trace_streamlines.h"
#include "parametrization/triangulate.h"

#include "extraction/assign_pixels_to_qedges.h"
#include "extraction/assign_samples_to_chains.h"
#include "extraction/determine_chain_adjacency.h"
#include "extraction/export_svg.h"
#include "extraction/extract_qedges.h"
#include "extraction/extract_qports.h"
#include "extraction/extract_qvertices.h"
#include "extraction/fit_chains.h"
#include "extraction/init_base_labels.h"
#include "extraction/init_chains.h"
#include "extraction/parametrize_chains.h"
#include "extraction/remove_short_chains.h"
#include "extraction/store_orientation.h"

#endif //IGSV_BASE_IGSV_H
