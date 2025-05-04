#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/poisson_disk_sampler.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

struct TraceStreamlineOptions {
  size_t maxSegments = 10;                                 // maximum number of segments in streamline
  double maxLen = std::numeric_limits<double>::infinity(); // maximum length of streamline curve on mesh
  FaceData<int>* nVisits = nullptr;         // if set, increments nVisits[f] each time streamline crosses face f
  const FaceData<int>* maxVisits = nullptr; // if maxVisits and nVisits are set, streamline terminates if it hits a face
                                            // f with nVisits[f] >= maxVisits[f]
  int verbosity = 0;                        // if set > 0, print out debug information
};
// Trace a streamline of field starting from surface point pStart.
// Warning: If nSym > 1, field should contain a representative vector (i.e. one of the nSym many vectors in the face),
// _not_ the "power vector" used e.g. by the methods in direction_fields.h
std::vector<SurfacePoint> traceStreamline(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, SurfacePoint pStart,
                                          const FaceData<Vector2>& field, size_t nSym = 1,
                                          TraceStreamlineOptions opt = TraceStreamlineOptions());

// Convenience function to trace streamlines starting from a collection of seed points sampled on the mesh via Poisson
// disk sampling
std::vector<std::vector<SurfacePoint>> traceManyStreamlines(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                                                            const FaceData<Vector2>& field, size_t nSym = 1,
                                                            TraceStreamlineOptions opt = TraceStreamlineOptions(),
                                                            PoissonDiskOptions pdOpt = PoissonDiskOptions());

//=== Streamline output helpers
struct SvgCurveOptions {
  double worldspaceStrokeWidth = .05;   // set stroke widths so that curves appear approximately this large on the mesh
  std::string backgroundColor = "#fff"; // svg background color
  std::string curveColor = "#000000";   // svg curve color
  std::function<std::string(size_t, size_t)> curveColorFunction; // if set, i'th curve gets color curveColorFunction(i,
                                                                 // curves.size()). Overrides curveColor
  double imageSize = 500;                                        // size of output svg viewBox
};
// Write streamlines to an SVG texture located at `filepath` using texture coordinates `param` on the input mesh
void draw_mesh_curves_to_svg(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom, const CornerData<Vector2>& param,
                             const std::vector<std::vector<SurfacePoint>>& curves, std::string filepath,
                             SvgCurveOptions opt = SvgCurveOptions());
} // namespace surface
} // namespace geometrycentral
