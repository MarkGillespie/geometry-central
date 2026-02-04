#pragma once

#include "geometrycentral/volume/manifold_volume_mesh.h"
#include "geometrycentral/volume/vertex_position_geometry.h"

namespace geometrycentral {
namespace volume {
Vector3 weightedEdgeCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Edge e);
Vector3 weightedFaceCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Face f);
Vector3 weightedCellCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Cell c);

EdgeData<Vector3> weightedEdgeCenters(VertexPositionGeometry& geom, const VertexData<double>& weights);
FaceData<Vector3> weightedFaceCenters(VertexPositionGeometry& geom, const VertexData<double>& weights);
CellData<Vector3> weightedCellCenters(VertexPositionGeometry& geom, const VertexData<double>& weights);

// if grad != nullptr, the gradient (with respect to the weights) is added to its contents
double weightedEdgeDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Vertex v, Edge e,
                        VertexData<double>* grad = nullptr);
double weightedFaceDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Edge e, Face f,
                        VertexData<double>* grad = nullptr);
double weightedCellDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Face f, Cell c,
                        VertexData<double>* grad = nullptr);

//==== Hodge-Optimized Triangulation energies. If grad != nullptr, the gradient (with respect to the weights) is added
// to its contents
double star0HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights,
                        VertexData<double>* grad = nullptr);
double star1HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights,
                        VertexData<double>* grad = nullptr);
double star2HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights,
                        VertexData<double>* grad = nullptr);
double star3HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights,
                        VertexData<double>* grad = nullptr);

VertexData<double> optimalStar3Hot22Weights(VertexPositionGeometry& geom);
} // namespace volume
} // namespace geometrycentral
