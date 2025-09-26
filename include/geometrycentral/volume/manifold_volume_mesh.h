#pragma once

#include "geometrycentral/combinatorial-maps/combinatorial_map.h"

namespace geometrycentral {
namespace volume {
using ManifoldVolumeMesh = combinatorial_map::CombinatorialMap<3>;

using Vertex = combinatorial_map::Cell<0, 3>;
using Edge = combinatorial_map::Cell<1, 3>;
using Face = combinatorial_map::Cell<2, 3>;
using Cell = combinatorial_map::Cell<3, 3>;

using Dart = combinatorial_map::Dart<3>;
using VertexCorner = combinatorial_map::Incidence<0, 3, 3>;
using EdgeCorner = combinatorial_map::Incidence<1, 3, 3>;
using FaceCorner = combinatorial_map::Incidence<0, 2, 3>;

template <typename T>
using VertexData = MeshData<Vertex, T>;

template <typename T>
using EdgeData = MeshData<Edge, T>;

template <typename T>
using FaceData = MeshData<Face, T>;

template <typename T>
using CellData = MeshData<Cell, T>;

template <typename T>
using VertexCornerData = MeshData<VertexCorner, T>;

template <typename T>
using EdgeCornerData = MeshData<EdgeCorner, T>;

template <typename T>
using FaceCornerData = MeshData<FaceCorner, T>;

} // namespace volume
} // namespace geometrycentral
