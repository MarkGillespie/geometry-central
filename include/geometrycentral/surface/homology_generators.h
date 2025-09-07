#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

namespace geometrycentral {
namespace surface {

enum class HomologyType { Absolute, Relative };

struct HomologyGenerators {
  std::vector<std::vector<Halfedge>> primalGenerators;
  HomologyType primalType;
  std::vector<std::vector<Halfedge>> dualGenerators;
  HomologyType dualType;
};

// HarmonicGenerators::primalGenerators are normalized to integrate to 1 along HomologyGenerators::primalGenerators, and
// HarmonicGenerators::dualGenerators are normalized to integrate to 1 along HomologyGenerators::dualGenerators
struct HarmonicGenerators {
  std::vector<EdgeData<double>> primalGenerators;
  HomologyType primalType;
  std::vector<EdgeData<double>> dualGenerators;
  HomologyType dualType;
};

enum class HomologyGeneratorType {
  AbsolutePrimal,
  RelativePrimal,
  AbsoluteDual,
  RelativeDual,
  AbsolutePrimalRelativeDual,
  RelativePrimalAbsoluteDual
};

struct HomologyGeneratorOptions {
  HomologyGeneratorType generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
  Vertex primalRoot; // if not set, use mesh.vertex(0)
  Face dualRoot;     // if not set, use mesh.face(0)
};
// Note: primalRoot and dualRoot may not be respected when constructing relative generators on meshes with boundary

extern const HomologyGeneratorOptions defaultHomologyGeneratorOptions;

// Returns primal and/or dual homology generators as specified by opt
HomologyGenerators computeHomologyGenerators(ManifoldSurfaceMesh& mesh,
                                             HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions);

// Returns harmonic 1-forms dual to the given generators
// HarmonicGenerators::primalGenerators are normalized to integrate to 1 along HomologyGenerators::primalGenerators, and
// HarmonicGenerators::dualGenerators are normalized to integrate to 1 along HomologyGenerators::dualGenerators
HarmonicGenerators computeHarmonicGenerators(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                             const HomologyGenerators& generators, bool computePrimal = true,
                                             bool computeDual = true);

// Returns harmonic 1-forms dual to the generators specified by opt
// HarmonicGenerators::primalGenerators are normalized to integrate to 1 along HomologyGenerators::primalGenerators, and
// HarmonicGenerators::dualGenerators are normalized to integrate to 1 along HomologyGenerators::dualGenerators
HarmonicGenerators computeHarmonicGenerators(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                             HomologyGeneratorType homologyType);

//=== helpers
namespace TreeCotree {
// Return tree encoded by mapping each vertex to the halfedge pointing to its parent. If a complementary tree is passed
// in, the new tree which we build is forbidden from using edges present in the input tree. (Storing halfedges rather
// than the parent vertex directly is useful when working with delta complexes with multi-edges.) We take the convention
// that tree[vertex].twin().vertex() is the parent vertex, and tree[face].twin().face() is the parent face.
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh, const FaceData<Halfedge>* dualTree = nullptr,
                                             HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions);
FaceData<Halfedge> buildDualSpanningTree(ManifoldSurfaceMesh& mesh, const VertexData<Halfedge>* primalTree = nullptr,
                                         HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions);

bool inDualTree(Halfedge ij, const FaceData<Halfedge>& dualTree);
bool inPrimalTree(Halfedge ij, const VertexData<Halfedge>& primalTree);

std::vector<Halfedge> extractPrimalGenerator(const VertexData<Halfedge>& primalTree, Halfedge ij);
std::vector<Halfedge> extractDualGenerator(const FaceData<Halfedge>& dualTree, Halfedge ij);
std::vector<Halfedge> walkUpPrimalTree(const VertexData<Halfedge>& primalTree, Vertex start);
std::vector<Halfedge> walkUpDualTree(const FaceData<Halfedge>& dualTree, Face start);
std::vector<Halfedge> glueAndTrimPaths(const std::vector<Halfedge>& forwardPath, Halfedge ij,
                                       const std::vector<Halfedge>& backwardPath);

HomologyType primalHomologyType(HomologyGeneratorType generatorType);
HomologyType dualHomologyType(HomologyGeneratorType generatorType);

} // namespace TreeCotree

std::string to_string(HomologyType homologyType);
std::string to_string(HomologyGeneratorType homologyGeneratorType);
std::ostream& operator<<(std::ostream& os, HomologyType homologyType);
std::ostream& operator<<(std::ostream& os, HomologyGeneratorType homologyGeneratorType);

} // namespace surface
} // namespace geometrycentral
