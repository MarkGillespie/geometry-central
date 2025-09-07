#include "geometrycentral/surface/homology_generators.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h" // blockDecomposeSquare


namespace geometrycentral {
namespace surface {
// The default options
const HomologyGeneratorOptions defaultHomologyGeneratorOptions;

HomologyGenerators computeHomologyGenerators(ManifoldSurfaceMesh& mesh, HomologyGeneratorOptions opt) {
  using namespace TreeCotree;

  HomologyGenerators result;

  if (opt.generatorType == HomologyGeneratorType::AbsolutePrimal ||
      opt.generatorType == HomologyGeneratorType::AbsolutePrimalRelativeDual ||
      opt.generatorType == HomologyGeneratorType::RelativeDual) {
    result.primalType = HomologyType::Absolute;
    result.dualType = HomologyType::Relative;

    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh, nullptr, opt);
    FaceData<Halfedge> dualTree = buildDualSpanningTree(mesh, &primalTree, opt);

    for (Edge e : mesh.edges()) {
      Halfedge ij = e.halfedge();
      if (!inPrimalTree(ij, primalTree) && !inDualTree(ij, dualTree) && !e.isBoundary()) {
        if (opt.generatorType != HomologyGeneratorType::RelativeDual) // if requested, store primal generator
          result.primalGenerators.push_back(extractPrimalGenerator(primalTree, ij));
        if (opt.generatorType != HomologyGeneratorType::AbsolutePrimal) // if requested, store dual generator
          result.dualGenerators.push_back(extractDualGenerator(dualTree, ij));
      }
    }
  } else if (opt.generatorType == HomologyGeneratorType::AbsoluteDual ||
             opt.generatorType == HomologyGeneratorType::RelativePrimalAbsoluteDual ||
             opt.generatorType == HomologyGeneratorType::RelativePrimal) {
    result.dualType = HomologyType::Absolute;
    result.primalType = HomologyType::Relative;

    FaceData<Halfedge> dualTree = buildDualSpanningTree(mesh, nullptr, opt);
    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh, &dualTree, opt);

    for (Edge e : mesh.edges()) {
      Halfedge ij = e.halfedge();
      if (!inDualTree(ij, dualTree) && !inPrimalTree(ij, primalTree) && !e.isBoundary()) {
        if (opt.generatorType != HomologyGeneratorType::RelativePrimal) // if requested, store dual generator
          result.dualGenerators.push_back(extractDualGenerator(dualTree, ij));
        if (opt.generatorType != HomologyGeneratorType::AbsoluteDual) // if requested, store primal generator
          result.primalGenerators.push_back(extractPrimalGenerator(primalTree, ij));
      }
    }
  }

  return result;
}

// Returns harmonic 1-forms dual to the given generators
HarmonicGenerators computeHarmonicGenerators(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                             const HomologyGenerators& generators, bool computePrimal,
                                             bool computeDual) {
  auto sign = [&](Halfedge ij) -> double { return ij.orientation() ? 1. : -1.; };
  HarmonicGenerators result;
  geom.requireCotanLaplacian();
  geom.requireDECOperators();
  SparseMatrix<double> L0 = geom.cotanLaplacian;
  const SparseMatrix<double>&d0 = geom.d0, &d1 = geom.d1;
  const SparseMatrix<double>&hodge1 = geom.hodge1, hodge1Inv = geom.hodge1Inverse;
  SparseMatrix<double> L2 = d1 * hodge1Inv * d1.transpose();

  std::vector<Eigen::Triplet<double>> d1DirichletTriplets, d1NeumannTriplets;
  geom.requireFaceIndices();
  geom.requireEdgeIndices();
  const FaceData<size_t>& fIdx = geom.faceIndices;
  const EdgeData<size_t>& eIdx = geom.edgeIndices;
  for (Face f : mesh.faces()) {
    size_t iF = fIdx[f];
    for (Halfedge h : f.adjacentHalfedges()) {
      size_t iE = eIdx[h.edge()];
      if (h.edge().isBoundary()) {
        d1DirichletTriplets.emplace_back(iF, iE, 1);
        // d1NeumannTriplets.emplace_back(iF, iE, 0); // don't need to set zero coefficient
      } else {
        d1DirichletTriplets.emplace_back(iF, iE, sign(h));
        d1NeumannTriplets.emplace_back(iF, iE, sign(h));
      }
    }
  }
  geom.unrequireEdgeIndices();
  geom.unrequireFaceIndices();

  size_t nE = mesh.nEdges(), nF = mesh.nFaces();
  SparseMatrix<double> d1Dirichlet(nF, nE), d1Neumann(nF, nE);
  d1Dirichlet.setFromTriplets(d1DirichletTriplets.begin(), d1DirichletTriplets.end());
  d1Neumann.setFromTriplets(d1NeumannTriplets.begin(), d1NeumannTriplets.end());

  SparseMatrix<double> L2Dirichlet = d1Dirichlet * hodge1Inv * d1Dirichlet.transpose();
  SparseMatrix<double> L2Neumann = d1Neumann * hodge1Inv * d1Neumann.transpose();

  VertexData<bool> isInteriorVertex(mesh, true);
  for (BoundaryLoop b : mesh.boundaryLoops()) {
    for (Vertex i : b.adjacentVertices()) isInteriorVertex[i] = false;
  }

  BlockDecompositionResult<double> decomp0 = blockDecomposeSquare(L0, isInteriorVertex.raw(), false);
  SparseMatrix<double> L0ii = decomp0.AA;
  const SparseMatrix<double>& L0ib = decomp0.AB;

  if (computePrimal) { //===== get primal harmonic generators by solving for jump across dual generators
    result.primalGenerators.reserve(generators.dualGenerators.size());
    result.primalType =
        (generators.primalType == HomologyType::Absolute) ? HomologyType::Relative : HomologyType::Absolute;
    for (const std::vector<Halfedge>& dualGenerator : generators.dualGenerators) {
      EdgeData<double> jump(mesh, 0);
      for (Halfedge ij : dualGenerator) jump[ij.edge()] += sign(ij);

      // Solve for a jump-harmonic function w/ given jump and appropriate boundary conditions
      Vector<double> alpha;
      switch (generators.primalType) {
      case HomologyType::Absolute: { // impose zero-Neumann boundary condition on potential
        Vector<double> rhs = d0.transpose() * hodge1 * jump.raw();
        alpha = solvePositiveDefinite(L0, rhs);
        break;
      }
      case HomologyType::Relative: { // impose zero-Dirichlet boundary condition on potential
        Vector<double> rhs = d0.transpose() * hodge1 * jump.raw();
        alpha = solvePositiveDefinite(L0, rhs);

        Vector<double> fullRHS = d0.transpose() * hodge1 * jump.raw(), iRHS, bRHS;
        decomposeVector(decomp0, fullRHS, iRHS, bRHS);

        Vector<double> iPotential = solvePositiveDefinite(L0ii, iRHS);
        Vector<double> bPotential = Vector<double>::Zero(bRHS.size());
        alpha = reassembleVector(decomp0, iPotential, bPotential);
        break;
      }
      }

      EdgeData<double> gamma(mesh, d0 * alpha);
      for (Halfedge ij : dualGenerator) gamma[ij.edge()] -= sign(ij);
      result.primalGenerators.push_back(gamma);
    }
  }

  if (computeDual) { //===== get dual harmonic generators by solving for jump across primal generators
    result.dualGenerators.reserve(generators.primalGenerators.size());
    result.dualType = (generators.dualType == HomologyType::Absolute) ? HomologyType::Relative : HomologyType::Absolute;
    for (const std::vector<Halfedge>& primalGenerator : generators.primalGenerators) {
      EdgeData<double> jump(mesh, 0);
      for (Halfedge ij : primalGenerator) jump[ij.edge()] += sign(ij);

      // Solve for a jump-harmonic function w/ given jump and appropriate boundary conditions
      EdgeData<double> gamma;
      switch (generators.dualType) {
      case HomologyType::Absolute: { // impose zero-Neumann boundary condition on potential
        Vector<double> rhs = d1Neumann * hodge1Inv * jump.raw();
        Vector<double> beta = solvePositiveDefinite(L2Neumann, rhs);
        gamma = EdgeData<double>(mesh, d1Neumann.transpose() * beta);
        break;
      }
      case HomologyType::Relative: { // impose zero-Dirichlet boundary condition on potential
        Vector<double> rhs = d1Dirichlet * hodge1Inv * jump.raw();
        Vector<double> beta = solvePositiveDefinite(L2Dirichlet, rhs);
        gamma = EdgeData<double>(mesh, d1Dirichlet.transpose() * beta);
        break;
      }
      }

      for (Halfedge ij : primalGenerator) gamma[ij.edge()] -= sign(ij);
      result.dualGenerators.push_back(gamma);
    }
  }

  geom.unrequireDECOperators();
  geom.unrequireCotanLaplacian();
  return result;
}

// Returns harmonic 1-forms dual to the generators specified by opt
HarmonicGenerators computeHarmonicGenerators(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
                                             HomologyGeneratorType homologyType) {
  // Set options for homology generator loops. Note that even if we only want to compute one set of harmonic generators,
  // the algorithm requires both sets of homology generators
  bool computePrimal, computeDual;
  HomologyGeneratorOptions opt;
  switch (homologyType) {
  case HomologyGeneratorType::AbsolutePrimalRelativeDual:
    computePrimal = true;
    computeDual = true;
    opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
    break;
  case HomologyGeneratorType::AbsolutePrimal:
    computePrimal = true;
    computeDual = false;
    opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
    break;
  case HomologyGeneratorType::RelativeDual:
    computePrimal = false;
    computeDual = true;
    opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
    break;
  case HomologyGeneratorType::RelativePrimalAbsoluteDual:
    computePrimal = true;
    computeDual = true;
    opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
    break;
  case HomologyGeneratorType::RelativePrimal:
    computePrimal = true;
    computeDual = false;
    opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
    break;
  case HomologyGeneratorType::AbsoluteDual:
    computePrimal = false;
    computeDual = true;
    opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
    break;
  }
  HomologyGenerators generators = computeHomologyGenerators(mesh, opt);
  return computeHarmonicGenerators(mesh, geom, generators, computePrimal, computeDual);
}

namespace TreeCotree {
// return tree encoded by mapping each vertex to the halfedge pointing to its parent if a
// complementary tree is passed in, the new tree which we build is forbidden from using edges
// present in the input tree.
// storing halfedges rather than the parent vertex directly is useful when working with delta
// complexes with multi-edges.
// we take the convention that tree[vertex].twin().vertex() is the parent vertex, and
// tree[face].twin().face() is the parent face
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh, const FaceData<Halfedge>* dualTree,
                                             HomologyGeneratorOptions opt) {
  HomologyType homologyType = primalHomologyType(opt.generatorType);

  VertexData<Halfedge> primalTree(mesh, Halfedge());
  VertexData<bool> visited(mesh, false);

  auto inDualTreePtr = [&](Halfedge ij) -> bool { return dualTree && inDualTree(ij, *dualTree); };

  std::deque<Vertex> toVisit;
  if (homologyType == HomologyType::Relative && mesh.hasBoundary()) {
    // if we're looking for relative generators, we should connect together all boundary vertices.
    // The easiest way to do that is just to start by pushing all boundary vertices onto the queue
    for (BoundaryLoop b : mesh.boundaryLoops()) {
      for (Vertex v : b.adjacentVertices()) {
        toVisit.push_back(v);
        visited[v] = true;
      }
    }
  } else {
    // if we want absolute generators, we should just pick an arbitrary root vertex and start there
    Vertex root = (opt.primalRoot == Vertex()) ? mesh.vertex(0) : opt.primalRoot;
    toVisit.push_back(root);
    visited[root] = true;
  }

  while (!toVisit.empty()) {
    Vertex i = toVisit.front();
    toVisit.pop_front();
    for (Halfedge ji : i.incomingHalfedges()) {
      Vertex j = ji.tailVertex();
      if (!inDualTreePtr(ji) && !visited[j]) {
        primalTree[j] = ji;
        toVisit.push_back(j);
        visited[j] = true;
      }
    }
  }
  return primalTree;
}

// we take the convention that tree[vertex].twin().vertex() is the parent vertex, and
// tree[face].twin().face() is the parent face
FaceData<Halfedge> buildDualSpanningTree(ManifoldSurfaceMesh& mesh, const VertexData<Halfedge>* primalTree,
                                         HomologyGeneratorOptions opt) {
  HomologyType homologyType = dualHomologyType(opt.generatorType);

  FaceData<Halfedge> dualTree(mesh, Halfedge());
  FaceData<bool> visited(mesh, false);

  auto inPrimalTreePtr = [&](Halfedge ij) -> bool { return primalTree && inPrimalTree(ij, *primalTree); };

  std::deque<Face> toVisit;
  if (homologyType == HomologyType::Relative && mesh.hasBoundary()) {
    // if we're looking for relative generators, we should connect together all boundary faces. The
    // easiest way to do that is just to start by pushing all boundary faces onto the queue
    for (BoundaryLoop b : mesh.boundaryLoops()) {
      for (Halfedge ij : b.adjacentHalfedges()) {
        if (inPrimalTreePtr(ij)) continue;
        toVisit.push_back(ij.twin().face());
        visited[ij.twin().face()] = true;
        dualTree[ij.twin().face()] = ij.twin(); // mark boundary edges as used in dual tree
      }
    }
  } else {
    // if we want absolute generators, we should just pick an arbitrary root vertex and start there
    Face root = (opt.dualRoot == Face()) ? mesh.face(0) : opt.dualRoot;
    toVisit.push_back(root);
    visited[root] = true;
  }

  while (!toVisit.empty()) {
    Face i = toVisit.front();
    toVisit.pop_front();
    for (Halfedge ij : i.adjacentHalfedges()) {
      Face j = ij.twin().face();
      if (!j.isBoundaryLoop() && !inPrimalTreePtr(ij) && !visited[j]) {
        dualTree[j] = ij.twin();
        toVisit.push_back(j);
        visited[j] = true;
      }
    }
  }

  return dualTree;
}

bool inPrimalTree(Halfedge ij, const VertexData<Halfedge>& primalTree) {
  return primalTree[ij.tailVertex()] == ij || primalTree[ij.tipVertex()] == ij.twin();
}
bool inDualTree(Halfedge ij, const FaceData<Halfedge>& dualTree) {
  return (!ij.face().isBoundaryLoop() && dualTree[ij.face()] == ij) ||
         (!ij.twin().face().isBoundaryLoop() && dualTree[ij.twin().face()] == ij.twin());
}

std::vector<Halfedge> extractPrimalGenerator(const VertexData<Halfedge>& primalTree, Halfedge ij) {
  std::vector<Halfedge> forwardPath = walkUpPrimalTree(primalTree, ij.tipVertex());
  std::vector<Halfedge> backwardPath = walkUpPrimalTree(primalTree, ij.tailVertex());
  return glueAndTrimPaths(forwardPath, ij, backwardPath);
}

std::vector<Halfedge> extractDualGenerator(const FaceData<Halfedge>& dualTree, Halfedge ij) {
  std::vector<Halfedge> forwardPath = walkUpDualTree(dualTree, ij.face());
  std::vector<Halfedge> backwardPath = walkUpDualTree(dualTree, ij.twin().face());
  return glueAndTrimPaths(forwardPath, ij.twin(), backwardPath);
}

std::vector<Halfedge> walkUpPrimalTree(const VertexData<Halfedge>& primalTree, Vertex start) {
  std::vector<Halfedge> path;
  Vertex curr = start;
  while (true) {
    if (primalTree[curr] == Halfedge()) {
      break;
    } else {
      path.push_back(primalTree[curr]);
      curr = primalTree[curr].tipVertex();
    }
  }
  return path;
}

std::vector<Halfedge> walkUpDualTree(const FaceData<Halfedge>& dualTree, Face start) {
  std::vector<Halfedge> path;
  Face curr = start;
  while (true) {
    if (curr.isBoundaryLoop() || dualTree[curr] == Halfedge()) {
      break;
    } else {
      path.push_back(dualTree[curr]);
      curr = dualTree[curr].twin().face();
    }
  }
  return path;
}

std::vector<Halfedge> glueAndTrimPaths(const std::vector<Halfedge>& forwardPath, Halfedge ij,
                                       const std::vector<Halfedge>& backwardPath) {
  // trim shared halfedges off ends of paths
  int m = forwardPath.size() - 1, n = backwardPath.size() - 1;
  while (m >= 0 && n >= 0 && forwardPath[m] == backwardPath[n]) m--, n--;

  // combine paths into gluedPath
  std::vector<Halfedge> gluedPath;
  for (int i = n; i >= 0; i--) gluedPath.push_back(backwardPath[i].twin());
  gluedPath.push_back(ij);
  for (int i = 0; i <= m; i++) gluedPath.push_back(forwardPath[i]);
  return gluedPath;
}

HomologyType primalHomologyType(HomologyGeneratorType generatorType) {
  switch (generatorType) {
  case HomologyGeneratorType::AbsolutePrimal:
  case HomologyGeneratorType::RelativeDual:
  case HomologyGeneratorType::AbsolutePrimalRelativeDual:
    return HomologyType::Absolute;
  case HomologyGeneratorType::RelativePrimal:
  case HomologyGeneratorType::AbsoluteDual:
  case HomologyGeneratorType::RelativePrimalAbsoluteDual:
    return HomologyType::Relative;
  }
  return HomologyType::Absolute; // should not be reachable
}

HomologyType dualHomologyType(HomologyGeneratorType generatorType) {
  switch (generatorType) {
  case HomologyGeneratorType::RelativePrimal:
  case HomologyGeneratorType::AbsoluteDual:
  case HomologyGeneratorType::RelativePrimalAbsoluteDual:
    return HomologyType::Absolute;
  case HomologyGeneratorType::AbsolutePrimal:
  case HomologyGeneratorType::RelativeDual:
  case HomologyGeneratorType::AbsolutePrimalRelativeDual:
    return HomologyType::Relative;
  }
  return HomologyType::Absolute; // should not be reachable
}
} // namespace TreeCotree

std::string to_string(HomologyType homologyType) {
  switch (homologyType) {
  case HomologyType::Absolute:
    return "HomologyType::Absolute";
  case HomologyType::Relative:
    return "HomologyType::Relative";
  }
  return ""; // should not be reachable
}

std::string to_string(HomologyGeneratorType homologyGeneratorType) {
  switch (homologyGeneratorType) {
  case HomologyGeneratorType::AbsolutePrimalRelativeDual:
    return "HomologyGeneratorType::AbsolutePrimalRelativeDual";
  case HomologyGeneratorType::AbsolutePrimal:
    return "HomologyGeneratorType::AbsolutePrimal";
  case HomologyGeneratorType::RelativeDual:
    return "HomologyGeneratorType::RelativeDual";
  case HomologyGeneratorType::RelativePrimalAbsoluteDual:
    return "HomologyGeneratorType::RelativePrimalAbsoluteDual";
  case HomologyGeneratorType::RelativePrimal:
    return "HomologyGeneratorType::RelativePrimal";
  case HomologyGeneratorType::AbsoluteDual:
    return "HomologyGeneratorType::AbsoluteDual";
  }
  return ""; // should not be reachable
}

std::ostream& operator<<(std::ostream& os, HomologyType homologyType) {
  os << to_string(homologyType);
  return os;
}

std::ostream& operator<<(std::ostream& os, HomologyGeneratorType homologyGeneratorType) {
  os << to_string(homologyGeneratorType);
  return os;
}

} // namespace surface
} // namespace geometrycentral
