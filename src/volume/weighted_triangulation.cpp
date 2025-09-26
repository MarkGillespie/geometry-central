#include "geometrycentral/volume/weighted_triangulation.h"

#include <Eigen/Sparse>

namespace geometrycentral {
namespace volume {
Vector3 weightedEdgeCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Edge e) {
  Vector3 pi = geom.vertexPositions[e.dart().tailVertex()], pj = geom.vertexPositions[e.dart().tipVertex()];
  double wi = weights[e.dart().tailVertex()], wj = weights[e.dart().tipVertex()];

  return (pj + pi) / 2. - (wj - wi) * (pj - pi) / (2. * (pj - pi).norm2());
}

Vector3 weightedFaceCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Face f) {
  geom.requireFaceCircumcenters();
  Vector3 result = geom.faceCircumcenters[f] - geom.faceGradient(weights, f) / 2.;
  geom.unrequireFaceCircumcenters();
  return result;
}

Vector3 weightedCellCenter(VertexPositionGeometry& geom, const VertexData<double>& weights, Cell c) {
  geom.requireCellCircumcenters();
  Vector3 result = geom.cellCircumcenters[c] - geom.gradient(weights, c) / 2.;
  geom.unrequireCellCircumcenters();
  return result;
}

EdgeData<Vector3> weightedEdgeCenters(VertexPositionGeometry& geom, const VertexData<double>& weights) {
  EdgeData<Vector3> result(geom.mesh);
  for (Edge e : geom.mesh.edges()) result[e] = weightedEdgeCenter(geom, weights, e);
  return result;
}

FaceData<Vector3> weightedFaceCenters(VertexPositionGeometry& geom, const VertexData<double>& weights) {
  FaceData<Vector3> result(geom.mesh);
  for (Face f : geom.mesh.faces()) result[f] = weightedFaceCenter(geom, weights, f);
  return result;
}

CellData<Vector3> weightedCellCenters(VertexPositionGeometry& geom, const VertexData<double>& weights) {
  CellData<Vector3> result(geom.mesh);
  for (Cell c : geom.mesh.cells()) result[c] = weightedCellCenter(geom, weights, c);
  return result;
}

// if grad != nullptr, the gradient (with respect to the weights) is added to its contents
double weightedEdgeDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Vertex v, Edge e,
                        VertexData<double>* grad) {
  Dart d = e.dartInCell(v);
  Vector3 pi = geom.vertexPositions[d.tailVertex()], pj = geom.vertexPositions[d.tipVertex()];
  double wi = weights[d.tailVertex()], wj = weights[d.tipVertex()];
  double lij2 = (pj - pi).norm2();
  double lij = sqrt(lij2);

  double result = (lij2 - (wj - wi)) / (2. * lij);
  if (grad) {
    (*grad)[d.tailVertex()] += 1. / (2. * lij); // d/di
    (*grad)[d.tipVertex()] -= 1. / (2. * lij);  // d/dj
  }
  return result;
}

double weightedFaceDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Edge e, Face f,
                        VertexData<double>* grad) {
  geom.requireEdgeLengths();
  geom.requireFaceAreas();
  geom.requireFaceCornerAngleCotans();
  Dart d = f.dartInCell(e);
  // we index the vertices so that edge e goes from vertex i -> vertex j, opposite vertex k in face f
  double wi = weights[d.tailVertex()], wj = weights[d.tipVertex()], wk = weights[d.next().tipVertex()];
  double cotbi = geom.faceCornerAngleCotans[d.faceCorner()], cotbj = geom.faceCornerAngleCotans[d.next().faceCorner()],
         cotbk = geom.faceCornerAngleCotans[d.next().next().faceCorner()];
  double lij = geom.edgeLengths[d.edge()], Aijk = geom.faceAreas[f];
  geom.unrequireFaceCornerAngleCotans();
  geom.unrequireFaceAreas();
  geom.unrequireEdgeLengths();
  double result = cotbk * lij / 2. + (cotbi * wj + cotbj * wi) / (2. * lij) - wk * lij / (4. * Aijk);
  if (grad) {
    (*grad)[d.tailVertex()] += cotbj / (2. * lij);      // d/di
    (*grad)[d.tipVertex()] += cotbi / (2. * lij);       // d/dj
    (*grad)[d.next().tipVertex()] -= lij / (4. * Aijk); // d/dk
  }
  return result;
}

double weightedCellDist(VertexPositionGeometry& geom, const VertexData<double>& weights, Face f, Cell c,
                        VertexData<double>* grad) {
  Dart d = f.dartInCell(c);
  Vector3 pi = geom.vertexPositions[d.tailVertex()], pj = geom.vertexPositions[d.tipVertex()],
          pk = geom.vertexPositions[d.next().tipVertex()], pl = geom.vertexPositions[d.partner(1).next().tipVertex()];

  Vector3 n = unit(cross(pj - pi, pk - pi));
  Vector3 inwardNormal = dot(n, pl - pi) > 0 ? n : -n;

  double result = dot(inwardNormal, weightedCellCenter(geom, weights, c) - weightedFaceCenter(geom, weights, f));
  if (grad) {
    // Since the weightedFaceCenter always stays on the face, the gradient is equal to
    // the gradient of dot(inwardNormal, weightedCellCenter). And since weightedCellCenter is equal to -.5 ∇w plus the
    // circumcenter (which is also independent of the weights), we just have to compute the derivative of
    // dot(inwardNormal, ∇w )
    geom.requireFaceNormals();
    geom.requireFaceAreas();
    geom.requireCellVolumes();
    double volume = geom.cellVolumes[c];
    for (Vertex i : c.adjacentVertices()) {
      Dart d = i.dartInCell(c);
      Face oppFace = d.next().partner(1).face(); // opposite face
      Vector3 n = geom.faceNormals[oppFace] * oppFace.signInCell(c);
      double area = geom.faceAreas[oppFace];

      (*grad)[i] -= .5 * area * dot(inwardNormal, n) / (3. * volume);
    }
    geom.unrequireCellVolumes();
    geom.unrequireFaceAreas();
    geom.unrequireFaceNormals();
  }
  return result;
}

//==== Hodge-Optimized Triangulation energies

// Implement HOT-2,2 energy with given denominators. Varying the denominators produces the energy for different Hodge
// stars. If grad != nullptr, adds the gradient (with respect to the weights) to grad
double HOT22EnergyHelper(VertexPositionGeometry& geom, const VertexData<double>& weights, double k0, double k1,
                         double k2, double k3, VertexData<double>* grad) {
  double result = 0;
  VertexData<double> gHf, ghe, gdv;
  for (Cell t : geom.mesh.cells()) {
    for (Face f : t.adjacentFaces()) {
      if (grad) gHf = VertexData<double>(geom.mesh, 0); // if we're computing gradients, save grad of Hf here
      double Hf = weightedCellDist(geom, weights, f, t, grad ? &gHf : nullptr);
      double Hf3 = pow(Hf, 3);
      for (Edge e : f.adjacentEdges()) {
        if (grad) ghe = VertexData<double>(geom.mesh, 0); // if we're computing gradients, save grad of He here
        double he = weightedFaceDist(geom, weights, e, f, grad ? &ghe : nullptr);
        double he3 = pow(he, 3);
        for (Vertex v : e.adjacentVertices()) {
          if (grad) gdv = VertexData<double>(geom.mesh, 0); // if we're computing gradients, save grad of dv here
          double dv = weightedEdgeDist(geom, weights, v, e, grad ? &gdv : nullptr);
          double dv3 = pow(dv, 3);
          result += (Hf3 * he * dv / k0 + Hf * he3 * dv / k1 + Hf * he * dv3 / k2) / k3;

          if (grad) { // if necessary, add gradient with respect to w to *grad
            *grad += (3 * pow(Hf, 2) * he * dv / k0 + he3 * dv / k1 + he * dv3 / k2) / k3 * gHf;
            *grad += (Hf3 * dv / k0 + 3 * Hf * pow(he, 2) * dv / k1 + Hf * dv3 / k2) / k3 * ghe;
            *grad += (Hf3 * he / k0 + Hf * he3 / k1 + 3 * Hf * he * pow(dv, 2) / k2) / k3 * gdv;
          }
        }
      }
    }
  }
  return result;
}

double star0HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights, VertexData<double>* grad) {
  double k0 = 12., k1 = 4., k2 = 2., k3 = 5. / 6.;
  return HOT22EnergyHelper(geom, weights, k0, k1, k2, k3, grad);
}

double star1HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights, VertexData<double>* grad) {
  double k0 = 12., k1 = 4., k2 = 6., k3 = 3. / 6.;
  return HOT22EnergyHelper(geom, weights, k0, k1, k2, k3, grad);
}

double star2HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights, VertexData<double>* grad) {
  double k0 = 6., k1 = 4., k2 = 12., k3 = 3. / 6.;
  return HOT22EnergyHelper(geom, weights, k0, k1, k2, k3, grad);
}

double star3HOT22Energy(VertexPositionGeometry& geom, const VertexData<double>& weights, VertexData<double>* grad) {
  double k0 = 2., k1 = 4., k2 = 12., k3 = 5. / 6.;
  return HOT22EnergyHelper(geom, weights, k0, k1, k2, k3, grad);
}

VertexData<double> optimalStar3Hot22Weights(VertexPositionGeometry& geom) {
  // The gradient for *3-HOT_{2,2} is a linear function of the weights, so they can be found by solving a linear system
  // See Appendix A of Hodge-optimized triangulations by Mullen et al. 2011, and Section 6.2 of Weighted Triangulations
  // for Geometry Processing by de Goes et al. 2014

  // Solve Lw = div(circumcenter - barycenter), where
  // div(circumcenter - barycenter)_i = \sum_{tets ijkl} dot(circumcenter_ijkl - barycenter_ijkl, A_jkl n_jkl)

  geom.requireCellCircumcenters();
  geom.requireCellBarycenters();
  geom.requireCellVolumes();
  geom.requireFaceNormals();
  geom.requireFaceAreas();
  geom.requireVertexIndices();
  const VertexData<size_t>& vIdx = geom.vertexIndices; // define shorter aliases for geom quantities
  const FaceData<Vector3>& n = geom.faceNormals;
  const FaceData<double>& A = geom.faceAreas;
  const CellData<Vector3>& cc = geom.cellCircumcenters;
  const CellData<Vector3>& bc = geom.cellBarycenters;
  const CellData<double>& vol = geom.cellVolumes;
  Vector<double> rhs = Vector<double>::Zero(geom.mesh.nVertices());
  for (Vertex i : geom.mesh.vertices()) {
    for (Cell t : i.adjacentCells()) {
      Dart d = t.dartInCell(i);
      Face fOpp = d.next().partner(1).face();
      rhs(vIdx[i]) += .5 * dot(cc[t] - bc[t], fOpp.signInCell(t) * n[fOpp] * A[fOpp]) / (3. * vol[t]);
    }
  }
  geom.unrequireVertexIndices();
  geom.unrequireFaceAreas();
  geom.unrequireFaceNormals();
  geom.unrequireCellVolumes();
  geom.unrequireCellBarycenters();
  geom.unrequireCellCircumcenters();

  geom.requireCotanLaplacian();

  // check that Lw - rhs = gradient of HOT energy
  VertexData<double> hotGrad(geom.mesh, 0);
  VertexData<double> w(geom.mesh, 0);
  VertexData<double> Lw(geom.mesh, geom.cotanLaplacian * w.raw());
  star3HOT22Energy(geom, w, &hotGrad);
  double err = (Lw.raw() - rhs - hotGrad.raw()).norm();
  if (err > 1e-5) {
    std::cout << "gradient error." << std::endl;
    for (size_t iV = 0; iV < 9; iV++) {
      std::cout << "  grad(" << iV << ") = " << hotGrad[iV] << std::endl;
      std::cout << "  -rhs(" << iV << ") = " << Lw[iV] - rhs(iV) << std::endl;
    }
  }

  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::DiagonalPreconditioner<double>> solver;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
  solver.setTolerance(1e-8);
  solver.setMaxIterations(1000);
  solver.compute(geom.cotanLaplacian);
  geom.unrequireCotanLaplacian();
  // Eigen::VectorXd x = solver.solve(rhs);
  return VertexData<double>(geom.mesh, solver.solve(rhs));
}

} // namespace volume
} // namespace geometrycentral
