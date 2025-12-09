#include "geometrycentral/volume/vertex_position_geometry.h"

namespace geometrycentral {
namespace volume {


// clang-format off
VertexPositionGeometry::VertexPositionGeometry(ManifoldVolumeMesh& mesh_)
    : mesh(mesh_),
  // Construct the dependency graph of managed quantities and their callbacks
  vertexIndicesQ           (&vertexIndices,          std::bind(&VertexPositionGeometry::computeVertexIndices, this),          quantities),
  edgeIndicesQ             (&edgeIndices,            std::bind(&VertexPositionGeometry::computeEdgeIndices, this),            quantities),
  faceIndicesQ             (&faceIndices,            std::bind(&VertexPositionGeometry::computeFaceIndices, this),            quantities),
  cellIndicesQ             (&cellIndices,            std::bind(&VertexPositionGeometry::computeCellIndices, this),            quantities),
  vertexCornerIndicesQ     (&vertexCornerIndices,    std::bind(&VertexPositionGeometry::computeVertexCornerIndices, this),    quantities),
  edgeCornerIndicesQ       (&edgeCornerIndices,      std::bind(&VertexPositionGeometry::computeEdgeCornerIndices, this),      quantities),
  faceCornerIndicesQ       (&faceCornerIndices,      std::bind(&VertexPositionGeometry::computeFaceCornerIndices, this),      quantities),
  faceNormalsQ             (&faceNormals,            std::bind(&VertexPositionGeometry::computeFaceNormals, this),            quantities),
  edgeLengthsQ             (&edgeLengths,            std::bind(&VertexPositionGeometry::computeEdgeLengths, this),            quantities),
  faceAreasQ               (&faceAreas,              std::bind(&VertexPositionGeometry::computeFaceAreas, this),              quantities),
  cellVolumesQ             (&cellVolumes,            std::bind(&VertexPositionGeometry::computeCellVolumes, this),            quantities),
  vertexDualVolumesQ       (&vertexDualVolumes,      std::bind(&VertexPositionGeometry::computeVertexDualVolumes, this),      quantities),
  faceDualEdgeLengthsQ     (&faceDualEdgeLengths,    std::bind(&VertexPositionGeometry::computeFaceDualEdgeLengths, this),    quantities),
  faceHodge2Q              (&faceHodge2,             std::bind(&VertexPositionGeometry::computeFaceHodge2, this),             quantities),
  faceCornerAnglesQ        (&faceCornerAngles,       std::bind(&VertexPositionGeometry::computeFaceCornerAngles, this),       quantities),
  faceCornerAngleCotansQ   (&faceCornerAngleCotans,  std::bind(&VertexPositionGeometry::computeFaceCornerAngleCotans, this),  quantities),
  dihedralAnglesQ          (&dihedralAngles,         std::bind(&VertexPositionGeometry::computeDihedralAngles, this),         quantities),
  dihedralAngleCotansQ     (&dihedralAngleCotans,    std::bind(&VertexPositionGeometry::computeDihedralAngleCotans, this),    quantities),
  faceBarycentersQ         (&faceBarycenters,        std::bind(&VertexPositionGeometry::computeFaceBarycenters, this),        quantities),
  faceCircumcentersQ       (&faceCircumcenters,      std::bind(&VertexPositionGeometry::computeFaceCircumcenters, this),      quantities),
  cellBarycentersQ         (&cellBarycenters,        std::bind(&VertexPositionGeometry::computeCellBarycenters, this),        quantities),
  cellCircumcentersQ       (&cellCircumcenters,      std::bind(&VertexPositionGeometry::computeCellCircumcenters, this),      quantities),
  cotanLaplacianQ          (&cotanLaplacian,         std::bind(&VertexPositionGeometry::computeCotanLaplacian, this),         quantities),
  vertexLumpedMassMatrixQ  (&vertexLumpedMassMatrix, std::bind(&VertexPositionGeometry::computeVertexLumpedMassMatrix, this), quantities)
  {
  }
// clang-format on

VertexPositionGeometry::VertexPositionGeometry(ManifoldVolumeMesh& mesh_,
                                               const VertexData<Vector3>& inputVertexPositions)
    : VertexPositionGeometry(mesh_) {
  vertexPositions = inputVertexPositions;
}

VertexPositionGeometry::VertexPositionGeometry(ManifoldVolumeMesh& mesh_,
                                               const std::vector<Vector3>& inputVertexPositions)
    : VertexPositionGeometry(mesh_) {
  vertexPositions = VertexData<Vector3>(mesh);
  for (Vertex i : mesh.vertices()) vertexPositions[i] = inputVertexPositions[i.getIndex()];
}

VertexPositionGeometry::~VertexPositionGeometry() {}

void VertexPositionGeometry::refreshQuantities() {
  for (DependentQuantity* q : quantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : quantities) {
    q->ensureHaveIfRequired();
  }
}

void VertexPositionGeometry::purgeQuantities() {
  for (DependentQuantity* q : quantities) {
    q->clearIfNotRequired();
  }
}

std::unique_ptr<VertexPositionGeometry> VertexPositionGeometry::reinterpretTo(ManifoldVolumeMesh& targetMesh) {
  std::unique_ptr<VertexPositionGeometry> newGeom(new VertexPositionGeometry(targetMesh));
  newGeom->vertexPositions = vertexPositions.reinterpretTo(targetMesh);
  return newGeom;
}

// === Quantity implementations

// == Indices
void VertexPositionGeometry::computeVertexIndices() { vertexIndices = mesh.getVertexIndices(); }
void VertexPositionGeometry::requireVertexIndices() { vertexIndicesQ.require(); }
void VertexPositionGeometry::unrequireVertexIndices() { vertexIndicesQ.unrequire(); }

void VertexPositionGeometry::computeEdgeIndices() { edgeIndices = mesh.getEdgeIndices(); }
void VertexPositionGeometry::requireEdgeIndices() { edgeIndicesQ.require(); }
void VertexPositionGeometry::unrequireEdgeIndices() { edgeIndicesQ.unrequire(); }

void VertexPositionGeometry::computeFaceIndices() { faceIndices = mesh.getFaceIndices(); }
void VertexPositionGeometry::requireFaceIndices() { faceIndicesQ.require(); }
void VertexPositionGeometry::unrequireFaceIndices() { faceIndicesQ.unrequire(); }

void VertexPositionGeometry::computeCellIndices() { cellIndices = mesh.getCellIndices(); }
void VertexPositionGeometry::requireCellIndices() { cellIndicesQ.require(); }
void VertexPositionGeometry::unrequireCellIndices() { cellIndicesQ.unrequire(); }

void VertexPositionGeometry::computeVertexCornerIndices() { vertexCornerIndices = mesh.getIncidenceIndices<0, 3>(); }
void VertexPositionGeometry::requireVertexCornerIndices() { vertexCornerIndicesQ.require(); }
void VertexPositionGeometry::unrequireVertexCornerIndices() { vertexCornerIndicesQ.unrequire(); }

void VertexPositionGeometry::computeEdgeCornerIndices() { edgeCornerIndices = mesh.getIncidenceIndices<1, 3>(); }
void VertexPositionGeometry::requireEdgeCornerIndices() { edgeCornerIndicesQ.require(); }
void VertexPositionGeometry::unrequireEdgeCornerIndices() { edgeCornerIndicesQ.unrequire(); }

void VertexPositionGeometry::computeFaceCornerIndices() { faceCornerIndices = mesh.getIncidenceIndices<0, 2>(); }
void VertexPositionGeometry::requireFaceCornerIndices() { faceCornerIndicesQ.require(); }
void VertexPositionGeometry::unrequireFaceCornerIndices() { faceCornerIndicesQ.unrequire(); }

// == Geometry
void VertexPositionGeometry::computeFaceNormals() {
  faceNormals = FaceData<Vector3>(mesh);

  for (Face f : mesh.faces()) {

    // Gather vertex positions for next three vertices
    Dart he = f.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == f.dart(), "faces must be triangular");

    faceNormals[f] = unit(cross(pB - pA, pC - pA));
  }
}
void VertexPositionGeometry::requireFaceNormals() { faceNormalsQ.require(); }
void VertexPositionGeometry::unrequireFaceNormals() { faceNormalsQ.unrequire(); }

void VertexPositionGeometry::computeEdgeLengths() {
  edgeLengths = EdgeData<double>(mesh);
  for (Edge e : mesh.edges()) {
    Dart he = e.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    Vector3 pB = vertexPositions[he.next().vertex()];
    edgeLengths[e] = (pB - pA).norm();
  }
}
void VertexPositionGeometry::requireEdgeLengths() { edgeLengthsQ.require(); }
void VertexPositionGeometry::unrequireEdgeLengths() { edgeLengthsQ.unrequire(); }

void VertexPositionGeometry::computeFaceAreas() {
  faceAreas = FaceData<double>(mesh);
  for (Face f : mesh.faces()) {
    // Gather vertex positions for next three vertices
    Dart he = f.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == f.dart(), "faces must be triangular");

    Vector3 N = cross(pB - pA, pC - pA);
    double area = 0.5 * norm(N);
    faceAreas[f] = area;
  }
}
void VertexPositionGeometry::requireFaceAreas() { faceAreasQ.require(); }
void VertexPositionGeometry::unrequireFaceAreas() { faceAreasQ.unrequire(); }

void VertexPositionGeometry::computeCellVolumes() {
  cellVolumes = CellData<double>(mesh);
  for (Cell f : mesh.cells()) {
    // WARNING: Logic duplicated between cached and immediate version

    Dart he = f.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == f.dart(), "2-cells must be triangular");

    he = he.partner(1).next().next();
    Vector3 pD = vertexPositions[he.vertex()];

    double det = dot(cross(pB - pA, pC - pA), pD - pA); // is the sign backwards?

    cellVolumes[f] = fabs(det) / 6.; // take absolute value to be sure of sign
  }
}
void VertexPositionGeometry::requireCellVolumes() { cellVolumesQ.require(); }
void VertexPositionGeometry::unrequireCellVolumes() { cellVolumesQ.unrequire(); }

void VertexPositionGeometry::computeVertexDualVolumes() {
  faceAreasQ.ensureHave();

  vertexDualVolumes = VertexData<double>(mesh, 0.);

  for (Cell c : mesh.cells()) {
    double V = cellVolumes[c];
    for (Vertex v : c.adjacentVertices()) {
      vertexDualVolumes[v] += V / 4.;
    }
  }
}
void VertexPositionGeometry::requireVertexDualVolumes() { vertexDualVolumesQ.require(); }
void VertexPositionGeometry::unrequireVertexDualVolumes() { vertexDualVolumesQ.unrequire(); }

void VertexPositionGeometry::computeFaceDualEdgeLengths() {
  cellCircumcentersQ.ensureHave();
  faceNormalsQ.ensureHave();

  faceDualEdgeLengths = FaceData<double>(mesh, 0.);

  for (Face f : mesh.faces()) {
    faceDualEdgeLengths[f] =
        dot(faceNormals[f], cellCircumcenters[f.dart().cell()] - cellCircumcenters[f.dart().partner(2).cell()]);
  }
}
void VertexPositionGeometry::requireFaceDualEdgeLengths() { faceDualEdgeLengthsQ.require(); }
void VertexPositionGeometry::unrequireFaceDualEdgeLengths() { faceDualEdgeLengthsQ.unrequire(); }

void VertexPositionGeometry::computeFaceHodge2() {
  faceAreasQ.ensureHave();
  faceDualEdgeLengthsQ.ensureHave();

  faceHodge2 = FaceData<double>(mesh, 0.);
  for (Face f : mesh.faces()) faceHodge2[f] = faceDualEdgeLengths[f] / faceAreas[f];
}
void VertexPositionGeometry::requireFaceHodge2() { faceHodge2Q.require(); }
void VertexPositionGeometry::unrequireFaceHodge2() { faceHodge2Q.unrequire(); }

void VertexPositionGeometry::computeFaceCornerAngles() {
  faceCornerAngles = FaceCornerData<double>(mesh);

  for (FaceCorner c : mesh.faceCorners()) {
    Dart he = c.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == c.dart(), "faces must be triangular");

    double q = dot(unit(pB - pA), unit(pC - pA));
    q = clamp(q, -1.0, 1.0);
    faceCornerAngles[c] = std::acos(q);
  }
}
void VertexPositionGeometry::requireFaceCornerAngles() { faceCornerAnglesQ.require(); }
void VertexPositionGeometry::unrequireFaceCornerAngles() { faceCornerAnglesQ.unrequire(); }

void VertexPositionGeometry::computeFaceCornerAngleCotans() {
  faceCornerAngleCotans = FaceCornerData<double>(mesh);

  for (FaceCorner c : mesh.faceCorners()) {
    Dart he = c.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == c.dart(), "faceCotans must be triangular");

    double cos = dot(pB - pA, pC - pA);
    double sin = cross(pB - pA, pC - pA).norm();
    faceCornerAngleCotans[c] = cos / sin;
  }
}
void VertexPositionGeometry::requireFaceCornerAngleCotans() { faceCornerAngleCotansQ.require(); }
void VertexPositionGeometry::unrequireFaceCornerAngleCotans() { faceCornerAngleCotansQ.unrequire(); }

void VertexPositionGeometry::computeDihedralAngles() {
  faceNormalsQ.ensureHave();

  dihedralAngles = EdgeCornerData<double>(mesh);

  for (EdgeCorner w : mesh.edgeCorners()) {
    Dart d = w.dart();
    Cell c = w.cell();
    Face f0 = d.face();
    Face f1 = d.partner(1).face(); // partner gets other face inside same 3-cell

    Vector3 n0 = faceNormals[f0] * f0.signInCell(c);
    Vector3 n1 = faceNormals[f1] * f1.signInCell(c);

    double q = dot(unit(n0), unit(n1));
    q = clamp(q, -1.0, 1.0);
    dihedralAngles[w] = M_PI - std::acos(q);
  }
}
void VertexPositionGeometry::requireDihedralAngles() { dihedralAnglesQ.require(); }
void VertexPositionGeometry::unrequireDihedralAngles() { dihedralAnglesQ.unrequire(); }

void VertexPositionGeometry::computeDihedralAngleCotans() {
  faceNormalsQ.ensureHave();

  dihedralAngleCotans = EdgeCornerData<double>(mesh);

  for (EdgeCorner w : mesh.edgeCorners()) {
    Dart d = w.dart();
    Cell c = w.cell();
    Face f0 = d.face();
    Face f1 = d.partner(1).face(); // partner gets other face inside same 3-cell

    Vector3 n0 = faceNormals[f0] * f0.signInCell(c);
    Vector3 n1 = faceNormals[f1] * f1.signInCell(c);

    double cos = -dot(n0, n1);         // cos(π-x) = -cos(x)
    double sin = cross(n0, n1).norm(); // sin(π-x) = sin(x)
    dihedralAngleCotans[w] = cos / sin;
  }
}
void VertexPositionGeometry::requireDihedralAngleCotans() { dihedralAngleCotansQ.require(); }
void VertexPositionGeometry::unrequireDihedralAngleCotans() { dihedralAngleCotansQ.unrequire(); }

void VertexPositionGeometry::computeFaceBarycenters() {
  faceBarycenters = FaceData<Vector3>(mesh);

  for (Face f : mesh.faces()) {
    Vector3 total = Vector3::zero();
    double degree = 0;
    for (Vertex i : f.adjacentVertices()) {
      total += vertexPositions[i];
      degree += 1;
    }

    faceBarycenters[f] = total / degree;
  }
}
void VertexPositionGeometry::requireFaceBarycenters() { faceBarycentersQ.require(); }
void VertexPositionGeometry::unrequireFaceBarycenters() { faceBarycentersQ.unrequire(); }

void VertexPositionGeometry::computeFaceCircumcenters() {
  faceCircumcenters = FaceData<Vector3>(mesh);

  for (Face f : mesh.faces()) {

    // Gather vertex positions for next three vertices
    Dart he = f.dart();
    Vector3 pi = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pj = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pk = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == f.dart(), "faces must be triangular");

    double lij2 = (pj - pi).norm2(), ljk2 = (pk - pj).norm2(), lki2 = (pi - pk).norm2();
    // see e.g. https://web.evanchen.cc/handouts/bary/bary-full.pdf, Appendix B.1
    Vector3 circumcenterHomog = {ljk2 * (lki2 + lij2 - ljk2), lki2 * (lij2 + ljk2 - lki2), lij2 * (ljk2 + lki2 - lij2)};
    Vector3 circumcenterBary = circumcenterHomog / (circumcenterHomog.x + circumcenterHomog.y + circumcenterHomog.z);

    faceCircumcenters[f] = circumcenterBary.x * pi + circumcenterBary.y * pj + circumcenterBary.z * pk;
  }
}
void VertexPositionGeometry::requireFaceCircumcenters() { faceCircumcentersQ.require(); }
void VertexPositionGeometry::unrequireFaceCircumcenters() { faceCircumcentersQ.unrequire(); }

void VertexPositionGeometry::computeCellBarycenters() {
  cellBarycenters = CellData<Vector3>(mesh);

  for (Cell c : mesh.cells()) {
    Vector3 total = Vector3::zero();
    double degree = 0;
    for (Vertex i : c.adjacentVertices()) {
      total += vertexPositions[i];
      degree += 1;
    }

    cellBarycenters[c] = total / degree;
  }
}
void VertexPositionGeometry::requireCellBarycenters() { cellBarycentersQ.require(); }
void VertexPositionGeometry::unrequireCellBarycenters() { cellBarycentersQ.unrequire(); }

void VertexPositionGeometry::computeCellCircumcenters() {
  cellCircumcenters = CellData<Vector3>(mesh);

  for (Cell c : mesh.cells()) {
    // Gather vertex positions for the tet's vertices
    Dart he = c.dart();
    Vector3 pA = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = vertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = vertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == c.dart(), "2-cells must be triangular");

    he = he.partner(1).next().next();
    Vector3 pD = vertexPositions[he.vertex()];

    // see, e.g. https://math.stackexchange.com/a/4481379
    // Lévy, Bruno; Liu, Yang (2010). "Lp Centroidal Voronoi Tessellation and
    // its applications". ACM: 119.
    Eigen::Matrix3d M;
    M.row(0) = Eigen::Vector3d(pB - pA);
    M.row(1) = Eigen::Vector3d(pC - pA);
    M.row(2) = Eigen::Vector3d(pD - pA);
    Eigen::Vector3d rhs;
    double a2 = pA.norm2(), b2 = pB.norm2(), c2 = pC.norm2(), d2 = pD.norm2();
    rhs << b2 - a2, c2 - a2, d2 - a2;
    rhs /= 2.;
    cellCircumcenters[c] = Vector3::fromEigen(M.colPivHouseholderQr().solve(rhs));
  }
}
void VertexPositionGeometry::requireCellCircumcenters() { cellCircumcentersQ.require(); }
void VertexPositionGeometry::unrequireCellCircumcenters() { cellCircumcentersQ.unrequire(); }

void VertexPositionGeometry::computeCotanLaplacian() {
  vertexIndicesQ.ensureHave();
  edgeLengthsQ.ensureHave();
  dihedralAngleCotansQ.ensureHave();
  std::vector<Eigen::Triplet<double>> triplets;
  for (EdgeCorner w : mesh.edgeCorners()) {
    Dart d = w.dart();
    Dart opp = d.next().partner(1).next();
    double lOpp = edgeLengths[opp.edge()];
    double cotAngleOpp = dihedralAngleCotans[opp.edgeCorner()];
    double weight = lOpp * cotAngleOpp / 6.;
    size_t i = vertexIndices[d.tailVertex()], j = vertexIndices[d.tipVertex()];
    triplets.emplace_back(i, i, weight);
    triplets.emplace_back(i, j, -weight);
    triplets.emplace_back(j, i, -weight);
    triplets.emplace_back(j, j, weight);
  }
  cotanLaplacian = SparseMatrix<double>(mesh.nVertices(), mesh.nVertices());
  cotanLaplacian.setFromTriplets(triplets.begin(), triplets.end());
}
void VertexPositionGeometry::requireCotanLaplacian() { cotanLaplacianQ.require(); }
void VertexPositionGeometry::unrequireCotanLaplacian() { cotanLaplacianQ.unrequire(); }

void VertexPositionGeometry::computeVertexLumpedMassMatrix() {
  vertexIndicesQ.ensureHave();
  vertexDualVolumesQ.ensureHave();
  Eigen::VectorXd hodge0V(mesh.nVertices());
  for (Vertex i : mesh.vertices()) hodge0V[vertexIndices[i]] = vertexDualVolumes[i];
  vertexLumpedMassMatrix = hodge0V.asDiagonal();
}
void VertexPositionGeometry::requireVertexLumpedMassMatrix() { vertexLumpedMassMatrixQ.require(); }
void VertexPositionGeometry::unrequireVertexLumpedMassMatrix() { vertexLumpedMassMatrixQ.unrequire(); }

// === Immediate computations
Vector3 VertexPositionGeometry::gradient(const VertexData<double>& u, Cell c) {
  cellVolumesQ.ensureHave();
  faceAreasQ.ensureHave();
  faceNormalsQ.ensureHave();

  Vector3 result = Vector3::zero();
  double volume = cellVolumes[c];
  for (Vertex i : c.adjacentVertices()) {
    Dart d = i.dartInCell(c);
    Face oppFace = d.next().partner(1).face(); // opposite face
    Vector3 n = faceNormals[oppFace] * oppFace.signInCell(c);
    double area = faceAreas[oppFace];

    result += u[i] * area * n / (3. * volume);
  }
  return result;
}

CellData<Vector3> VertexPositionGeometry::gradient(const VertexData<double>& u) {
  CellData<Vector3> result(mesh);
  for (Cell c : mesh.cells()) result[c] = gradient(u, c);
  return result;
}

Vector3 VertexPositionGeometry::faceGradient(const VertexData<double>& u, Face f) {
  const VertexData<Vector3>& p = vertexPositions; // shorter alias
  Vector3 result = Vector3::zero();
  Dart d = f.dart();
  Vector3 areaNormal = cross(p[d.tipVertex()] - p[d.tailVertex()], p[d.next().tipVertex()] - p[d.tailVertex()]);
  do { // explicit loop to ensure correct orientation
    result += u[d.next().tipVertex()] * (p[d.tipVertex()] - p[d.tailVertex()]);
    d = d.next();
  } while (d != f.dart());
  return cross(areaNormal, result) / areaNormal.norm2();
}

FaceData<Vector3> VertexPositionGeometry::faceGradient(const VertexData<double>& u) {
  FaceData<Vector3> result(mesh);
  for (Face f : mesh.faces()) result[f] = faceGradient(u, f);
  return result;
}

} // namespace volume
} // namespace geometrycentral
