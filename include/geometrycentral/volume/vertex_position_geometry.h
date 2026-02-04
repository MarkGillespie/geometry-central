#pragma once

#include "geometrycentral/utilities/dependent_quantity.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"
#include "geometrycentral/volume/manifold_volume_mesh.h"

namespace geometrycentral {
namespace volume {

class VertexPositionGeometry {

public:
  VertexPositionGeometry(ManifoldVolumeMesh& mesh);

  // Construct from positions
  VertexPositionGeometry(ManifoldVolumeMesh& mesh_, const VertexData<Vector3>& inputVertexPositions);
  VertexPositionGeometry(ManifoldVolumeMesh& mesh_, const std::vector<Vector3>& inputVertexPositions);

  // Construct from positions (stored in an Eigen matrix)
  template <typename T>
  VertexPositionGeometry(ManifoldVolumeMesh& mesh_, const Eigen::MatrixBase<T>& vertexPositions);

  ~VertexPositionGeometry();

  // == Members
  ManifoldVolumeMesh& mesh;

  // == Utility methods

  // Recompute all require'd quantities from input data. Call this after e.g. repositioning a vertex or mutating the
  // mesh
  void refreshQuantities();

  // Clear out any cached quantities which were previously computed but are not currently required.
  void purgeQuantities();

  // Construct a new geometry which is exactly the same as this one, on the same mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  std::unique_ptr<VertexPositionGeometry> copy();

  // Construct a new geometry which is exactly the same as this one, on another mesh.
  // This is a deep copy, no quantites are shared, etc. Require counts/computed quantities are not copied.
  // The meshes must be in correspondence (have the same connectivity).
  std::unique_ptr<VertexPositionGeometry> reinterpretTo(ManifoldVolumeMesh& targetMesh);

  // Hide copy and move constructors; users are more likely to use them accidentally than intentionally.
  // See the explicit copy() function in derived classes.
  VertexPositionGeometry(const VertexPositionGeometry& other) = delete;
  VertexPositionGeometry& operator=(const VertexPositionGeometry& other) = delete;
  VertexPositionGeometry(VertexPositionGeometry&& other) = delete;
  VertexPositionGeometry& operator=(VertexPositionGeometry&& other) = delete;

  // === Quantities

  // == Indices
  // Note: These don't depend on any geometric information, and are no different than the getVertexIndices() offered by
  // the mesh class. However, its useful to offer them here so they can be used with the caching system.

  // Vertex indices
  VertexData<size_t> vertexIndices;
  void requireVertexIndices();
  void unrequireVertexIndices();

  EdgeData<size_t> edgeIndices;
  void requireEdgeIndices();
  void unrequireEdgeIndices();

  FaceData<size_t> faceIndices;
  void requireFaceIndices();
  void unrequireFaceIndices();

  CellData<size_t> cellIndices;
  void requireCellIndices();
  void unrequireCellIndices();

  VertexCornerData<size_t> vertexCornerIndices;
  void requireVertexCornerIndices();
  void unrequireVertexCornerIndices();

  EdgeCornerData<size_t> edgeCornerIndices;
  void requireEdgeCornerIndices();
  void unrequireEdgeCornerIndices();

  FaceCornerData<size_t> faceCornerIndices;
  void requireFaceCornerIndices();
  void unrequireFaceCornerIndices();

  // == Geometry
  VertexData<Vector3> vertexPositions;

  FaceData<Vector3> faceAreaNormals;
  void requireFaceAreaNormals();
  void unrequireFaceAreaNormals();

  FaceData<Vector3> faceNormals;
  void requireFaceNormals();
  void unrequireFaceNormals();

  EdgeData<Vector3> edgeVectors;
  void requireEdgeVectors();
  void unrequireEdgeVectors();

  EdgeData<double> edgeLengths;
  void requireEdgeLengths();
  void unrequireEdgeLengths();

  FaceData<double> faceAreas;
  void requireFaceAreas();
  void unrequireFaceAreas();

  CellData<double> cellVolumes;
  void requireCellVolumes();
  void unrequireCellVolumes();

  VertexData<double> vertexDualVolumes;
  void requireVertexDualVolumes();
  void unrequireVertexDualVolumes();

  FaceData<double> faceDualEdgeLengths;
  void requireFaceDualEdgeLengths();
  void unrequireFaceDualEdgeLengths();

  FaceData<double> faceHodge2;
  void requireFaceHodge2();
  void unrequireFaceHodge2();

  FaceCornerData<double> faceCornerAngles;
  void requireFaceCornerAngles();
  void unrequireFaceCornerAngles();

  FaceCornerData<double> faceCornerAngleCotans;
  void requireFaceCornerAngleCotans();
  void unrequireFaceCornerAngleCotans();

  EdgeCornerData<double> dihedralAngles;
  void requireDihedralAngles();
  void unrequireDihedralAngles();

  EdgeCornerData<double> dihedralAngleCotans;
  void requireDihedralAngleCotans();
  void unrequireDihedralAngleCotans();

  FaceData<Vector3> faceBarycenters;
  void requireFaceBarycenters();
  void unrequireFaceBarycenters();

  FaceData<Vector3> faceCircumcenters;
  void requireFaceCircumcenters();
  void unrequireFaceCircumcenters();

  CellData<Vector3> cellBarycenters;
  void requireCellBarycenters();
  void unrequireCellBarycenters();

  CellData<Vector3> cellCircumcenters;
  void requireCellCircumcenters();
  void unrequireCellCircumcenters();

  SparseMatrix<double> cotanLaplacian;
  void requireCotanLaplacian();
  void unrequireCotanLaplacian();

  SparseMatrix<double> vertexLumpedMassMatrix;
  void requireVertexLumpedMassMatrix();
  void unrequireVertexLumpedMassMatrix();

  // === Immediate computations
  Vector3 gradient(const VertexData<double>& u, Cell c);
  CellData<Vector3> gradient(const VertexData<double>& u);

  Vector3 faceGradient(const VertexData<double>& u, Face f);
  FaceData<Vector3> faceGradient(const VertexData<double>& u);

protected:
  // All of the quantities available (subclasses will also add quantities to this list)
  // Note that this is a vector of non-owning pointers; the quantities are generally value members in the class, so
  // there is no need to delete these.
  std::vector<DependentQuantity*> quantities;

  // === Implementation details for quantities

  // == Indices

  DependentQuantityD<VertexData<size_t>> vertexIndicesQ;
  void computeVertexIndices();

  DependentQuantityD<EdgeData<size_t>> edgeIndicesQ;
  void computeEdgeIndices();

  DependentQuantityD<FaceData<size_t>> faceIndicesQ;
  void computeFaceIndices();

  DependentQuantityD<CellData<size_t>> cellIndicesQ;
  void computeCellIndices();

  DependentQuantityD<VertexCornerData<size_t>> vertexCornerIndicesQ;
  void computeVertexCornerIndices();

  DependentQuantityD<EdgeCornerData<size_t>> edgeCornerIndicesQ;
  void computeEdgeCornerIndices();

  DependentQuantityD<FaceCornerData<size_t>> faceCornerIndicesQ;
  void computeFaceCornerIndices();

  // == Geometry
  DependentQuantityD<FaceData<Vector3>> faceAreaNormalsQ;
  void computeFaceAreaNormals();

  DependentQuantityD<FaceData<Vector3>> faceNormalsQ;
  void computeFaceNormals();

  DependentQuantityD<EdgeData<Vector3>> edgeVectorsQ;
  void computeEdgeVectors();

  DependentQuantityD<EdgeData<double>> edgeLengthsQ;
  void computeEdgeLengths();

  DependentQuantityD<FaceData<double>> faceAreasQ;
  void computeFaceAreas();

  DependentQuantityD<CellData<double>> cellVolumesQ;
  void computeCellVolumes();

  DependentQuantityD<VertexData<double>> vertexDualVolumesQ;
  void computeVertexDualVolumes();

  DependentQuantityD<FaceData<double>> faceDualEdgeLengthsQ;
  void computeFaceDualEdgeLengths();

  DependentQuantityD<FaceData<double>> faceHodge2Q;
  void computeFaceHodge2();

  DependentQuantityD<FaceCornerData<double>> faceCornerAnglesQ;
  void computeFaceCornerAngles();

  DependentQuantityD<FaceCornerData<double>> faceCornerAngleCotansQ;
  void computeFaceCornerAngleCotans();

  DependentQuantityD<EdgeCornerData<double>> dihedralAnglesQ;
  void computeDihedralAngles();

  DependentQuantityD<EdgeCornerData<double>> dihedralAngleCotansQ;
  void computeDihedralAngleCotans();

  DependentQuantityD<FaceData<Vector3>> faceBarycentersQ;
  void computeFaceBarycenters();

  DependentQuantityD<FaceData<Vector3>> faceCircumcentersQ;
  void computeFaceCircumcenters();

  DependentQuantityD<CellData<Vector3>> cellBarycentersQ;
  void computeCellBarycenters();

  DependentQuantityD<CellData<Vector3>> cellCircumcentersQ;
  void computeCellCircumcenters();

  DependentQuantityD<SparseMatrix<double>> cotanLaplacianQ;
  void computeCotanLaplacian();

  DependentQuantityD<SparseMatrix<double>> vertexLumpedMassMatrixQ;
  void computeVertexLumpedMassMatrix();
};

} // namespace volume
} // namespace geometrycentral
