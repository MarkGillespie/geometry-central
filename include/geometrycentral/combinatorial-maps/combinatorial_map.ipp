#pragma once

namespace geometrycentral {
namespace combinatorial_map {

template <size_t D>
CombinatorialMap<D>::CombinatorialMap() {}

// Methods for getting number of mesh elements
template <size_t D>
inline size_t CombinatorialMap<D>::nDarts() const {
  return nDartsCount;
}

template <size_t D>
template <size_t k>
inline size_t CombinatorialMap<D>::nCells() const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return nCellsCount[k];
}

template <size_t D>
inline size_t CombinatorialMap<D>::nVertices() const {
  return nCells<0>();
}

template <size_t D>
inline size_t CombinatorialMap<D>::nEdges() const {
  return nCells<1>();
}

template <size_t D>
inline size_t CombinatorialMap<D>::nFaces() const {
  return nCells<2>();
}

// Capacities
template <size_t D>
inline size_t CombinatorialMap<D>::nDartsCapacity() const {
  return nDartsCapacityCount;
}

template <size_t D>
template <size_t k>
inline size_t CombinatorialMap<D>::nCellsCapacity() const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return nCellsCapacityCount[k];
}

// Connectivity
template <size_t D>
inline size_t CombinatorialMap<D>::dartPartner(size_t iD, size_t dim) const {
  return dartMap[dim][iD];
}

// template <size_t D>
// inline size_t CombinatorialMap<D>::heNextIncomingNeighbor(size_t iD)  const {
//   return usesImplicitTwin() ? heTwinImplicit(heNextArr[iD]) : heVertInNextArr[iD];
// }

// template <size_t D>
// inline size_t CombinatorialMap<D>::heNextOutgoingNeighbor(size_t iD) const {
//   return usesImplicitTwin() ? heNextArr[heTwinImplicit(iD)] : heVertOutNextArr[iD];
// }

template <size_t D>
Dart<D> CombinatorialMap<D>::getNewDart() {

  // The boring case, when no resize is needed
  if (nDartsFillCount < nDartsCapacityCount) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newDartCapacity = std::max(nDartsCapacityCount * 2, (size_t)1);

    // Resize internal arrays
    for (size_t iD = 0; iD < D; ++iD) {
      dartMap[iD].resize(newDartCapacity);
    }
    for (size_t iD = 0; iD <= D; ++iD) {
      dCellArr[iD].resize(newDartCapacity);
    }

    nDartsCapacityCount = newDartCapacity;

    // Invoke relevant callback functions
    for (auto& f : dartExpandCallbackList) {
      f(newDartCapacity);
    }
  }

  nDartsFillCount++;
  nDartsCount++;

  modificationTick++;
  return Dart<D>(this, nDartsFillCount - 1);
}

template <size_t D>
template <size_t k>
Cell<k, D> CombinatorialMap<D>::getNewCell() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");

  // The boring case, when no resize is needed
  if (nCellsFillCount[k] < nCellsCapacityCount[k]) {
    // No work needed
  }
  // The intesting case, where vectors resize
  else {
    size_t newCellCapacity = std::max(nCellsCapacityCount[k] * 2, (size_t)1);

    // Resize internal arrays
    cDartArr[k].resize(newCellCapacity);

    nCellsCapacityCount[k] = newCellCapacity;

    // Invoke relevant callback functions
    for (auto& f : cellExpandCallbackList[k]) {
      f(newCellCapacity);
    }
  }

  nCellsFillCount[k]++;
  nCellsCount[k]++;

  modificationTick++;
  return Cell<k, D>(this, nCellsFillCount[k] - 1);
}

template <size_t D>
template <size_t k>
inline bool CombinatorialMap<D>::cellIsDead(size_t iC) const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return cDartArr[k][iC] == INVALID_IND;
}

template <size_t D>
inline bool CombinatorialMap<D>::dartIsDead(size_t iD) const {
  return dartMap[0][iD] == INVALID_IND;
}

// Methods for iterating over mesh elements w/ range-based for loops ===========

template <size_t D>
inline DartSet<D> CombinatorialMap<D>::darts() {
  return DartSet<D>(this, 0, nDartsFillCount);
}

template <size_t D>
inline VertexSet<D> CombinatorialMap<D>::vertices() {
  return cells<0>();
}

template <size_t D>
inline EdgeSet<D> CombinatorialMap<D>::edges() {
  return cells<1>();
}

template <size_t D>
inline FaceSet<D> CombinatorialMap<D>::faces() {
  return cells<2>();
}

template <size_t D>
template <size_t k>
inline CellSet<k, D> CombinatorialMap<D>::cells() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return CellSet<k, D>(this, 0, nCellsFillCount[k]);
}

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.

template <size_t D>
inline Dart<D> CombinatorialMap<D>::dart(size_t index) {
  return Dart<D>(this, index);
}

template <size_t D>
template <size_t k>
CellData<k, D, size_t> CombinatorialMap<D>::getCellIndices() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  CellData<k, D, size_t> indices(*this);
  size_t i = 0;
  for (Cell<k, D> c : cells<k>()) {
    indices[c] = i;
    i++;
  }
  return indices;
}

template <size_t D>
DartData<D, size_t> CombinatorialMap<D>::getDartIndices() {
  DartData<D, size_t> indices(*this);
  size_t i = 0;
  for (Dart<D> dart : darts()) {
    indices[dart] = i;
    i++;
  }
  return indices;
}

template <size_t D>
template <size_t E>
std::vector<std::vector<size_t>> CombinatorialMap<D>::getCellVertexList() {
  std::vector<std::vector<size_t>> cellVertexList;
  // VertexData<D, size_t> vIdx = getVertexIndices();

  DartData<D, char> visited(*this, false);
  DartData<D, char> onStack(*this, false);

  for (Dart<D> d : darts()) {
    if (visited[d]) continue;

    std::vector<size_t> cellVertices;
    std::deque<Dart<D>> dartsToVisit;
    dartsToVisit.push_back(d);
    onStack[d] = true;

    while (!dartsToVisit.empty()) {
      Dart<D> curr = dartsToVisit.back();
      dartsToVisit.pop_back();

      size_t currIdx = curr.vertex().getIndex();
      if (std::find(cellVertices.begin(), cellVertices.end(), currIdx) == cellVertices.end()) {
        cellVertices.push_back(currIdx);
      }

      visited[curr] = true;

      // You need to go in descending order to orient tets properly
      for (size_t iD = D; iD > 0; --iD) {
        if (iD != E) {
          Dart<D> neighbor = curr.partner(iD - 1);
          if (!visited[neighbor] && !onStack[neighbor]) {
            dartsToVisit.push_back(neighbor);
            onStack[neighbor] = true;
          }
        }
      }
    }
    cellVertexList.push_back(cellVertices);
  }
  return cellVertexList;
}

// Misc utility methods =====================================

template <size_t D>
inline bool CombinatorialMap<D>::isCompressed() const {
  return isCompressedFlag;
}

template <size_t D>
CombinatorialMap<D>::~CombinatorialMap() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}

template <size_t D>
std::unique_ptr<CombinatorialMap<D>> CombinatorialMap<D>::copy() const {
  return copyToCombinatorialMap();
}

template <size_t D>
std::unique_ptr<CombinatorialMap<D>> CombinatorialMap<D>::copyToCombinatorialMap() const {
  CombinatorialMap<D>* newMesh = new CombinatorialMap<D>();
  copyInternalFields(*newMesh);
  return std::unique_ptr<CombinatorialMap<D>>(newMesh);
}

template <size_t D>
void CombinatorialMap<D>::copyInternalFields(CombinatorialMap<D>& target) const {
  // == Copy _all_ the fields!

  // Raw data buffers (underlying std::vectors duplicate storage automatically)
  // TODO: does this still do a deep copy now that this is an std::array?
  target.dartMap = dartMap;
  target.dCellArr = dCellArr;
  target.cDartArr = cDartArr;

  // counts and flags
  target.nDartsCount = nDartsCount;
  target.nDartsCapacityCount = nDartsCapacityCount;
  target.nCellsCapacityCount = nCellsCapacityCount;
  target.nDartsFillCount = nDartsFillCount;
  target.nCellsFillCount = nCellsFillCount;

  target.isCompressedFlag = isCompressedFlag;

  // Note: _don't_ copy callbacks lists! New mesh has new callbacks
}

// index k-cells and fill cDartArr[k] and dCellArr[k] based off of dartMap
template <size_t D>
template <size_t k>
void CombinatorialMap<D>::indexCells() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");

  // use union-find to identify k-cells as subgroups generated by compositions of dart maps
  // 0-cells are generated by <map[0].map[1], map[0].map[2], ..., map[0].map[D-1]>
  // 1-cells are generated by <map[1], ..., map[D-1]>
  // 2-cells are generated by <map[0], map[2], ..., map[D-1]>

  std::vector<size_t> parent;
  parent.reserve(nDarts()); // initialize every dart as its own parent
  for (size_t i = 0; i < nDarts(); i++) parent.push_back(i);

  std::vector<size_t> rank(nDarts(), 0); // initialize every dart to rank 0

  auto findRoot = [&parent](size_t x) -> size_t {
    std::vector<size_t> visitedNodes;
    while (parent[x] != x) {
      visitedNodes.push_back(x);
      x = parent[x];
    }
    for (size_t n : visitedNodes) parent[n] = x;
    return x;
  };

  auto unite = [&parent, &rank, &findRoot](size_t x, size_t y) -> void {
    size_t rootX = findRoot(x), rootY = findRoot(y);
    if (rootX == rootY) return;

    // Union by rank
    if (rank[rootX] < rank[rootY]) {
      parent[rootX] = rootY;
    } else if (rank[rootX] > rank[rootY]) {
      parent[rootY] = rootX;
    } else {
      parent[rootY] = rootX;
      rank[rootX]++;
    }
  };

  // identify subgroups
  for (size_t iDart = 0; iDart < nDarts(); iDart++) {
    // special case for k=0
    if (k == 0) {
      for (size_t iMap = 1; iMap < D; iMap++) {
        size_t jDart = dartMap[0][dartMap[iMap][iDart]];
        if (jDart != INVALID_IND) unite(iDart, jDart);
      }
    } else {
      for (size_t iMap = 0; iMap < D; iMap++) {
        if (iMap == k) continue;
        size_t jDart = dartMap[iMap][iDart];
        if (jDart != INVALID_IND) unite(iDart, jDart);
      }
    }
  }

  // Clear any existing k-cells
  std::fill(dCellArr[k].begin(), dCellArr[k].end(), INVALID_IND);
  nCellsFillCount[k] = 0;
  nCellsCount[k] = 0;

  // Assign new k-cell indices
  for (size_t iDart = 0; iDart < nDarts(); iDart++) {
    size_t iRoot = findRoot(iDart);
    if (dCellArr[k][iRoot] == INVALID_IND) {
      dCellArr[k][iRoot] = getNewCell<k>().getIndex();
      cDartArr[k][dCellArr[k][iRoot]] = iRoot;
    }
    dCellArr[k][iDart] = dCellArr[k][iRoot];
  }
}


// Builds a halfedge mesh
template <>
CombinatorialMap<2>::CombinatorialMap(const std::vector<std::vector<size_t>>& polygons) {
  surface::ManifoldSurfaceMesh mesh(polygons);

  nCellsCount[0] = mesh.nVertices();
  nDartsCount = mesh.nHalfedges();
  nCellsCapacityCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nCellsFillCount[0] = nCellsCount[0];
  nDartsFillCount = nDartsCount;

  // TODO: copy over arrays?
  surface::VertexData<size_t> vIdx = mesh.getVertexIndices();
  surface::EdgeData<size_t> eIdx = mesh.getEdgeIndices();
  surface::FaceData<size_t> fIdx = mesh.getFaceIndices();
  surface::HalfedgeData<size_t> hIdx = mesh.getHalfedgeIndices();
  dartMap[0].reserve(mesh.nHalfedges());
  dartMap[1].reserve(mesh.nHalfedges());
  dCellArr[0].reserve(mesh.nHalfedges());
  dCellArr[1].reserve(mesh.nHalfedges());
  dCellArr[2].reserve(mesh.nHalfedges());
  cDartArr[0].reserve(mesh.nVertices());
  for (surface::Halfedge he : mesh.halfedges()) {
    dartMap[0].push_back(hIdx[he.next()]);
    dartMap[1].push_back(hIdx[he.twin()]);
    dCellArr[0].push_back(vIdx[he.vertex()]);
    dCellArr[1].push_back(eIdx[he.edge()]);
    dCellArr[2].push_back(fIdx[he.face()]);
  }
  for (surface::Vertex v : mesh.vertices()) {
    cDartArr[0].push_back(hIdx[v.halfedge()]);
  }
}

// Builds a tet mesh
template <>
CombinatorialMap<3>::CombinatorialMap(const std::vector<std::vector<size_t>>& tets) {
  nCellsCount[0] = 0;
  for (const std::vector<size_t>& tet : tets) {
    GC_SAFETY_ASSERT(tet.size() == 4, "CombinatorialMap<3> can only construct tet meshes from a list of cell vertices");
    for (size_t i : tet) {
      nCellsCount[0] = std::max(nCellsCount[0], i);
    }
  }
  nCellsCount[0]++; // 0-based means count is max + 1

  cDartArr[0] = std::vector<size_t>(nCellsCount[0], INVALID_IND);

  std::map<std::array<size_t, 3>, size_t> createdDarts;

  auto shift = [&](std::array<size_t, 3> key) -> std::array<size_t, 3> { return {key[1], key[2], key[0]}; };

  auto createdDartLookup = [&](std::array<size_t, 3> key) -> size_t {
    auto keyIter = createdDarts.find(key);
    if (keyIter != createdDarts.end()) {
      return keyIter->second;
    }
    keyIter = createdDarts.find(shift(key));
    if (keyIter != createdDarts.end()) {
      return dartMap[0][dartMap[0][keyIter->second]];
    }
    keyIter = createdDarts.find(shift(shift(key)));
    if (keyIter != createdDarts.end()) {
      return dartMap[0][keyIter->second];
    }
    createdDarts[key] = INVALID_IND;
    return INVALID_IND;
  };

  // === Walk the tets, creating darts. Hook up dartMap[0] and dartMap[1] pointers (halfedges on tet surfaces), but
  // don't hook up dartMap[2] yet (gluing tets together).

  for (size_t iTet = 0; iTet < tets.size(); iTet++) {
    const std::vector<size_t>& tet = tets[iTet];
    size_t iCell3 = getNewCell<3>().getIndex();

    // The oriented faces of tet {0, 1, 2, 3} are given by {{0, 1, 2}, {0, 2, 3}, {1, 3, 2}, {0, 3, 1}}
    // We index the tet's halfedges as 0 1 2, 3 4 5, 6 7 8, 9 10 11
    // The next array is 1 2 0, 4 5 3, 7 8 6, 10 11 9
    // The twin array is 11 8 3, 2 7 9, 10 4 1, 5 6 0
    const std::array<std::array<size_t, 3>, 4> tetFaceIndices{
        std::array<size_t, 3>{0, 1, 2}, std::array<size_t, 3>{0, 2, 3}, std::array<size_t, 3>{1, 3, 2},
        std::array<size_t, 3>{0, 3, 1}};
    const std::array<std::array<size_t, 3>, 4> tetFaces{
        std::array<size_t, 3>{tet[0], tet[1], tet[2]}, std::array<size_t, 3>{tet[0], tet[2], tet[3]},
        std::array<size_t, 3>{tet[1], tet[3], tet[2]}, std::array<size_t, 3>{tet[0], tet[3], tet[1]}};
    const std::array<size_t, 12> next{1, 2, 0, 4, 5, 3, 7, 8, 6, 10, 11, 9};
    const std::array<size_t, 12> twin{11, 8, 3, 2, 7, 9, 10, 4, 1, 5, 6, 0};

    std::array<size_t, 12> createdDarts;
    for (size_t iDart = 0; iDart < 12; ++iDart) {
      size_t newDart = getNewDart().getIndex();
      createdDarts[iDart] = newDart;

      size_t iV = tetFaces[iDart / 3][iDart % 3];
      dCellArr[0][newDart] = iV;
      cDartArr[0][iV] = newDart;
    }

    for (size_t iDart = 0; iDart < 12; ++iDart) {
      dartMap[0][createdDarts[iDart]] = createdDarts[next[iDart]];
      dartMap[1][createdDarts[iDart]] = createdDarts[twin[iDart]];
      dCellArr[3][createdDarts[iDart]] = iCell3;
    }
    cDartArr[3][iCell3] = createdDarts[0];

    for (size_t iF = 0; iF < 4; ++iF) {
      const std::array<size_t, 3> face = tetFaces[iF];
      size_t twinFaceDart = createdDartLookup(face);
      if (twinFaceDart == INVALID_IND) {
        // if the opposite face has not been created, set the appropriate pointers to empty
        for (size_t iDart : tetFaceIndices[iF]) {
          dartMap[2][createdDarts[iDart]] = INVALID_IND;
        }
      } else {
        // if the opposite face has already created, hook up the appropriate pointers
        dartMap[2][createdDarts[tetFaceIndices[iF][0]]] = twinFaceDart;
        dartMap[2][createdDarts[tetFaceIndices[iF][2]]] = dartMap[0][twinFaceDart];
        dartMap[2][createdDarts[tetFaceIndices[iF][1]]] = dartMap[0][dartMap[0][twinFaceDart]];
      }
    }
  }

  // TODO: do something about boundary

  nCellsCapacityCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nCellsFillCount[0] = nCellsCount[0];
  nDartsFillCount = nDartsCount;

  indexCells<1>();
  indexCells<2>();
} // namespace combinatorial_map

} // namespace combinatorial_map
} // namespace geometrycentral
