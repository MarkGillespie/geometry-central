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
template <size_t k>
std::vector<Dart<D>> CombinatorialMap<D>::adjacentDarts(Cell<k, D> cell) const {
  static_assert(k <= D, "input cell dimension k must be less than or equal to complex dimension D");
  static_assert(k != 0, "vertex adjacent cells not implemented yet"); // TODO: implement this

  // to find all adjacent darts, we express the input cell as an orbit of dart maps,
  std::vector<Dart<D>> neighbors;

  std::set<Dart<D>> seenDarts;
  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(cell.dart());
  seenDarts.insert(cell.dart());

  while (!dartsToVisit.empty()) {
    Dart<D> curr = dartsToVisit.back();
    dartsToVisit.pop_back();
    neighbors.push_back(cell.dart());

    // You need to go in descending order to orient tets properly
    for (size_t iD = D; iD > 0; --iD) {
      if (iD != k) {
        Dart<D> next = curr.partner(iD - 1);
        if (seenDarts.find(next) == seenDarts.end()) {
          dartsToVisit.push_back(next);
          seenDarts.insert(next);
        }
      }
    }
  }

  return neighbors;
}

template <size_t D>
template <size_t k1, size_t k2>
std::vector<Cell<k2, D>> CombinatorialMap<D>::adjacentCells(Cell<k1, D> cell) const {
  static_assert(k1 <= D, "input cell dimension k1 must be less than or equal to complex dimension D");
  static_assert(k2 <= D, "output cell dimension k2 must be less than or equal to complex dimension D");
  static_assert(k1 != 0, "vertex adjacent cells not implemented yet"); // TODO: implement this

  // to find all adjacent k2-cells, we express the input cell as an orbit of dart maps,
  // and call d.cell<k2>() for each of these darts
  std::vector<Cell<k2, D>> neighbors;

  std::set<Dart<D>> seenDarts;
  std::set<Cell<k2, D>> seenCells; // TODO: profile against comparing with neighbors list?
  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(cell.dart());
  seenDarts.insert(cell.dart());

  while (!dartsToVisit.empty()) {
    Dart<D> curr = dartsToVisit.back();
    dartsToVisit.pop_back();

    Cell<k2, D> currCell = curr.template cell<k2>();
    if (seenCells.find(currCell) == seenCells.end()) {
      neighbors.push_back(currCell);
      seenCells.insert(currCell);
    }

    // You need to go in descending order to orient tets properly
    for (size_t iD = D; iD > 0; --iD) {
      if (iD != k1) {
        Dart<D> neighbor = curr.partner(iD - 1);
        if (seenDarts.find(neighbor) == seenDarts.end()) { // if we haven't seen neighbor yet
          dartsToVisit.push_back(neighbor);
          seenDarts.insert(neighbor);
        }
      }
    }
  }

  return neighbors;
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
inline Vertex<D> CombinatorialMap<D>::vertex(size_t index) {
  return cell<0>(index);
}

template <size_t D>
inline Edge<D> CombinatorialMap<D>::edge(size_t index) {
  return cell<1>(index);
}

template <size_t D>
inline Face<D> CombinatorialMap<D>::face(size_t index) {
  return cell<2>(index);
}

template <size_t D>
template <size_t k>
inline Cell<k, D> CombinatorialMap<D>::cell(size_t index) {
  return Cell<k, D>(this, index);
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
VertexData<D, size_t> CombinatorialMap<D>::getVertexIndices() {
  return getCellIndices<0>();
}

template <size_t D>
EdgeData<D, size_t> CombinatorialMap<D>::getEdgeIndices() {
  return getCellIndices<1>();
}

template <size_t D>
FaceData<D, size_t> CombinatorialMap<D>::getFaceIndices() {
  return getCellIndices<2>();
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
template <size_t k>
std::vector<std::vector<size_t>> CombinatorialMap<D>::getCellVertexList() {
  VertexData<D, size_t> vIdx = getVertexIndices();
  std::vector<std::vector<size_t>> cellVertexList;
  for (Cell<k, D> cell : cells<k>()) {
    cellVertexList.push_back({});
    for (Vertex<D> v : cell.adjacentVertices()) cellVertexList.back().push_back(vIdx[v]);
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
  const bool DEBUG_PRINT = false;

  // TODO : we should be able to figure out the orientations as well

  // use union-find to identify k-cells as orbits generated by compositions of dart maps
  // 0-cells are generated by <map[0].map[1], map[0].map[2], ..., map[D-2].map[D-1]>
  // 1-cells are generated by <map[1], ..., map[D-1]>
  // 2-cells are generated by <map[0], map[2], ..., map[D-1]>
  // see e.g. https://doc.cgal.org/latest/Combinatorial_map/index.html#title3

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
        for (size_t jMap = 0; jMap < iMap; jMap++) {
          size_t jDart = dartMap[iMap][iDart];
          if (jDart == INVALID_IND) continue;
          size_t kDart = dartMap[jMap][jDart];
          if (kDart == INVALID_IND) continue;
          unite(iDart, kDart);
        }
      }
    } else {
      for (size_t iMap = 0; iMap < D; iMap++) {
        if (iMap + 1 == k) continue;

        size_t jDart = dartMap[iMap][iDart];
        if (jDart == INVALID_IND) continue;
        unite(iDart, jDart);

        if (DEBUG_PRINT) {
          std::cout << "uniting dart " << iDart << " with dart " << jDart << " via map " << iMap << std::endl;
        }
      }
    }
  }

  if (DEBUG_PRINT) {
    std::cout << std::endl << "Final " << k << "-cell orbits: " << std::endl;
    for (size_t iRoot = 0; iRoot < nDarts(); iRoot++) {
      if (findRoot(iRoot) != iRoot) continue;
      std::cout << "  root " << iRoot << std::endl;
      for (size_t iDart = 0; iDart < nDarts(); iDart++) {
        if (findRoot(iDart) == iRoot) {
          std::cout << "     dart " << iDart << std::endl;
        }
      }
      std::cout << std::endl;
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
  const bool DEBUG_PRINT = false;
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
  auto flip = [&](std::array<size_t, 3> key) -> std::array<size_t, 3> { return {key[1], key[0], key[2]}; };

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
    return INVALID_IND;
  };

  // === Walk the tets, creating darts. Hook up dartMap[0] and dartMap[1] pointers (halfedges on tet surfaces), but
  // don't hook up dartMap[2] yet (gluing tets together).

  // The oriented faces of tet {0, 1, 2, 3} are given by {{0, 1, 2}, {0, 2, 3}, {1, 3, 2}, {0, 3, 1}}
  // We index the tet's halfedges as 0 1 2, 3 4 5, 6 7 8, 9 10 11
  // The next array is 1 2 0, 4 5 3, 7 8 6, 10 11 9
  // The twin array is 11 8 3, 2 7 9, 10 4 1, 5 6 0
  const std::array<std::array<size_t, 3>, 4> tetFaceIndices{
      std::array<size_t, 3>{0, 1, 2}, std::array<size_t, 3>{0, 2, 3}, std::array<size_t, 3>{1, 3, 2},
      std::array<size_t, 3>{0, 3, 1}};
  const std::array<size_t, 12> next{1, 2, 0, 4, 5, 3, 7, 8, 6, 10, 11, 9};
  const std::array<size_t, 12> twin{11, 8, 3, 2, 7, 9, 10, 4, 1, 5, 6, 0};

  for (size_t iTet = 0; iTet < tets.size(); iTet++) {
    const std::vector<size_t>& tet = tets[iTet];
    size_t iCell3 = getNewCell<3>().getIndex();

    const std::array<std::array<size_t, 3>, 4> tetFaces{
        std::array<size_t, 3>{tet[0], tet[1], tet[2]}, std::array<size_t, 3>{tet[0], tet[2], tet[3]},
        std::array<size_t, 3>{tet[1], tet[3], tet[2]}, std::array<size_t, 3>{tet[0], tet[3], tet[1]}};

    std::array<size_t, 12> newDartIndices;
    for (size_t iDart = 0; iDart < 12; ++iDart) {
      size_t newDart = getNewDart().getIndex();
      newDartIndices[iDart] = newDart;

      size_t iV = tetFaces[iDart / 3][iDart % 3];
      dCellArr[0][newDart] = iV;
      cDartArr[0][iV] = newDart;
    }

    for (size_t iDart = 0; iDart < 12; ++iDart) {
      dartMap[0][newDartIndices[iDart]] = newDartIndices[next[iDart]];
      dartMap[1][newDartIndices[iDart]] = newDartIndices[twin[iDart]];
      dCellArr[3][newDartIndices[iDart]] = iCell3;

      // initialize dartMap[2] to INVALID_IND (nothing glued together)
      dartMap[2][newDartIndices[iDart]] = INVALID_IND;
    }
    cDartArr[3][iCell3] = newDartIndices[0];

    if (DEBUG_PRINT) {
      for (size_t iDart = 0; iDart < 12; ++iDart) {
        std::cout << "Dart " << newDartIndices[iDart] << " : " << dCellArr[0][newDartIndices[iDart]] << "->"
                  << dCellArr[0][dartMap[0][newDartIndices[iDart]]] << std::endl;
      }
    }

    // glue together opposite faces
    for (size_t iF = 0; iF < 4; ++iF) {
      const std::array<size_t, 3> face = tetFaces[iF];
      size_t twinFaceDart = createdDartLookup(flip(face));
      if (twinFaceDart == INVALID_IND) {
        // if the opposite face has not been created, leave the dartMap[2] pointers empty
        createdDarts[face] = newDartIndices[3 * iF];
      } else {
        if (DEBUG_PRINT) {
          size_t myDart = newDartIndices[3 * iF + 0];
          size_t oppDart = twinFaceDart;
          std::cout << " ----- gluing " << dCellArr[0][myDart] << "->" << dCellArr[0][dartMap[0][myDart]] << "[dart "
                    << myDart << "] to " << dCellArr[0][oppDart] << "-> " << dCellArr[0][dartMap[0][oppDart]]
                    << "[dart " << oppDart << "]" << std::endl;
          myDart = newDartIndices[3 * iF + 2];
          oppDart = dartMap[0][twinFaceDart];
          std::cout << " ----- gluing " << dCellArr[0][myDart] << "->" << dCellArr[0][dartMap[0][myDart]] << "[dart "
                    << myDart << "] to " << dCellArr[0][oppDart] << "-> " << dCellArr[0][dartMap[0][oppDart]]
                    << "[dart " << oppDart << "]" << std::endl;
          myDart = newDartIndices[3 * iF + 1];
          oppDart = dartMap[0][dartMap[0][twinFaceDart]];
          std::cout << " ----- gluing " << dCellArr[0][myDart] << "->" << dCellArr[0][dartMap[0][myDart]] << "[dart "
                    << myDart << "] to " << dCellArr[0][oppDart] << "-> " << dCellArr[0][dartMap[0][oppDart]]
                    << "[dart " << oppDart << "]" << std::endl;
        }

        // if the opposite face has already created, hook up the appropriate pointers
        dartMap[2][newDartIndices[3 * iF + 0]] = twinFaceDart;
        dartMap[2][newDartIndices[3 * iF + 2]] = dartMap[0][twinFaceDart];
        dartMap[2][newDartIndices[3 * iF + 1]] = dartMap[0][dartMap[0][twinFaceDart]];
      }
    }
  }


  // TODO: do something about boundary

  nCellsCapacityCount[0] = nCellsCount[0];
  nCellsFillCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nDartsFillCount = nDartsCount;

  indexCells<1>();
  indexCells<2>();
} // namespace combinatorial_map

} // namespace combinatorial_map
} // namespace geometrycentral
