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
std::vector<Dart<D>> CombinatorialMap<D>::adjacentDarts(Cell<k, D> cell) {
  static_assert(k <= D, "input cell dimension k must be less than or equal to complex dimension D");

  // to find all adjacent darts, we express the input cell as an orbit of dart maps,
  std::vector<Dart<D>> neighbors;
  neighbors.reserve(16); // reserve some amount of space

  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(cell.dart());
  neighbors.push_back(cell.dart());

  while (!dartsToVisit.empty()) {
    // for some reason, iterating in DFS order is important for orienting cells
    Dart<D> curr = dartsToVisit.back();
    dartsToVisit.pop_back();

    for (std::pair<Dart<D>, bool> n : orbitNeighbors<k>(curr)) {
      if (std::find(neighbors.begin(), neighbors.end(), n.first) == neighbors.end()) {
        neighbors.push_back(n.first);
        dartsToVisit.push_back(n.first);
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

  // to find all adjacent k2-cells, we express the input cell as an orbit of dart maps,
  // and call d.cell<k2>() for each of these darts
  std::vector<Cell<k2, D>> neighbors;
  std::vector<Dart<D>> seenDarts;
  neighbors.reserve(8);  // reserve some amount of space
  seenDarts.reserve(16); // reserve some amount of space

  std::deque<std::pair<Dart<D>, bool>> dartsToVisit;
  dartsToVisit.push_back(std::make_pair(cell.dart(), true));
  seenDarts.push_back(cell.dart());

  while (!dartsToVisit.empty()) {
    // for some reason, iterating in DFS order is important for orienting cells
    Dart<D> currDart;
    bool currOrientation;
    std::tie(currDart, currOrientation) = dartsToVisit.back();
    dartsToVisit.pop_back();

    Cell<k2, D> currCell = currDart.template cell<k2>();
    currCell.setOrientation(currCell.orientation() == currOrientation);
    if (std::find(neighbors.begin(), neighbors.end(), currCell) == neighbors.end()) {
      neighbors.push_back(currCell);
    }

    for (std::pair<Dart<D>, bool> n : orbitNeighbors<k1>(currDart)) {
      if (std::find(seenDarts.begin(), seenDarts.end(), n.first) == seenDarts.end()) {
        dartsToVisit.push_back(std::make_pair(n.first, currOrientation == n.second));
        seenDarts.push_back(n.first);
      }
    }
  }

  return neighbors;
}

template <size_t D>
template <size_t k>
std::vector<Vertex<D>> CombinatorialMap<D>::adjacentVertices(Cell<k, D> cell) const {
  return adjacentCells<k, 0>(cell);
}

template <size_t D>
template <size_t k>
std::vector<Edge<D>> CombinatorialMap<D>::adjacentEdges(Cell<k, D> cell) const {
  return adjacentCells<k, 1>(cell);
}

template <size_t D>
template <size_t k>
std::vector<Face<D>> CombinatorialMap<D>::adjacentFaces(Cell<k, D> cell) const {
  return adjacentCells<k, 2>(cell);
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
inline size_t CombinatorialMap<D>::nVerticesCapacity() const {
  return nCellsCapacity<0>();
}
template <size_t D>
inline size_t CombinatorialMap<D>::nEdgesCapacity() const {
  return nCellsCapacity<1>();
}
template <size_t D>
inline size_t CombinatorialMap<D>::nFacesCapacity() const {
  return nCellsCapacity<2>();
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
  assert(dim <= D);
  assert(iD < nDarts());
  return dartMap[dim][iD];
}

template <size_t D>
size_t CombinatorialMap<D>::dartIndexSize() const {
  return nDartsFillCount;
}
template <size_t D>
size_t CombinatorialMap<D>::vertexIndexSize() const {
  return cellIndexSize<0>();
}
template <size_t D>
size_t CombinatorialMap<D>::edgeIndexSize() const {
  return cellIndexSize<1>();
}
template <size_t D>
size_t CombinatorialMap<D>::faceIndexSize() const {
  return cellIndexSize<2>();
}
template <size_t D>
template <size_t k>
size_t CombinatorialMap<D>::cellIndexSize() const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return nCellsFillCount[k];
}

template <size_t D>
template <size_t k>
OrbitNeighborhood<D, k> CombinatorialMap<D>::orbitNeighbors(Dart<D> d) const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return OrbitNeighborhood<D, k>(d);
}

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
      dCellSgn[iD].resize(newDartCapacity);
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
inline bool CombinatorialMap<D>::cellIsDead(size_t k, size_t iC) const {
  return k > D || cDartArr[k][iC] == INVALID_IND;
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
SparseMatrix<int> CombinatorialMap<D>::getBoundaryMatrix() {
  static_assert(k > 0, "Boundary_0 matrix not defined");
  static_assert(k <= D, "Boundary_k matrix not defined for k > complex dimension D");

  std::vector<Eigen::Triplet<int>> triplets;

  CellData<k, D, size_t> kIndices = getCellIndices<k>();
  CellData<k - 1, D, size_t> bdyIndices = getCellIndices<k - 1>();

  for (Cell<k, D> cell : cells<k>()) {
    for (Cell<k - 1, D> bdyCell : cell.template adjacentCells<k - 1>()) {
      triplets.emplace_back(bdyIndices[bdyCell], kIndices[cell], bdyCell.sign());
    }
  }

  SparseMatrix<int> bdy(nCells<k - 1>(), nCells<k>());
  bdy.setFromTriplets(triplets.begin(), triplets.end());
  return bdy;
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
  target.dCellSgn = dCellSgn;
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

  // use union-find to identify k-cells as orbits generated by compositions of dart maps
  // 0-cells are generated by <map[0].map[1], map[0].map[2], ..., map[D-2].map[D-1]>
  // 1-cells are generated by <map[1], ..., map[D-1]>
  // 2-cells are generated by <map[0], map[2], ..., map[D-1]>
  // see e.g. https://doc.cgal.org/latest/Combinatorial_map/index.html#title3

  std::vector<size_t> parent;
  parent.reserve(nDarts()); // initialize every dart as its own parent
  for (size_t i = 0; i < nDarts(); i++) parent.push_back(i);

  std::vector<size_t> rank(nDarts(), 0); // initialize every dart to rank 0
  std::vector<bool> sharesParentSign(nDarts(), true);

  // find root, and update all nodes in path to point to root, updating their `sharesParentSign` fields as necessary
  auto findRoot = [&parent, &sharesParentSign](size_t x) -> size_t {
    // early return if x or its parent is the root, in which case we don't have to update anything
    if (parent[x] == parent[parent[x]]) return parent[x];

    // otherwise, update all nodes between x and the root to be direct children of the root
    std::vector<size_t> visitedNodes;
    while (parent[x] != x) {
      visitedNodes.push_back(x);
      x = parent[x];
    }

    // iterate through in "last in first out" order, updating sharesParentSign
    // note that for booleans, a==b is the same as a xor b, i.e. sign multiplication
    bool runningSign = true;
    for (int iN = visitedNodes.size() - 1; iN >= 0; iN--) {
      parent[visitedNodes[iN]] = x;
      sharesParentSign[visitedNodes[iN]] = (sharesParentSign[visitedNodes[iN]] == runningSign);
      runningSign = sharesParentSign[visitedNodes[iN]];
    }

    return x;
  };

  // join together x and y, updating their `sharesParentSign` fields as necessary
  auto unite = [&parent, &rank, &sharesParentSign, &findRoot](size_t x, size_t y, bool samesign) -> void {
    size_t rootX = findRoot(x), rootY = findRoot(y);
    if (rootX == rootY) return;

    // Union by rank
    // note that for booleans, a==b is the same as a xor b, i.e. sign multiplication
    if (rank[rootX] < rank[rootY]) {
      parent[rootX] = rootY;
      sharesParentSign[rootX] = ((sharesParentSign[x] == sharesParentSign[y]) == samesign);
    } else if (rank[rootX] > rank[rootY]) {
      parent[rootY] = rootX;
      sharesParentSign[rootY] = ((sharesParentSign[x] == sharesParentSign[y]) == samesign);
    } else {
      parent[rootY] = rootX;
      sharesParentSign[rootY] = ((sharesParentSign[x] == sharesParentSign[y]) == samesign);
      rank[rootX]++;
    }
  };

  for (Dart<3> d : darts()) {
    for (std::pair<Dart<3>, bool> n : orbitNeighbors<k>(d)) {
      unite(d.getIndex(), n.first.getIndex(), n.second);
      if (DEBUG_PRINT) {
        std::cout << "uniting dart " << d.getIndex() << " with dart " << n.first.getIndex()
                  << " | orientationPreserving: " << (n.second ? "true" : "false") << std::endl;
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
          std::cout << "     dart " << iDart << " : " << dCellArr[0][iDart] << "->" << dCellArr[0][dartMap[0][iDart]]
                    << "\tpositive orientation: " << (sharesParentSign[iDart] ? "true" : "false") << std::endl;
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
      dCellSgn[k][iRoot] = true;
      cDartArr[k][dCellArr[k][iRoot]] = iRoot;
    }
    dCellArr[k][iDart] = dCellArr[k][iRoot];
    dCellSgn[k][iDart] = sharesParentSign[iDart];
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
  dCellSgn[0].reserve(mesh.nHalfedges());
  dCellSgn[1].reserve(mesh.nHalfedges());
  dCellSgn[2].reserve(mesh.nHalfedges());
  cDartArr[0].reserve(mesh.nVertices());
  for (surface::Halfedge he : mesh.halfedges()) {
    dartMap[0].push_back(hIdx[he.next()]);
    dartMap[1].push_back(hIdx[he.twin()]);
    dCellArr[0].push_back(vIdx[he.vertex()]);
    dCellArr[1].push_back(eIdx[he.edge()]);
    dCellArr[2].push_back(fIdx[he.face()]);
    dCellSgn[0].push_back(true);
    dCellSgn[1].push_back(he.orientation());
    dCellSgn[2].push_back(true);
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

  auto attachDartMap2 = [&](size_t iDart, size_t jDart) -> void {
    dartMap[2][iDart] = jDart;
    dartMap[2][jDart] = iDart;
  };

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
      dCellSgn[0][newDart] = true;
      cDartArr[0][iV] = newDart;
    }

    for (size_t iDart = 0; iDart < 12; ++iDart) {
      dartMap[0][newDartIndices[iDart]] = newDartIndices[next[iDart]];
      dartMap[1][newDartIndices[iDart]] = newDartIndices[twin[iDart]];
      dCellArr[3][newDartIndices[iDart]] = iCell3;
      dCellSgn[3][newDartIndices[iDart]] = true;
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
        // if the opposite face has not been created, set the dartMap[2] pointers to INVALID_IND
        for (size_t iD = 0; iD < 3; iD++) dartMap[2][newDartIndices[3 * iF + iD]] = INVALID_IND;
        createdDarts[face] = newDartIndices[3 * iF];
      } else {
        // if the opposite face has already created, hook up the appropriate pointers
        attachDartMap2(newDartIndices[3 * iF + 0], twinFaceDart);
        attachDartMap2(newDartIndices[3 * iF + 2], dartMap[0][twinFaceDart]);
        attachDartMap2(newDartIndices[3 * iF + 1], dartMap[0][dartMap[0][twinFaceDart]]);

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
      }
    }
  }


  nCellsCapacityCount[0] = nCellsCount[0];
  nCellsFillCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nDartsFillCount = nDartsCount;

  // construct 1-cells and 2-cells
  indexCells<1>();
  indexCells<2>();
}

template <size_t D>
void CombinatorialMap<D>::validateConnectivity() {
  // Sanity check sizes and counts
  if (nDartsCount > nDartsFillCount) throw std::logic_error("dart count > dart fill");
  if (nDartsFillCount > nDartsCapacityCount) throw std::logic_error("dart fill > dart capacity");

  for (size_t dim = 0; dim <= D; dim++) {
    if (nCellsCount[dim] > nCellsFillCount[dim])
      throw std::logic_error(std::to_string(dim) + "-cell count > " + std::to_string(dim) + "-cell fill");
    if (nCellsFillCount[dim] > nCellsCapacityCount[dim])
      throw std::logic_error(std::to_string(dim) + "-cell fill > " + std::to_string(dim) + "-cell capacity");
  }

  // Check for overflow / other unreasonable values
  if (nDartsCount > std::numeric_limits<uint64_t>::max() / 2) throw std::logic_error("dart count overflow");

  for (size_t dim = 0; dim <= D; dim++) {
    if (nCellsCount[dim] > std::numeric_limits<uint64_t>::max() / 2)
      throw std::logic_error(std::to_string(dim) + "-cell overflow");
  }

  // Helpers to check the validity of references
  auto validateDart = [&](size_t iDart, std::string msg, bool allowInvalidInd) {
    if ((iDart == INVALID_IND) && allowInvalidInd) return;
    if (iDart >= nDartsFillCount || dartIsDead(iDart))
      throw std::logic_error(msg + " | " + std::to_string(iDart) + " - bad dart reference");
  };
  auto validateCell = [&](size_t k, size_t iC, std::string msg) {
    if (iC >= nCellsFillCount[k] || cellIsDead(k, iC))
      throw std::logic_error(msg + " - bad " + std::to_string(k) + "-cell reference");
  };

  // == Darts
  // Check valid pointers
  // Note: we intentionally mostly avoid using iterators here, because they can be hard to debug when things are broken.
  for (size_t iDart = 0; iDart < nDartsFillCount; iDart++) {
    assert(!dartIsDead(iDart)); // no darts should be dead yet
    if (dartIsDead(iDart)) continue;
    // check partner, only allow INVALID_IND partner for dim>0 (i.e. not `next` map)
    for (size_t dim = 0; dim < D; dim++) {
      validateDart(dartMap[dim][iDart], "dart.partner(" + std::to_string(dim) + ")", dim > 0);
    }
    for (size_t k = 0; k <= D; k++) validateCell(k, dCellArr[k][iDart], "he.cell<" + std::to_string(k) + ">()");
  }
  for (size_t k = 0; k <= D; k++) {
    for (size_t iC = 0; iC < nCellsFillCount[k]; iC++) {
      if (cellIsDead(k, iC)) continue;
      validateDart(cDartArr[k][iC], "cell.dart()", false);
    }
  }

  // check adjacency sanity
  // TODO: loop over all dimensions?
  for (Vertex<D> v : vertices()) {
    for (Dart<D> d : v.adjacentDarts()) {
      if (v != d.vertex()) throw std::logic_error("vertex dart doesn't match dart.vertex");
    }
  }
  for (Edge<D> e : edges()) {
    for (Dart<D> d : e.adjacentDarts()) {
      if (e != d.edge()) throw std::logic_error("edge dart doesn't match dart.edge");
    }
  }
  for (Face<D> e : faces()) {
    for (Dart<D> d : e.adjacentDarts()) {
      if (e != d.face()) throw std::logic_error("face dart doesn't match dart.face");
    }
  }
  for (Cell<3, D> e : cells<3>()) {
    for (Dart<D> d : e.adjacentDarts()) {
      if (e != d.template cell<3>()) throw std::logic_error("3-cell dart doesn't match dart.cell<3>");
    }
  }
}

// ==========================================================
// ================= Special  Iterators   ===================
// ==========================================================

// iterate over k-cells, viewed as orbits of darts generated by compositions of dart maps
// 0-cells are generated by <map[0].map[1], map[0].map[2], ..., map[D-2].map[D-1]>
// 1-cells are generated by <map[1], ..., map[D-1]>
// 2-cells are generated by <map[0], map[2], ..., map[D-1]>
// ...

// when k=0, orbit over all map[i].map[j] for i < j
// when k>0, orbit over all map[i] for i != k-1
// more explicitly, these iterators are roughly equivalent to the following loops

// // ========= k = 0
// for (size_t iMap = 1; iMap < D; ++iMap) {
//   for (size_t jMap = 0; jMap < iMap; ++jMap) {
//     Dart<D> iPartner = d.partner(iMap);
//     if (iPartner == d) continue; // skip (INVALID_IND)
//     Dart<D> next = iPartner.partner(jMap);
//     if (next == iPartner) continue; // skip (INVALID_IND)
//     bool orientationPreserving = true;
//     neighbors.push_back(std::make_pair(next, orientationPreserving));
//   }
// }
//
// // ========= k > 0
// for (size_t iMap = 0; iMap < D; ++iMap) {
//   if (iMap + 1 != k) { // iMap != k-1
//     Dart<D> next = d.partner(iMap);
//     if (next == d) continue; // skip (INVALID_IND)
//     bool orientationPreserving = size_t(iMap + 1) < k;
//     neighbors.push_back(std::make_pair(next, orientationPreserving));
//   }
// }

template <size_t D, size_t k>
struct OrbitNeighborhoodIterator {
  static_assert(k > 0 && k <= D, "k must be in range (0, D]");

  Dart<D> startDart;
  size_t iMap;

  OrbitNeighborhoodIterator(Dart<D> d) : startDart(d), iMap(0) {
    while (!isValid() && !finished()) iMap++;
  }

  bool isValid() const {
    bool result = (iMap + 1 != k) && ((**this).first != startDart);
    return result;
  }

  bool finished() const { return iMap >= D; }

  const OrbitNeighborhoodIterator<D, k>& operator++() {
    do {
      iMap++;
    } while (!isValid() && !finished());
    return *this;
  }

  // any two finished iterators are equal, otherwise compare internals
  bool operator==(const OrbitNeighborhoodIterator<D, k>& other) const {
    return (finished() && other.finished()) || (startDart == other.startDart && iMap == other.iMap);
  }

  bool operator!=(const OrbitNeighborhoodIterator<D, k>& other) const { return !(*this == other); }

  std::pair<Dart<D>, bool> operator*() const {
    bool orientationPreserving = (iMap + 1) < k;

    Dart<D> currE;
    if (iMap < D) {
      currE = startDart.partner(iMap);
    } else {
      currE = startDart;
    }

    return std::make_pair(currE, orientationPreserving);
  }
};

// Special case for k = 0
template <size_t D>
struct OrbitNeighborhoodIterator<D, 0> {
  Dart<D> startDart;
  size_t iMap, jMap;
  Dart<D> iPartner;

  OrbitNeighborhoodIterator(Dart<D> d) : startDart(d), iMap(1), jMap(0) {
    iPartner = startDart.partner(iMap);
    while (!isValid() && !finished()) advance();
  }

  void advance() {
    jMap++;
    if (jMap >= iMap) {
      jMap = 0;
      iMap++;
      if (iMap < D) iPartner = startDart.partner(iMap);
    }
  }

  bool isValid() const { return (**this).first != startDart; }

  // any two finished iterators are equal, otherwise compare internals
  bool finished() const { return iMap >= D; }

  const OrbitNeighborhoodIterator<D, 0>& operator++() {
    do {
      advance();
    } while (!isValid() && !finished());
    return *this;
  }

  bool operator==(const OrbitNeighborhoodIterator<D, 0>& other) const {
    return (finished() && other.finished()) ||
           (startDart == other.startDart && iMap == other.iMap && jMap == other.jMap);
  }

  bool operator!=(const OrbitNeighborhoodIterator<D, 0>& other) const { return !(*this == other); }

  std::pair<Dart<D>, bool> operator*() const {
    Dart<D> currE;
    if (iPartner != startDart) { // not INVALID_IND
      currE = iPartner.partner(jMap);
    } else {
      currE = startDart;
    }
    return std::make_pair(currE, true); // orientation always preserving for k=0
  }
};

template <size_t D, size_t k>
class OrbitNeighborhood {
public:
  OrbitNeighborhood(Dart<D> d) : dStart(d), cachedEnd(d) { cachedEnd.iMap = D + 1; }

  OrbitNeighborhoodIterator<D, k> begin() const { return OrbitNeighborhoodIterator<D, k>(dStart); }

  // since cachedEnd.finished() == true, checking equality with cachedEnd checks if an iterator is finished()
  OrbitNeighborhoodIterator<D, k> end() const { return cachedEnd; }

private:
  Dart<D> dStart;
  OrbitNeighborhoodIterator<D, k> cachedEnd;
};

} // namespace combinatorial_map
} // namespace geometrycentral
