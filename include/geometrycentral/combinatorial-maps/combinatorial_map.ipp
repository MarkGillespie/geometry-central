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
inline size_t CombinatorialMap<D>::nCells(size_t k) const {
  return nCellsCount[k];
}


template <size_t D>
template <size_t k1, size_t k2>
inline size_t CombinatorialMap<D>::nIncidences() const { // WARNING: if incidences have not been used, returns 0
  static_assert(k1 < k2, "an incidence must have cell dimensions k1 < k2");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");
  auto it = nIncidencesCount.find(std::make_pair(k1, k2));
  return it == nIncidencesCount.end() ? 0 : it->second;
}

template <size_t D>
template <size_t k>
std::vector<Dart<D>> CombinatorialMap<D>::adjacentDarts(Cell<k, D> cell) const {
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

    for (std::pair<Dart<D>, bool> n : orbitNeighbors(curr, k)) {
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
std::vector<Dart<D>> CombinatorialMap<D>::adjacentDarts(Incidence<k1, k2, D> incidence) const {
  // to find all adjacent darts, we express the input cell as an orbit of dart maps,
  std::vector<Dart<D>> neighbors;
  neighbors.reserve(16); // reserve some amount of space

  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(incidence.dart());
  neighbors.push_back(incidence.dart());

  while (!dartsToVisit.empty()) {
    // for some reason, iterating in DFS order is important for orienting cells
    Dart<D> curr = dartsToVisit.back();
    dartsToVisit.pop_back();

    for (Dart<D> n : incidenceNeighboringDarts(curr, k1, k2)) {
      if (std::find(neighbors.begin(), neighbors.end(), n) == neighbors.end()) {
        neighbors.push_back(n);
        dartsToVisit.push_back(n);
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
  dartsToVisit.push_back(std::make_pair(cell.dart(), cell.orientation()));
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

    for (std::pair<Dart<D>, bool> n : orbitNeighbors(currDart, k1)) {
      if (std::find(seenDarts.begin(), seenDarts.end(), n.first) == seenDarts.end()) {
        dartsToVisit.push_back(std::make_pair(n.first, currOrientation == n.second));
        seenDarts.push_back(n.first);
      }
    }
  }

  return neighbors;
}

template <size_t D>
template <size_t k1, size_t k2, size_t k>
std::vector<Cell<k, D>> CombinatorialMap<D>::adjacentCells(Incidence<k1, k2, D> incidence) const {
  // to find all adjacent darts, we express the input cell as an orbit of dart maps,
  std::vector<Cell<k, D>> neighbors;
  std::vector<Dart<D>> seenDarts;
  neighbors.reserve(16); // reserve some amount of space
  seenDarts.reserve(16); // reserve some amount of space

  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(incidence.dart());
  seenDarts.push_back(incidence.dart());

  while (!dartsToVisit.empty()) {
    // for some reason, iterating in DFS order is important for orienting cells
    Dart<D> curr = dartsToVisit.back();
    dartsToVisit.pop_back();

    Cell<k, D> currCell = curr.template cell<k>();
    if (std::find(neighbors.begin(), neighbors.end(), currCell) == neighbors.end()) {
      neighbors.push_back(currCell);
    }

    for (Dart<D> n : incidenceNeighboringDarts(curr, k1, k2)) {
      if (std::find(seenDarts.begin(), seenDarts.end(), n) == seenDarts.end()) {
        dartsToVisit.push_back(n);
        seenDarts.push_back(n);
      }
    }
  }

  if (k == 0 && k1 == 1) { // since we don't use implicit twin, this one case of the tip vertex of a wedge can be missed
    Cell<k, D> tipVertex = incidence.dart().next().template cell<k>();
    if (std::find(neighbors.begin(), neighbors.end(), tipVertex) == neighbors.end()) neighbors.push_back(tipVertex);
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

// Returns a dart in c1 which is also in c2, or Dart<D>() if no such dart can be found
template <size_t D>
template <size_t k1, size_t k2>
Dart<D> CombinatorialMap<D>::adjacentDartInCell(Cell<k1, D> c1, Cell<k2, D> c2) const {
  static_assert(k1 <= D, "cell dimension k1 must be less than or equal to complex dimension D");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");

  // to find all adjacent k2-cells, we express the input cell as an orbit of dart maps,
  // and call d.cell<k2>() for each of these darts
  std::vector<Dart<D>> seenDarts;
  seenDarts.reserve(16); // reserve some amount of space

  std::deque<Dart<D>> dartsToVisit;
  dartsToVisit.push_back(c1.dart());
  seenDarts.push_back(c1.dart());

  while (!dartsToVisit.empty()) {
    Dart<D> currDart = dartsToVisit.back();
    dartsToVisit.pop_back();

    if (currDart.template cell<k2>() == c2) return currDart;

    for (std::pair<Dart<D>, bool> n : orbitNeighbors(currDart, k1)) {
      if (std::find(seenDarts.begin(), seenDarts.end(), n.first) == seenDarts.end()) {
        dartsToVisit.push_back(n.first);
        seenDarts.push_back(n.first);
      }
    }
  }

  return Dart<D>();
}

// OrderedIncidence<a, b> is just an ordinary incidence, but with a and b ordered properly, i.e. Incidence<a,b> if a <
// b and Incidence<b,a> otherwise
template <size_t D>
template <size_t k1, size_t k2>
std::vector<OrderedIncidence<k1, k2, D>> CombinatorialMap<D>::adjacentIncidences(Cell<k1, D> cell) {
  static_assert(k1 <= D, "input cell dimension k1 must be less than or equal to complex dimension D");
  static_assert(k2 <= D, "incident cell dimension k2 must be less than or equal to complex dimension D");

  std::pair<size_t, size_t> key = std::minmax(k1, k2);
  ensureHaveIncidences(k1, k2);

  // to find all adjacent k2-cells, we express the input cell as an orbit of dart maps,
  // and call d.cell<k2>() for each of these darts
  std::vector<OrderedIncidence<k1, k2, D>> neighbors;
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
    OrderedIncidence<k1, k2, D> neighbor(this, dIncidenceArr[key][currDart.getIndex()]);
    if (std::find(neighbors.begin(), neighbors.end(), neighbor) == neighbors.end()) {
      neighbors.push_back(neighbor);
    }

    for (std::pair<Dart<D>, bool> n : orbitNeighbors(currDart, k1)) {
      if (std::find(seenDarts.begin(), seenDarts.end(), n.first) == seenDarts.end()) {
        dartsToVisit.push_back(std::make_pair(n.first, currOrientation == n.second));
        seenDarts.push_back(n.first);
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

template <size_t D>
inline size_t CombinatorialMap<D>::nCells() const {
  return nCells<3>();
}

template <size_t D>
inline size_t CombinatorialMap<D>::nVertexCorners() const {
  return nIncidences<0, D>();
}

template <size_t D>
inline size_t CombinatorialMap<D>::nEdgeCorners() const {
  return nIncidences<1, D>();
}

template <size_t D>
inline size_t CombinatorialMap<D>::nFaceCorners() const {
  return nIncidences<0, 2>();
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

template <size_t D>
template <size_t k1, size_t k2>
inline size_t CombinatorialMap<D>::nIncidencesCapacity() const {
  auto it = nIncidencesCapacityCount.find(std::make_pair(k1, k2));
  return it == nIncidencesCapacityCount.end() ? 0 : it->second;
}

// Connectivity
template <size_t D>
inline size_t CombinatorialMap<D>::dartPartner(size_t iDart, size_t dim) const {
  assert(dim <= D);
  return iDart == INVALID_IND ? INVALID_IND : dartMap[dim][iDart];
}

template <size_t D>
inline void CombinatorialMap<D>::ensureHaveIncidences(size_t k1, size_t k2) {
  if (dIncidenceArr.find(std::make_pair(k1, k2)) == dIncidenceArr.end()) {
    indexIncidences(k1, k2);
  }
}

template <size_t D>
const std::array<std::vector<size_t>, D>& CombinatorialMap<D>::getDartMap() const {
  return dartMap;
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
OrbitNeighborhood<D> CombinatorialMap<D>::orbitNeighbors(Dart<D> d, size_t k) const {
  return OrbitNeighborhood<D>(d, k);
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

template <size_t D> // ensure we have space for n more darts
void CombinatorialMap<D>::allocateDarts(size_t n) {
  // The boring case, when no resize is needed
  while (nDartsFillCount + n >= nDartsCapacityCount) {
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

  nDartsFillCount += n;
  nDartsCount += n;

  modificationTick++;
}

template <size_t D>
template <size_t k>
Cell<k, D> CombinatorialMap<D>::getNewCell() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  if (nCellsFillCount[k] < nCellsCapacityCount[k]) { // The boring case, when no resize is needed
  } else {                                           // The intesting case, where vectors resize
    size_t newCellCapacity = std::max(nCellsCapacityCount[k] * 2, (size_t)1);

    // Resize internal arrays
    cDartArr[k].resize(newCellCapacity);
    nCellsCapacityCount[k] = newCellCapacity;

    // Invoke relevant callback functions
    for (auto& f : cellExpandCallbackList[k]) f(newCellCapacity);
  }

  nCellsFillCount[k]++;
  nCellsCount[k]++;

  modificationTick++;
  return Cell<D, k>(this, nCellsFillCount[k] - 1);
}

template <size_t D>
size_t CombinatorialMap<D>::getNewCellIndex(size_t k) {
  if (nCellsFillCount[k] < nCellsCapacityCount[k]) { // The boring case, when no resize is needed
  } else {                                           // The intesting case, where vectors resize
    size_t newCellCapacity = std::max(nCellsCapacityCount[k] * 2, (size_t)1);

    // Resize internal arrays
    cDartArr[k].resize(newCellCapacity);
    nCellsCapacityCount[k] = newCellCapacity;

    // Invoke relevant callback functions
    for (auto& f : cellExpandCallbackList[k]) f(newCellCapacity);
  }

  nCellsFillCount[k]++;
  nCellsCount[k]++;

  modificationTick++;
  return nCellsFillCount[k] - 1;
}

template <size_t D>
size_t CombinatorialMap<D>::getNewIncidenceIndex(std::pair<size_t, size_t> k1k2) {
  if (nIncidencesFillCount[k1k2] < nIncidencesCapacityCount[k1k2]) { // The boring case, when no resize is needed
  } else {                                                           // The intesting case, where vectors resize
    size_t newIncidenceCapacity = std::max(nIncidencesCapacityCount[k1k2] * 2, (size_t)1);

    // Resize internal arrays
    iDartArr[k1k2].resize(newIncidenceCapacity);
    nIncidencesCapacityCount[k1k2] = newIncidenceCapacity;

    // Invoke relevant callback functions
    for (auto& f : incidenceExpandCallbackList[k1k2]) f(newIncidenceCapacity);
  }

  nIncidencesFillCount[k1k2]++;
  nIncidencesCount[k1k2]++;

  modificationTick++;
  return nIncidencesFillCount[k1k2] - 1;
}

template <size_t D> // ensure we have space for n more k-cells
void CombinatorialMap<D>::allocateCells(size_t k, size_t n) {
  while (nCellsFillCount[k] + n > nCellsCapacityCount[k]) { // Resize vectors if necessary
    size_t newCellCapacity = std::max(nCellsCapacityCount[k] * 2, (size_t)1);

    // Resize internal arrays
    cDartArr[k].resize(newCellCapacity);
    nCellsCapacityCount[k] = newCellCapacity;

    // Invoke relevant callback functions
    for (auto& f : cellExpandCallbackList[k]) f(newCellCapacity);
  }

  nCellsFillCount[k] += n;
  nCellsCount[k] += n;

  modificationTick++;
}

template <size_t D>
template <size_t k> // ensure we have space for n more k-cells
void CombinatorialMap<D>::allocateCells(size_t n) {
  while (nCellsFillCount[k] + n > nCellsCapacityCount[k]) { // Resize vectors if necessary
    size_t newCellCapacity = std::max(nCellsCapacityCount[k] * 2, (size_t)1);

    // Resize internal arrays
    cDartArr[k].resize(newCellCapacity);
    nCellsCapacityCount[k] = newCellCapacity;

    // Invoke relevant callback functions
    for (auto& f : cellExpandCallbackList[k]) f(newCellCapacity);
  }

  nCellsFillCount[k] += n;
  nCellsCount[k] += n;

  modificationTick++;
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
template <size_t k1, size_t k2>
inline bool CombinatorialMap<D>::incidenceIsDead(size_t iI) const {
  static_assert(k1 < k2, "an incidence must have cell dimensions k1 < k2");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");
  auto it = iDartArr.find(std::make_pair(k1, k2));
  return it == iDartArr.end() || it->second[iI] == INVALID_IND;
}

template <size_t D>
inline bool CombinatorialMap<D>::dartIsDead(size_t iDart) const {
  // // a dart is dead if its lowest and highest dimension partner are both INVALID_IND
  // // because we allow duals of meshes with boundary, dartMap[0] might be INVALID_IND even for a live dart
  // return dartMap[0][iD] == INVALID_IND && dartMap[D - 1][iD] == INVALID_IND;
  return iDart == INVALID_IND || dartMap[0][iDart] == INVALID_IND;
}

// Methods for iterating over mesh elements w/ range-based for loops ===========

template <size_t D>
inline DartSet<D> CombinatorialMap<D>::darts() {
  return DartSet<D>(this, 0, nDartsFillCount);
}

template <size_t D>
template <size_t k>
inline CellSet<k, D> CombinatorialMap<D>::cells() {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return CellSet<k, D>(this, 0, nCellsFillCount[k]);
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
inline CellSet<3, D> CombinatorialMap<D>::cells() {
  return cells<3>();
}

template <size_t D>
template <size_t k1, size_t k2>
inline IncidenceSet<k1, k2, D> CombinatorialMap<D>::incidences() {
  static_assert(k1 < k2, "an incidence must have cell dimensions k1 < k2");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");
  ensureHaveIncidences(k1, k2);
  return IncidenceSet<k1, k2, D>(this, 0, nIncidencesFillCount[std::make_pair(k1, k2)]);
}

template <size_t D>
inline IncidenceSet<0, D, D> CombinatorialMap<D>::vertexCorners() {
  return incidences<0, D>();
}

template <size_t D>
inline IncidenceSet<1, D, D> CombinatorialMap<D>::edgeCorners() {
  return incidences<1, D>();
}

template <size_t D>
inline IncidenceSet<0, 2, D> CombinatorialMap<D>::faceCorners() {
  return incidences<0, 2>();
}

// Methods for accessing elements by index =====================================
// Note that these are only valid when the mesh is compressed.

template <size_t D>
inline Dart<D> CombinatorialMap<D>::dart(size_t index) {
  return Dart<D>(this, index);
}

template <size_t D>
template <size_t k>
inline Cell<k, D> CombinatorialMap<D>::cell(size_t index) {
  return Cell<k, D>(this, index);
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
inline Cell<3, D> CombinatorialMap<D>::cell(size_t index) {
  return cell<3>(index);
}

template <size_t D>
template <size_t k1, size_t k2>
inline Incidence<k1, k2, D> CombinatorialMap<D>::incidence(size_t index) {
  return Incidence<k1, k2, D>(this, index);
}

template <size_t D>
inline Incidence<0, D, D> CombinatorialMap<D>::vertexCorner(size_t index) {
  return incidence<0, D>(this, index);
}
template <size_t D>
inline Incidence<1, D, D> CombinatorialMap<D>::edgeCorner(size_t index) {
  return incidence<1, D>(this, index);
}
template <size_t D>
inline Incidence<0, 2, D> CombinatorialMap<D>::faceCorner(size_t index) {
  return incidence<0, 2>(this, index);
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
CellData<3, D, size_t> CombinatorialMap<D>::getCellIndices() {
  return getCellIndices<3>();
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
template <size_t k1, size_t k2>
IncidenceData<k1, k2, D, size_t> CombinatorialMap<D>::getIncidenceIndices() {
  static_assert(k1 < k2, "an incidence must have cell dimensions k1 < k2");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");
  ensureHaveIncidences(k1, k2);
  IncidenceData<k1, k2, D, size_t> indices(*this);
  size_t i = 0;
  for (Incidence<k1, k2, D> inc : incidences<k1, k2>()) {
    indices[inc] = i;
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

// TODO: implement me
// template <size_t D>
// void CombinatorialMap<D>::compress() {}

template <size_t D>
CombinatorialMap<D>::~CombinatorialMap() {
  for (auto& f : meshDeleteCallbackList) {
    f();
  }
}

template <size_t D>
std::unique_ptr<CombinatorialMap<D>> CombinatorialMap<D>::copy() const {
  CombinatorialMap<D>* newMesh = new CombinatorialMap<D>();
  copyInternalFields(*newMesh);
  return std::unique_ptr<CombinatorialMap<D>>(newMesh);
}

template <size_t D>
std::unique_ptr<CombinatorialMap<D>> CombinatorialMap<D>::dual() const {
  const bool DEBUG_PRINT = false;

  if (DEBUG_PRINT) {
    std::cout << "====== computing dual of: " << std::endl;
    std::cout << "#\tnext\ttwin" << std::endl;
    for (size_t iDart = 0; iDart < this->nDartsCount; iDart++) {
      std::cout << iDart << "\t" << this->dartPartner(iDart, 0) << "\t" << this->dartPartner(iDart, 1) << std::endl;
    }
  }

  std::unique_ptr<CombinatorialMap<D>> result(new CombinatorialMap<D>());

  CombinatorialMap& newMesh = *result;

  for (size_t k = 0; k <= D; k++) {
    newMesh.nCellsCount[k] = this->nCellsCount[D - k];
    newMesh.nCellsCapacityCount[k] = newMesh.nCellsCount[k];
    newMesh.nCellsFillCount[k] = newMesh.nCellsCount[k];
  }
  newMesh.nDartsCount = this->nDartsCount;
  newMesh.nDartsCapacityCount = this->nDartsCount;
  newMesh.nDartsFillCount = this->nDartsCount;

  for (size_t k = 0; k <= D; k++) {
    newMesh.dartMap[k].reserve(newMesh.nDartsCount);
    newMesh.dCellArr[k].reserve(newMesh.nDartsCount);
    newMesh.dCellSgn[k].reserve(newMesh.nDartsCount);
    // initialize cDartArr arrays, to be filled while constructing darts
    newMesh.cDartArr[k] = std::vector<size_t>(newMesh.nCellsCount[k], INVALID_IND);
  }

  // Identify boundary vertices, which will get deleted in the dual mesh
  // This code would be simpler if we could loop over this->vertices(), but for annoying technical reasons that function
  // is not const, so we work directly with indices here
  std::vector<bool> isBdyVtx(this->nCellsCount[0], false);
  for (size_t iDart = 0; iDart < this->nDartsCount; iDart++) {
    if (this->dartMap[D - 1][iDart] == INVALID_IND) isBdyVtx[this->dCellArr[0][iDart]] = true;
  }

  if (DEBUG_PRINT) {
    std::cout << "====== isBdyVtx:" << std::endl;
    for (size_t iV = 0; iV < this->nCellsCount[0]; iV++) {
      std::cout << "  " << iV << " : " << (isBdyVtx[iV] ? "true" : "false") << std::endl;
    }
  }

  // Helper function for computing dart maps in dual mesh
  // If partner lies in an interior cell, return partner, otherwise return INVALID_IND
  // If iDart is INVALID_IND, return INVALID_IND (just like this->dartPartner does)
  auto validatedPartner = [&](size_t iDart, size_t dim) -> size_t {
    size_t partner = this->dartPartner(iDart, dim);
    return (partner == INVALID_IND || isBdyVtx[this->dCellArr[0][partner]]) ? INVALID_IND : partner;
  };

  // Construct dart and cell indices/signs on dual mesh
  for (size_t iDart = 0; iDart < newMesh.nDartsCount; iDart++) {
    if (isBdyVtx[dCellArr[0][iDart]]) { // mark darts on boundary vertices as dead
      for (size_t k = 0; k < D; k++) newMesh.dartMap[k].push_back(INVALID_IND);
      for (size_t k = 0; k <= D; k++) {
        newMesh.dCellArr[k].push_back(INVALID_IND);
        newMesh.dCellSgn[k].push_back(1);
      }
      continue;
    }

    // Otherwise, apply standard construction of dual maps:

    // On the dual mesh, we set
    // β*_0 = β_{D-2} o β_{D-1}
    // β*_k = β_{D-k-2} o β_{D-1}
    // β*_{D-1} = β_{D-1}
    newMesh.dartMap[0].push_back(validatedPartner(this->dartPartner(iDart, D - 1), D - 2));
    for (size_t k = 1; k < D - 1; k++) {
      newMesh.dartMap[k].push_back(validatedPartner(this->dartPartner(iDart, D - 1), D - k - 2));
    }
    newMesh.dartMap[D - 1].push_back(validatedPartner(iDart, D - 1));

    if (DEBUG_PRINT) {
      std::cout << "processing dart " << iDart << std::endl;
      for (size_t k = 0; k < D; k++) {
        std::cout << "   set dartMap[" << k << "] : " << newMesh.dartMap[k][iDart] << std::endl;
      }
      std::cout << "         > β_{" << D << "-1} = dart " << this->dartPartner(iDart, D - 1) << std::endl;
      std::cout << "         > β_{" << D << "-2} o β_{" << D << "-1} = dart "
                << this->dartPartner(this->dartPartner(iDart, D - 1), D - 2) << std::endl;
      std::cout << "         > β_{" << D << "-2} = dart " << this->dartPartner(iDart, D - 2) << std::endl;
      for (size_t k = 2; k < D; k++) {
        std::cout << "       > β_{" << D << "-1} = dart " << this->dartPartner(iDart, D - 1) << std::endl;
        std::cout << "       > β_{" << D << "-k-2} o β_{" << D << "-1} = dart "
                  << this->dartPartner(this->dartPartner(iDart, D - 1), D - k - 2) << std::endl;
      }
    }

    for (size_t k = 0; k <= D; k++) {
      size_t iCell = this->dCellArr[D - k][iDart];
      int sgn = this->dCellSgn[D - k][iDart];
      newMesh.dCellArr[k].push_back(iCell);
      newMesh.dCellSgn[k].push_back(sgn);
      if (sgn == 1 || newMesh.cDartArr[k][iCell] == INVALID_IND) newMesh.cDartArr[k][iCell] = iDart;
    }
  }

  for (size_t k = 0; k <= D; k++) {
    for (size_t iC = 0; iC < this->nCellsCount[D - k]; iC++) {
      // If possible, copy over cell.dart() from primal mesh. But if that dart is dead, leave the valid dart that we
      // found earlier
      size_t oldCellDart = this->cDartArr[D - k][iC];
      if (!newMesh.dartIsDead(oldCellDart)) newMesh.cDartArr[k][iC] = oldCellDart;
      if (newMesh.dartIsDead(newMesh.cDartArr[k][iC])) { // filter out dead darts (from boundary cells)
        newMesh.nCellsCount[k]--;
        newMesh.cDartArr[k][iC] = INVALID_IND;
      }
    }
  }

  return std::move(result);
}

template <size_t D>
void CombinatorialMap<D>::copyInternalFields(CombinatorialMap<D>& target) const {
  // == Copy _all_ the fields!

  // Raw data buffers (underlying std::vectors duplicate storage automatically)
  target.dartMap = dartMap;
  target.dCellArr = dCellArr;
  target.dCellSgn = dCellSgn;
  target.cDartArr = cDartArr;

  // counts and flags
  target.nDartsCount = nDartsCount;
  target.nDartsCapacityCount = nDartsCapacityCount;
  target.nDartsFillCount = nDartsFillCount;
  target.nCellsCount = nCellsCount;
  target.nCellsCapacityCount = nCellsCapacityCount;
  target.nCellsFillCount = nCellsFillCount;

  target.isCompressedFlag = isCompressedFlag;

  // Note: _don't_ copy callbacks lists! New mesh has new callbacks
}

// index k-cells and fill cDartArr[k] and dCellArr[k] based off of dartMap
template <size_t D>
void CombinatorialMap<D>::indexCells(size_t k) {
  using namespace unionfind;
  if (k > D) return;
  // static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
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

  for (Dart<D> d : darts()) {
    for (std::pair<Dart<D>, bool> n : orbitNeighbors(d, k)) {
      unite(d.getIndex(), n.first.getIndex(), n.second, parent, sharesParentSign, rank);
      if (DEBUG_PRINT) {
        std::cout << "uniting dart " << d.getIndex() << " with dart " << n.first.getIndex()
                  << " | orientationPreserving: " << (n.second ? "true" : "false") << std::endl;
      }
    }
  }

  if (DEBUG_PRINT) {
    std::cout << std::endl << "Final " << k << "-cell orbits: " << std::endl;
    for (size_t iRoot = 0; iRoot < nDarts(); iRoot++) {
      if (findRoot(iRoot, parent, sharesParentSign) != iRoot) continue;
      std::cout << "  root " << iRoot << std::endl;
      for (size_t iDart = 0; iDart < nDarts(); iDart++) {
        if (findRoot(iDart, parent, sharesParentSign) == iRoot) {
          std::cout << "     dart " << iDart << " : " << dCellArr[0][iDart] << "->" << dCellArr[0][dartMap[0][iDart]]
                    << "\tpositive orientation: " << (sharesParentSign[iDart] ? "true" : "false") << std::endl;
        }
      }
      std::cout << std::endl;
    }
  }

  // Clear any existing k-cells
  dCellArr[k].resize(nDarts());
  std::fill(dCellArr[k].begin(), dCellArr[k].end(), INVALID_IND);
  dCellSgn[k].resize(nDarts());
  std::fill(dCellSgn[k].begin(), dCellSgn[k].end(), true);
  nCellsFillCount[k] = 0;
  nCellsCount[k] = 0;

  // Assign new k-cell indices
  for (size_t iDart = 0; iDart < nDarts(); iDart++) {
    size_t iRoot = findRoot(iDart, parent, sharesParentSign);
    if (dCellArr[k][iRoot] == INVALID_IND) {
      dCellArr[k][iRoot] = getNewCellIndex(k);
      dCellSgn[k][iRoot] = true;
      cDartArr[k][dCellArr[k][iRoot]] = iRoot;
    }
    dCellArr[k][iDart] = dCellArr[k][iRoot];
    dCellSgn[k][iDart] = sharesParentSign[iDart];
  }

  // Shrink internal arrays
  cDartArr[k].resize(nCellsCount[k]);
  nCellsCapacityCount[k] = nCellsCount[k];
}

// index k-cells and fill cDartArr[k] and dCellArr[k] based off of dartMap
template <size_t D>
void CombinatorialMap<D>::indexIncidences(size_t k1, size_t k2) {
  using namespace unionfind;
  if (k1 > D || k2 > D) return;
  const bool DEBUG_PRINT = false;

  std::vector<size_t> parent;
  parent.reserve(nDarts()); // initialize every dart as its own parent
  for (size_t i = 0; i < nDarts(); i++) parent.push_back(i);

  std::vector<size_t> rank(nDarts(), 0); // initialize every dart to rank 0
  std::vector<bool> sharesParentSign(nDarts(), true);

  for (Dart<D> d : darts()) {
    for (Dart<D> n : incidenceNeighboringDarts(d, k1, k2)) {
      unite(d.getIndex(), n.getIndex(), true, parent, sharesParentSign, rank);
      if (DEBUG_PRINT) {
        std::cout << "uniting dart " << d.getIndex() << " with dart " << n.getIndex() << std::endl;
      }
    }
  }

  if (DEBUG_PRINT) {
    std::cout << std::endl << "Final (" << k1 << ", " << k2 << ")-incidence orbits: " << std::endl;
    for (size_t iRoot = 0; iRoot < nDarts(); iRoot++) {
      if (findRoot(iRoot, parent, sharesParentSign) != iRoot) continue;
      std::cout << "  root " << iRoot << std::endl;
      for (size_t iDart = 0; iDart < nDarts(); iDart++) {
        if (findRoot(iDart, parent, sharesParentSign) == iRoot) {
          std::cout << "     dart " << iDart << " : " << dCellArr[0][iDart] << "->" << dCellArr[0][dartMap[0][iDart]]
                    << std::endl;
        }
      }
      std::cout << std::endl;
    }
  }

  // Clear any existing k1,k2-incidences
  std::pair<size_t, size_t> key = std::make_pair(k1, k2);
  dIncidenceArr[key] = std::vector<size_t>(nDartsCapacity(), INVALID_IND);
  nIncidencesFillCount[key] = 0;
  nIncidencesCount[key] = 0;

  // Assign new incidence indices
  for (size_t iDart = 0; iDart < nDarts(); iDart++) {
    size_t iRoot = findRoot(iDart, parent, sharesParentSign);
    if (dIncidenceArr[key][iRoot] == INVALID_IND) {
      dIncidenceArr[key][iRoot] = getNewIncidenceIndex(key);
      iDartArr[key][dIncidenceArr[key][iRoot]] = iRoot;
    }
    dIncidenceArr[key][iDart] = dIncidenceArr[key][iRoot];
  }

  // Shrink internal arrays
  iDartArr[key].resize(nIncidencesCount[key]);
  nIncidencesCapacityCount[key] = nIncidencesCount[key];
}

// computes n! / 2, used as a helper function for simplicial complex constructor
constexpr size_t halfFactorial(size_t n) { return (n <= 2) ? 1 : n * halfFactorial(n - 1); }

// generate a list of the n!/2 positive permutations on [0, ..., n-1], listed in lexicographic order
template <size_t n>
std::array<std::array<size_t, n>, halfFactorial(n)> listPositivePermutations() {
  std::array<std::array<size_t, n>, halfFactorial(n)> result;

  std::array<size_t, n> perm; // Create initial permutation [0, ..., n-1]
  for (size_t i = 0; i < n; ++i) perm[i] = i;

  // Helper function to compute sign of permutation
  auto computeSign = [](const std::array<size_t, n>& p) -> int {
    size_t inversions = 0; // Count inversions: pairs (i,j) where i < j but p[i] > p[j]
    for (size_t i = 0; i < n - 1; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (p[i] > p[j]) ++inversions;
      }
    }

    return (inversions % 2 == 0) ? 1 : -1; // Sign is +1 if even number of inversions, -1 if odd
  };

  // List all permutations and keep only positive ones
  size_t iP = 0; // index of current permutation
  do {
    if (computeSign(perm) == 1) result[iP++] = perm;
  } while (std::next_permutation(perm.begin(), perm.end()));

  return result;
}

// construct dart maps for D-simplex with our chosen dart ordering
// D-simplex has D-1 dart maps for its (D+1)!/2 darts
template <size_t D>
std::array<std::array<size_t, halfFactorial(D + 1)>, D - 1>
simplexDartMaps(const std::array<std::array<size_t, D + 1>, halfFactorial(D + 1)>& positivePermutations) {
  auto permIndex = [&positivePermutations](const std::array<size_t, D + 1>& perm) -> size_t {
    // Use std::lower_bound with lexicographic comparison
    auto it = std::lower_bound(positivePermutations.begin(), positivePermutations.end(), perm);

    // Check if we found the permutation
    if (it != positivePermutations.end() && *it == perm) return std::distance(positivePermutations.begin(), it);

    return INVALID_IND; // Permutation not found
  };

  auto shiftFirstThreeIndices = [](const std::array<size_t, D + 1>& perm) -> std::array<size_t, D + 1> {
    std::array<size_t, D + 1> result = perm;
    result[0] = perm[1], result[1] = perm[2], result[2] = perm[0];
    return result;
  };

  auto swapIndexPairs = [](std::array<size_t, D + 1> perm, std::array<size_t, 2> i,
                           std::array<size_t, 2> j) -> std::array<size_t, D + 1> {
    std::swap(perm[i[0]], perm[i[1]]);
    std::swap(perm[j[0]], perm[j[1]]);
    return perm;
  };

  std::array<std::array<size_t, halfFactorial(D + 1)>, D - 1> dartMaps;

  for (size_t iDart = 0; iDart < halfFactorial(D + 1); iDart++) {
    // next map: cyclic shift first 3 indices
    dartMaps[0][iDart] = permIndex(shiftFirstThreeIndices(positivePermutations[iDart]));

    // other maps: swap two pairs of indices
    for (size_t dim = 1; dim < D - 1; dim++) {
      dartMaps[dim][iDart] = permIndex(swapIndexPairs(positivePermutations[iDart], {0, 1}, {dim + 1, dim + 2}));
    }
  }

  return dartMaps;
}

//=== Explicit specializations of index arrays for 2- and 3-dimensional combinatorial maps
// More can be generated using the helper code in the CombinatorialMap<D> constructor below
// clang-format off
template <> // listPositivePermutations<3>(), used in CombinatorialMap<2>
constexpr std::array<std::array<size_t, 3>,3> listPositivePermutations<3>() {
  return {{ {0, 1, 2}, {1, 2, 0}, {2, 0, 1} }};
}

template <> // listPositivePermutations<4>(), used in CombinatorialMap<3>
constexpr std::array<std::array<size_t, 4>,12> listPositivePermutations<4>() {
  return {{ {0, 1, 2, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 3, 2, 0}, {2, 0, 1, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 2, 1, 0} }};
}

template <> // simplexDartMaps<2>(), used in CombinatorialMap<2>
constexpr std::array<std::array<size_t, 3>, 1> simplexDartMaps<2>(const std ::array<std::array<size_t, 3>, 3>& _) {
  (void)_; // tell the compiler not to complain that _ is unused
  return {{ // for some reason, C++11 wants double braces for std::array
      {1, 2, 0} // dartMap[0]
  }};
}

template <> // simplexDartMaps<3>(), used in CombinatorialMap<3>
constexpr std::array<std::array<size_t, 12>, 2> simplexDartMaps<3>(const std::array<std::array<size_t, 4>, 12>& _) {
  (void)_; // tell the compiler not to complain that _ is unused
  return {{  // for some reason, C++11 wants double braces for std::array
      {4, 8, 10, 2, 6, 11, 0, 5, 9, 1, 3, 7}, // dartMap[0]
      {3, 6, 9, 0, 7, 10, 1, 4, 11, 2, 5, 8}  // dartMap[1]
  }};
}
// clang-format on

template <size_t D>
CombinatorialMap<D>::CombinatorialMap(const std::vector<std::array<size_t, D + 1>>& simplices) {
  constexpr bool DEBUG_PRINT = false;

  // darts on a D-simplex are in 1-1 correspondence with positive permutations on D+1 elements
  // the number of darts in a D-simplex is (D+1)! / 2:
  // a k-simplex has (k+1) top-dimensional faces, each of which has (k-1) top-dimensional faces, all the way down to a
  // 1-simplex which has 1 dart
  constexpr size_t nSimplexDarts = halfFactorial(D + 1);
  std::array<std::array<size_t, D + 1>, nSimplexDarts> positivePermutations = listPositivePermutations<D + 1>();
  if (DEBUG_PRINT) {
    std::cout << "===== constructing " << D
              << "-dimensional simplicial complex. Dart maps are constructed from the following list of positive "
                 "permutations:"
              << std::endl;
    for (size_t iP = 0; iP < positivePermutations.size(); iP++) {
      std::cout << "   permutation " << std::setw(3) << iP << ": { ";
      for (size_t i = 0; i < positivePermutations[iP].size(); i++) {
        std::cout << positivePermutations[iP][i] << ((i + 1 < positivePermutations[iP].size()) ? ", " : " ");
      }
      std::cout << "}" << std::endl;
    }
  }

  nCellsCount[0] = 0;
  for (const std::array<size_t, D + 1>& simplex : simplices) {
    for (size_t i : simplex) {
      nCellsCount[0] = std::max(nCellsCount[0], i);
    }
  }
  nCellsCount[0]++; // 0-based means count is max + 1

  cDartArr[0] = std::vector<size_t>(nCellsCount[0], INVALID_IND);
  const std::array<std::array<size_t, nSimplexDarts>, D - 1> faceDartMaps = simplexDartMaps<D>(positivePermutations);

  constexpr bool printPositivePermutations = false; // if true, print template specialization so it can be used directly
  constexpr bool printSimplexMaps = false;          // if true, print template specialization so it can be used directly
  if (printPositivePermutations) {
    std::cout << "template <> // listPositivePermutations<" << (D + 1) << ">(), used in CombinatorialMap<" << D << ">"
              << std::endl;
    std::cout << "constexpr std::array<std::array<size_t, " << (D + 1) << ">," << nSimplexDarts
              << "> listPositivePermutations<" << D + 1 << ">() {" << std::endl;
    std::cout << "\treturn {{ ";
    for (size_t iP = 0; iP < nSimplexDarts; iP++) {
      std::cout << "{";
      for (size_t i = 0; i < D + 1; i++) {
        std::cout << positivePermutations[iP][i];
        if (i + 1 < D + 1) std::cout << ", ";
      }
      std::cout << "}";
      if (iP + 1 < nSimplexDarts) std::cout << ", ";
    }
    std::cout << " }};" << std::endl << "}" << std::endl;
  }
  if (printSimplexMaps) {
    std::cout << "template <> // simplexDartMaps<" << D << ">(), used in CombinatorialMap<" << D << ">" << std::endl;
    std::cout << "constexpr std::array<std::array<size_t, " << nSimplexDarts << ">," << (D - 1) << "> simplexDartMaps<"
              << D << ">(const std::array<std::array<size_t, " << D + 1 << ">, " << nSimplexDarts << ">& _) {"
              << std::endl;
    std::cout << "\treturn {{ // for some reason, C++11 wants double braces for std::array" << std::endl;
    for (size_t dim = 0; dim < D - 1; dim++) {
      std::cout << "\t\t{";
      for (size_t iDart = 0; iDart < nSimplexDarts; iDart++) {
        std::cout << faceDartMaps[dim][iDart];
        if (iDart + 1 < nSimplexDarts) std::cout << ", ";
      }
      std::cout << "}";
      if (dim + 1 < D - 1) std::cout << ",";
      std::cout << " // dartMap[" << dim << "]";
      std::cout << std::endl;
    }
    std::cout << "\t}};" << std::endl << "}" << std::endl;
  }

  //=== Walk over simplices, attaching dart maps
  std::map<std::array<size_t, D>, size_t> createdDarts; // cache darts in codimenision-1 faces
  auto attachTopDartMap = [&](size_t iDart, size_t jDart) -> void {
    dartMap[D - 1][iDart] = jDart;
    dartMap[D - 1][jDart] = iDart;
  };

  for (const std::array<size_t, D + 1>& simplex : simplices) {
    size_t iCellD = getNewCell<D>().getIndex();
    if (DEBUG_PRINT) std::cout << " ... processing input simplex " << iCellD << std::endl;
    std::array<size_t, nSimplexDarts> newDartIndices; // index new darts
    // set dartMap[0..D-2], dCellArr[0 and D], cDartArr[0 and D] for new dart, and dCellSgn[0 and D]
    for (size_t iDart = 0; iDart < nSimplexDarts; iDart++) newDartIndices[iDart] = getNewDart().getIndex();
    for (size_t iDart = 0; iDart < nSimplexDarts; iDart++) {
      size_t newDart = newDartIndices[iDart];
      const std::array<size_t, D + 1>& perm = positivePermutations[iDart]; // permutation representation of dart
      // set dart-cell index maps
      dCellArr[0][newDart] = simplex[perm[0]];
      dCellSgn[0][newDart] = true;
      cDartArr[0][dCellArr[0][newDart]] = newDart;
      dCellArr[D][newDart] = iCellD;
      dCellSgn[D][newDart] = true;
      cDartArr[D][dCellArr[D][newDart]] = newDart;
      dartMap[D - 1][newDart] = INVALID_IND; // initialize top map to INVALID_IND

      for (size_t dim = 0; dim < D - 1; dim++) // set dart maps 0 .. D-2 from faceDartMaps
        dartMap[dim][newDart] = newDartIndices[faceDartMaps[dim][iDart]];
    }
    // take care of dartMap[D-1] in a separate loop so we can use values of dartMap[0]
    for (size_t iDart = 0; iDart < nSimplexDarts; iDart++) {
      size_t newDart = newDartIndices[iDart];
      if (dartMap[D - 1][newDart] == INVALID_IND) {
        const std::array<size_t, D + 1>& perm = positivePermutations[iDart]; // permutation representation of dart
        std::array<size_t, D> key;                                           // search for codimension-1 face to glue to
        for (size_t dim = 0; dim < D; dim++) key[dim] = simplex[perm[dim]];
        auto topTwin = createdDarts.find(key);
        if (topTwin == createdDarts.end()) { // if partner(D-1) has not been created, save dart to map
          std::swap(key[0], key[1]);
          // make sure face has not already been created with same orientation
          auto repeatedFaceDart = createdDarts.find(key);
          if (repeatedFaceDart != createdDarts.end()) {
            std::stringstream ss;
            ss << "Input simplex orientation error: face { ";
            for (size_t i = 0; i < key.size(); i++) ss << key[i] << (i + 1 < key.size() ? ", " : "");
            ss << " } found in both simplex " << iCellD << " and simplex " << dCellArr[D][repeatedFaceDart->second]
               << " with same orientation";

            throw std::runtime_error(ss.str());
          }
          // save dart to map in case we find its partner later
          createdDarts[key] = newDart;
        } else { // otherwise fill in partner(D-1) for this dart, and its next and next.next
          attachTopDartMap(newDart, topTwin->second);
          if (D > 2) { // for D > 2, we can glue the rest of the darts in the 2-face as well
            attachTopDartMap(dartMap[0][newDart], dartMap[0][dartMap[0][topTwin->second]]);
            attachTopDartMap(dartMap[0][dartMap[0][newDart]], dartMap[0][topTwin->second]);
          }
        }
      }
    }
  }


  if (DEBUG_PRINT) {
    for (size_t iD = 0; iD < nDarts(); iD++) {
      std::cout << "Dart " << iD << " : " << dCellArr[0][iD] << "->" << dCellArr[0][dartMap[0][iD]] << std::endl;
    }
    for (size_t dim = 0; dim < D; dim++) {
      for (size_t iDart = 0; iDart < nDarts(); iDart++) {
        std::cout << "dartMap[" << dim << "][" << iDart << "] = " << dartMap[dim][iDart] << std::endl;
      }
      std::cout << std::endl;
    }
  }

  // trim excess capacity
  nCellsCapacityCount[0] = nCellsCount[0];
  nCellsFillCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nDartsFillCount = nDartsCount;
  for (size_t k = 0; k < D; k++) dartMap[k].resize(nDartsCount);

  // Shrink internal arrays for 3-cells (which may be over-sized since they're allocated by doubling)
  cDartArr[D].resize(nCellsCount[D]);
  nCellsCapacityCount[D] = nCellsCount[D];

  // construct intermediate k-cells
  for (size_t k = 1; k < D; k++) indexCells(k);

  if (DEBUG_PRINT) {
    for (size_t k = 0; k <= D; k++) std::cout << "# " << k << "-cells : " << nCellsCount[k] << std::endl;
  }
}

// Construct directly from dart map
template <size_t D>
CombinatorialMap<D>::CombinatorialMap(const std::array<std::vector<size_t>, D>& dartMap_) {
  dartMap = dartMap_;

  nDartsCount = dartMap[0].size();
  nDartsFillCount = nDartsCount;
  nDartsCapacityCount = nDartsCount;

  for (size_t k = 0; k <= D; k++) indexCells(k);

  isCompressedFlag = true;
}

template <size_t D>
static CombinatorialMap<D> CombinatorialMap<D>::Random(size_t nDarts) {
  if (nDarts % 2 == 1)
    throw std::logic_error("CombinatorialMap<D>::Random error: number of darts in a combinatorial map must be even");
  std::array<std::vector<size_t>, D> dartMap;
  dartMap[0] = std::vector<size_t>(nDarts); // Without loss of generality, set adjacent darts to be twins
  for (size_t iDart = 0; iDart < nDarts; iDart += 2) {
    dartMap[0][iDart] = iDart + 1;
    dartMap[0][iDart + 1] = iDart;
  }

  // generate remaining permutations randomly
}

// // TODO: finish this, maybe by constructing boundary matrices
// template <size_t D>
// CombinatorialMap<D>::CombinatorialMap(const NestedVector<D, size_t>& cells) {}

#ifdef SPECIALIZATIONS

// Builds a 2D polygon mesh
template <>
CombinatorialMap<2>::CombinatorialMap(const std::vector<std::vector<size_t>>& polygons) {
  surface::ManifoldSurfaceMesh mesh(polygons);

  nCellsCount[0] = mesh.nVertices(), nCellsCount[1] = mesh.nEdges(), nCellsCount[2] = mesh.nFaces();
  nDartsCount = mesh.nHalfedges();
  for (size_t dim = 0; dim <= 2; dim++) {
    nCellsCapacityCount[dim] = nCellsCount[dim];
    nCellsFillCount[dim] = nCellsCount[dim];
  }
  nDartsCapacityCount = nDartsCount, nDartsFillCount = nDartsCount;

  surface::VertexData<size_t> vIdx = mesh.getVertexIndices();
  surface::EdgeData<size_t> eIdx = mesh.getEdgeIndices();
  surface::FaceData<size_t> fIdx = mesh.getFaceIndices();
  surface::HalfedgeData<size_t> hIdx = mesh.getHalfedgeIndices();

  dartMap[0].reserve(mesh.nHalfedges());
  dartMap[1].reserve(mesh.nHalfedges());
  for (size_t dim = 0; dim <= 2; dim++) {
    dCellArr[dim].reserve(mesh.nHalfedges());
    dCellSgn[dim].reserve(mesh.nHalfedges());
  }
  cDartArr[0].reserve(mesh.nVertices());
  cDartArr[1].reserve(mesh.nEdges());
  cDartArr[2].reserve(mesh.nFaces());
  for (surface::Halfedge he : mesh.halfedges()) {
    // set exterior halfedges to INVALID_IND since our combinatorial maps don't use the implicit twin convention
    dartMap[0].push_back(he.next().isInterior() ? hIdx[he.next()] : INVALID_IND);
    dartMap[1].push_back(he.twin().isInterior() ? hIdx[he.twin()] : INVALID_IND);
    dCellArr[0].push_back(he.isInterior() ? vIdx[he.vertex()] : INVALID_IND);
    dCellArr[1].push_back(he.isInterior() ? eIdx[he.edge()] : INVALID_IND);
    dCellArr[2].push_back(he.isInterior() ? fIdx[he.face()] : INVALID_IND);
    dCellSgn[0].push_back(true);
    dCellSgn[1].push_back(he.orientation());
    dCellSgn[2].push_back(true);
  }
  for (surface::Vertex v : mesh.vertices()) cDartArr[0].push_back(hIdx[v.halfedge()]);
  for (surface::Edge e : mesh.edges()) cDartArr[1].push_back(hIdx[e.halfedge()]);
  for (surface::Face f : mesh.faces()) cDartArr[2].push_back(hIdx[f.halfedge()]);
}


// Builds a 3D volume mesh
template <>
CombinatorialMap<3>::CombinatorialMap(const std::vector<std::vector<std::vector<size_t>>>& cells) {
  const bool DEBUG_PRINT = true;
  nCellsCount[0] = 0;
  for (const std::vector<std::vector<size_t>>& cell : cells) {
    for (const std::vector<size_t>& face : cell) {
      for (size_t i : face) nCellsCount[0] = std::max(nCellsCount[0], i);
    }
  }
  nCellsCount[0]++; // 0-based means count is max + 1

  cDartArr[0] = std::vector<size_t>(nCellsCount[0], INVALID_IND);

  // take in canonicalized vertex list, return first dart index, number of shifts to canonicalize, orientation
  std::map<std::vector<size_t>, std::tuple<size_t, int, bool>> createdDarts;

  // turn input list into canonicalized vertex list, number of shifts, and orientation
  struct CanonicalVertexList {
    std::vector<size_t> vertices; // canonical ordering
    int rotation;                 // how many times to rotate to get from input to canonical
    bool orientation;             // false <=> input was flipped during canonicalization
  };
  auto canonicalize = [](const std::vector<size_t>& face) -> CanonicalVertexList {
    size_t degree = face.size();

    size_t minIdx = 0; // find minimum vertex index
    for (size_t i = 1; i < degree; ++i) {
      if (face[i] < face[minIdx]) minIdx = i;
    }

    // there are two possible orderings starting with the minimum vertex
    // choose the unique ordering where the second vertex is smaller than the last
    std::vector<size_t> result;
    result.reserve(degree);
    if (face[(minIdx + 1) % degree] < face[(minIdx + degree - 1) % degree]) {
      for (size_t i = 0; i < degree; i++) result.push_back(face[(minIdx + i) % degree]); // f[minIdx], f[minIdx+1],...
      return {result, static_cast<int>(minIdx), true};
    } else {
      for (size_t i = 0; i < degree; i++)
        result.push_back(face[(minIdx + degree - i) % degree]); // f[minIdx], f[minIdx-1],...
      return {result, static_cast<int>(minIdx), false};
    }
  };

  auto oppDartLookup = [&](const std::vector<size_t>& face, const CanonicalVertexList& canon) -> size_t {
    size_t degree = face.size();
    auto dartIt = createdDarts.find(canon.vertices);
    if (dartIt == createdDarts.end()) {
      return INVALID_IND; // never seen this face with any orientation
    } else {
      // make sure this face hasn't already appeared with the same orientation
      GC_SAFETY_ASSERT(canon.orientation != std::get<2>(dartIt->second),
                       "tet mesh orientation problem: duplicate face {" + std::to_string(face[0]) + ", " +
                           std::to_string(face[1]) + ", " + std::to_string(face[2]) + "}");
      int storedRotation = std::get<1>(dartIt->second), inputRotation = canon.rotation;
      int rotationDiff = (inputRotation + storedRotation + degree - 1) % degree;
      size_t resultDart = std::get<0>(dartIt->second); // repeatedly apply dartMap[0] to reference dart stored in map
      for (int i = 0; i < rotationDiff; ++i) resultDart = dartMap[0][resultDart];

      return resultDart;
    }
  };

  auto attachDartMap2 = [&](size_t iDart, size_t jDart) -> void {
    dartMap[2][iDart] = jDart;
    dartMap[2][jDart] = iDart;
  };

  for (size_t iC = 0; iC < cells.size(); iC++) {
    const std::vector<std::vector<size_t>>& cell = cells[iC];
    size_t iCell3 = getNewCell<3>().getIndex();

    // construct new darts, set cDartArr[0], dCellArr/dCellSign[0 and 3]
    std::vector<std::vector<size_t>> newDartIndices(cell.size());
    for (size_t iF = 0; iF < cell.size(); iF++) {
      const std::vector<size_t>& face = cell[iF];
      newDartIndices[iF].reserve(face.size());
      for (size_t iD = 0; iD < face.size(); iD++) {
        size_t newDart = getNewDart().getIndex();
        newDartIndices[iF].push_back(newDart);

        cDartArr[0][face[iD]] = newDart;
        dCellArr[0][newDart] = face[iD];
        dCellSgn[0][newDart] = true;
        dCellArr[3][newDart] = iCell3;
        dCellSgn[3][newDart] = true;
      }
    }
    cDartArr[3][iCell3] = newDartIndices[0][0];

    // fill in dartMap[0][newDart] and dartMap[1][newDart] using cell's next and twin maps
    // map std::minmax(dart.tailvertex, dart.tipvertex) -> (dart id, dart orientation), helps construct twin map
    std::map<std::pair<size_t, size_t>, std::pair<size_t, bool>> cellDarts;
    for (size_t iF = 0; iF < cell.size(); iF++) {
      const std::vector<size_t>& face = cell[iF];
      for (size_t iD = 0; iD < face.size(); iD++) {
        size_t dart = newDartIndices[iF][iD], next = newDartIndices[iF][(iD + 1) % face.size()];
        dartMap[0][dart] = next; // set next() to next dart around face
        size_t vTail = dCellArr[0][dart], vTip = dCellArr[0][next];
        std::pair<size_t, size_t> key = std::minmax(vTail, vTip);
        auto itTwin = cellDarts.find(key);
        if (itTwin != cellDarts.end()) { // found another dart along this edge
          GC_SAFETY_ASSERT((vTail < vTip) != itTwin->second.second, "cell orientation error: " + std::to_string(vTail) +
                                                                        "->" + std::to_string(vTip) +
                                                                        " appears twice in cell " + std::to_string(iC));
          dartMap[1][dart] = itTwin->second.first;
          dartMap[1][itTwin->second.first] = dart;
        } else { // first time seeing this edge, save dart
          cellDarts[key] = std::make_pair(dart, (vTail < vTip));
        }
      }
    }

    if (DEBUG_PRINT) {
      for (size_t iF = 0; iF < cell.size(); iF++) {
        const std::vector<size_t>& face = cell[iF];
        for (size_t iD = 0; iD < face.size(); ++iD) {
          std::cout << "Dart " << newDartIndices[iF][iD] << " : " << dCellArr[0][newDartIndices[iF][iD]] << "->"
                    << dCellArr[0][dartMap[0][newDartIndices[iF][iD]]] << std::endl;
        }
      }
    }

    // glue together opposite faces
    for (size_t iF = 0; iF < cell.size(); ++iF) {
      const std::vector<size_t> face = cell[iF];
      CanonicalVertexList canon = canonicalize(face);
      size_t twinFaceDart = oppDartLookup(face, canon);
      if (twinFaceDart == INVALID_IND) { // if the opposite face does not exist yet, set dartMap[2] to INVALID_IND
        for (size_t iD = 0; iD < face.size(); iD++) dartMap[2][newDartIndices[iF][iD]] = INVALID_IND;
        createdDarts[canon.vertices] = std::make_tuple(newDartIndices[iF][0], canon.rotation, canon.orientation);
      } else { // if the opposite face has already created, hook up the appropriate pointers
        for (size_t i = face.size(); i > 0; i--) {
          attachDartMap2(newDartIndices[iF][i % face.size()], twinFaceDart);
          twinFaceDart = dartMap[0][twinFaceDart];
        }
        size_t myDart = newDartIndices[iF][0], oppDart = twinFaceDart;
        GC_SAFETY_ASSERT(dCellArr[0][myDart] == dCellArr[0][dartMap[0][oppDart]], "dart gluing misaligned");
      }
    }
  }

  nCellsCapacityCount[0] = nCellsCount[0];
  nCellsFillCount[0] = nCellsCount[0];
  nDartsCapacityCount = nDartsCount;
  nDartsFillCount = nDartsCount;

  // construct 1-cells and 2-cells
  indexCells(1);
  indexCells(2);
}
#endif

// template <size_t D> // Construct a cell complex given as an array of boundary matrices
// CombinatorialMap<D>::CombinatorialMap(const std::array<SparseMatrix<int>, D>& boundaryMatrices) {
//   const bool DEBUG_PRINT = true;

//   // reserve space for k-cells for each dimension 0 <= k <= D
//   for (size_t dim = 0; dim < D; dim++) {
//     const SparseMatrix<int>& B = boundaryMatrices[dim]; // boundary_{dim+1} matrix
//     allocateCells(dim, B.cols());
//     if (dim + 1 == D) { // for top matrix, read off top-dimensional cell count too
//       allocaateCells(D, B.rows());
//     } else { // make sure that matrix products B_dim . B_{dim+1} make sense when B_{dim+1} exists
//       GC_SAFETY_ASSERT(B.rows() == boundaryMatrices[dim + 1].cols(), "boundary matrix dimensions do not match");
//     }
//   }
// }

// TODO: construct d-dimensional version properly
// // Builds a 2D polygon mesh
// template <>
// CombinatorialMap<2>::CombinatorialMap(const std::vector<std::vector<size_t>>& polygons) {
//   std::map<std::pair<size_t, size_t>, size_t> edgeIndices;
//   size_t nEdges = 0;

//   std::array<std::vector<std::vector<std::pair<size_t, bool>>>, 2> boundaryMaps;

//   for (const std::vector<size_t>& face : polygons) {
//     boundaryMaps[1].push_back(std::vector<std::pair<size_t, bool>>{});
//     for (size_t iE = 0; iE < face.size(); iE++) {
//       size_t vi = face[iE], vj = face[(iE + 1) % face.size()];
//       std::pair<size_t, size_t> key = std::minmax(vi, vj);
//       bool orientation = vi < vj;
//       if (edgeIndices.find(key) == edgeIndices.end()) { // first time seeing this edge, add to boundary map
//         edgeIndices[key] = nEdges;
//         boundaryMaps[0].push_back(
//             std::vector<std::pair<size_t, bool>>{std::make_pair(vj, orientation), std::make_pair(vi,
//             !orientation)});
//         nEdges++;
//       }
//       boundaryMaps[1].back().push_back(std::make_pair(edgeIndices[key], orientation));
//     }
//   }

//   constructFromBoundaryMaps(boundaryMaps);
// }

template <size_t D> // Construct a cell complex given as an array of boundary matrices
// CombinatorialMap<D>::CombinatorialMap(const std::array<SparseMatrix<int>, D>& boundaryMatrices) {
CombinatorialMap<D>::CombinatorialMap(
    const std::array<std::vector<std::vector<std::pair<size_t, bool>>>, D>& boundaryMaps) {
  constructFromBoundaryMaps(boundaryMaps);
}

template <size_t D> // Construct a cell complex given as an array of boundary matrices
                    // CombinatorialMap<D>::CombinatorialMap(const std::array<SparseMatrix<int>, D>& boundaryMatrices)
                    // {
void CombinatorialMap<D>::constructFromBoundaryMaps(
    const std::array<std::vector<std::vector<std::pair<size_t, bool>>>, D>& boundaryMaps) {
  const bool DEBUG_PRINT = false;

  if (DEBUG_PRINT) {
    std::cout << "Constructing combinatorial map from boundary maps" << std::endl;
    for (size_t k = 0; k < boundaryMaps.size(); k++) {
      std::cout << "Boundary map " << k << " ===========" << std::endl;
      for (size_t iC = 0; iC < boundaryMaps[k].size(); iC++) {
        std::cout << "  Cell " << iC << ":" << std::endl << "      ";
        for (size_t iB = 0; iB < boundaryMaps[k][iC].size(); iB++) {
          std::cout << " " << (boundaryMaps[k][iC][iB].second ? "+" : "-") << boundaryMaps[k][iC][iB].first;
        }
        std::cout << std::endl;
      }
    }
  }

  // Sanity check that each edge has length 2, with one positive vertex and one negative vertex
  for (size_t iE = 0; iE < boundaryMaps[0].size(); iE++) {
    GC_SAFETY_ASSERT(boundaryMaps[0][iE].size() == 2, "edge " + std::to_string(iE) + " has " +
                                                          std::to_string(boundaryMaps[0][iE].size()) +
                                                          " vertices instead of the expected 2");
    GC_SAFETY_ASSERT(boundaryMaps[0][iE][0].second != boundaryMaps[0][iE][1].second,
                     "edge " + std::to_string(iE) + " has two vertices with the same sign");
  }

  // count 0-cells
  size_t nVertices = 0;
  for (const std::vector<std::pair<size_t, bool>>& edge : boundaryMaps[0]) {
    for (const std::pair<size_t, bool>& vertex : edge) {
      nVertices = std::max(nVertices, vertex.first);
    }
  }
  nVertices++; // 0-based means count is max + 1
  allocateCells(0, nVertices);

  // reserve space for k-cells for each dimension 1 <= k <= D
  for (size_t dim = 0; dim < D; dim++) allocateCells(dim + 1, boundaryMaps[dim].size());

  // allocateDarts(nCellDarts);

  // Represent darts as oriented flags. The dart {(i, oi), (j, oj), (k, ok), ...} corresponds to the flag consisting
  // of edge i with orientation oi, face j with orientation oj, 3-cell k with orientation ok, etc. Since we store the
  // orientation explicitly, we do not store the 0-dimensional vertex of the flag, which is uniquely determined by the
  // edge + orientation cellDarts[k][i] lists the darts making up k-cell i.
  std::array<std::vector<std::vector<std::pair<size_t, bool>>>, D + 1> cellDarts;

  // fill in edge darts directly
  for (size_t iE = 0; iE < boundaryMaps[0].size(); iE++) {
    return cellDarts[1].push_back({std::vector<std::pair<size_t, bool>>{std::make_pair(iE, true)}});
  }

  throw std::runtime_error("constructor from boundary maps not implemented yet");

  // // recursively construct darts for higher-dimensional cells
  // for (size_t dim = 2; dim <= D; dim++) {
  //   for (size_t iC = 0; iC < boundaryMaps[dim - 1].size(); iC++) {
  //     cellDarts[dim].push_back(std::vector<std::vector<std::pair<size_t, bool>>>{});
  //     for (std::pair<size_t, bool> bdyFace : boundaryMaps[dim - 1][iC]) {
  //     }
  //   }
  // }


  // std::vector<std::vector<std::pair<size_t, bool>>> darts;

  // for (size_t iCell = 0; iCell < boundaryMaps[D - 1].size(); iCell++) {
  // }


  // // TODO: keep indices from input mesh
  // for (size_t k = 0; k <= D; k++) indexCells(k);
}

template <size_t D>
void CombinatorialMap<D>::validateConnectivity(bool allowDeadDarts) {
  const bool DEBUG_PRINT = false;

  // Sanity check sizes and counts
  if (nDartsCount > nDartsFillCount) throw std::logic_error("dart count > dart fill");
  if (nDartsFillCount > nDartsCapacityCount) throw std::logic_error("dart fill > dart capacity");

  for (size_t dim = 0; dim <= D; dim++) {
    if (nCellsCount[dim] > nCellsFillCount[dim])
      throw std::logic_error(std::to_string(dim) + "-cell count (" + std::to_string(nCellsCount[dim]) + ") > " +
                             std::to_string(dim) + "-cell fill (" + std::to_string(nCellsFillCount[dim]) + ")");
    if (nCellsFillCount[dim] > nCellsCapacityCount[dim])
      throw std::logic_error(std::to_string(dim) + "-cell fill (" + std::to_string(nCellsFillCount[dim]) + ") > " +
                             std::to_string(dim) + "-cell capacity (" + std::to_string(nCellsCapacityCount[dim]) + ")");
    if (dim < D && dartMap[dim].size() != nDartsCapacityCount)
      throw std::logic_error("dart map[" + std::to_string(dim) + "] has size " + std::to_string(dartMap[dim].size()) +
                             " even though nDartsCapacityCount is " + std::to_string(nDartsCapacityCount));
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
    if (iC >= nCellsFillCount[k]) {
      throw std::logic_error(msg + " - bad " + std::to_string(k) + "-cell reference " + std::to_string(iC) +
                             " is larger than nCellsFillCOunt[k]: " + std::to_string(nCellsFillCount[k]));
    } else if (cellIsDead(k, iC)) {
      throw std::logic_error(msg + " - bad " + std::to_string(k) + "-cell reference is dead");
    }
  };

  if (DEBUG_PRINT) {
    std::cout << "====== validating connectivity of: " << std::endl;
    std::cout << "#\tnext\ttwin" << std::endl;
    for (size_t iDart = 0; iDart < this->nDartsCount; iDart++) {
      std::cout << iDart << "\t" << this->dartPartner(iDart, 0) << "\t" << this->dartPartner(iDart, 1) << std::endl;
    }
  }


  // == Darts
  // Check dartMap invariants and pointer validity
  // Note: we intentionally mostly avoid using iterators here, because they can be hard to debug when things are
  // broken.
  for (size_t iDart = 0; iDart < nDartsFillCount; iDart++) {
    GC_SAFETY_ASSERT(allowDeadDarts || !dartIsDead(iDart),
                     "invalid mesh -- dead dart"); // no darts should be dead yet
    if (dartIsDead(iDart)) continue;

    // check partner, only allow INVALID_IND partner for dim = D-1
    // check that darts point opposite to each other
    size_t vTail = dCellArr[0][iDart], vTip = dCellArr[0][dartMap[0][iDart]];
    for (size_t dim = 0; dim < D; dim++) {
      validateDart(dartMap[dim][iDart], "dart " + std::to_string(iDart) + ".partner(" + std::to_string(dim) + ")",
                   dim == D - 1);
      // if (dartMap[0][iDart] == INVALID_IND && dartMap[D - 1][iDart] == INVALID_IND) {
      //   throw std::logic_error("dart " + std::to_string(iDart) + ".partner() and .partner(" + std::to_string(D - 1)
      //   +
      //                          ") are both INVALID_IND - bad dart reference");
      // }
      if (dartMap[dim][iDart] != INVALID_IND && dim > 0) {
        size_t vTailOp = dCellArr[0][dartMap[dim][iDart]], vTipOp = dCellArr[0][dartMap[0][dartMap[dim][iDart]]];
        if ((vTailOp != vTip || vTipOp != vTail)) {
          throw std::logic_error("dart " + std::to_string(iDart) + " : " + std::to_string(vTail) + "->" +
                                 std::to_string(vTip) + " is glued to dart " + std::to_string(dartMap[dim][iDart]) +
                                 " : " + std::to_string(vTailOp) + "->" + std::to_string(vTipOp) + " by dart map " +
                                 std::to_string(dim) + " but the two are misaligned");
        }
      }
    }

    // == Check that dartMap[1] .. dartMap[D-1] are all involutions and commute properly with dartMap[0]
    // Note that we already verified that dartMap[k] != INVALID_IND for k < D-1
    for (size_t k = 1; k < D; k++) {
      size_t doublePartner = dartPartner(dartPartner(iDart, k), k);
      if ((k < D - 1 && doublePartner != iDart) ||
          (k == D - 1 && !(doublePartner == iDart || doublePartner == INVALID_IND))) {
        throw std::logic_error("dartMap[" + std::to_string(k) + "] is not an involution. Applying twice to dart " +
                               std::to_string(iDart) + " yields dart " + std::to_string(doublePartner));
      }
      // we require that .next().twin().next() = id for k > 1
      if (k > 1) {
        size_t nextTwinNextTwin = dartPartner(dartPartner(dartPartner(dartPartner(iDart, 0), k), 0), k);
        if ((k < D - 1 && nextTwinNextTwin != iDart) ||
            (k == D - 1 && !(nextTwinNextTwin == iDart || nextTwinNextTwin == INVALID_IND))) {
          std::cout << "dartMap[0](dart " << iDart << ") = dart " << dartPartner(iDart, 0) << std::endl;
          std::cout << "dartMap[" << k << "](..) = dart " << dartPartner(dartPartner(iDart, 0), k) << std::endl;
          std::cout << "dartMap[0](..) = dart " << dartPartner(dartPartner(dartPartner(iDart, 0), k), 0) << std::endl;
          std::cout << "dartMap[" << k << "](..) = dart " << nextTwinNextTwin << std::endl;
          throw std::logic_error("dartMap[0] and dartMap[" + std::to_string(k) +
                                 "] do not satisfy commutation relation at dart " + std::to_string(iDart));
        }
      }
    }

    // validate dart.cell<k>()
    for (size_t k = 0; k <= D; k++) validateCell(k, dCellArr[k][iDart], "dart.cell<" + std::to_string(k) + ">()");
  }

  for (size_t k = 0; k <= D; k++) {
    for (size_t iC = 0; iC < nCellsFillCount[k]; iC++) {
      if (cellIsDead(k, iC)) continue;
      validateDart(cDartArr[k][iC], "cell.dart()", false);
    }
  }

  // check adjacency sanity
  for (size_t k = 0; k <= D; k++) {                      // check k-cells
    for (size_t iC = 0; iC < nCellsFillCount[k]; iC++) { // check that cell = cell.dart().cell()
      if (cellIsDead(k, iC)) continue;
      if (dCellArr[k][cDartArr[k][iC]] != iC) {
        std::cout << k << "-cell " << iC << ", dart " << cDartArr[k][iC] << " | dart.cell<" << k
                  << "() = " << dCellArr[k][cDartArr[k][iC]] << std::endl;
        throw std::logic_error(std::to_string(k) + "-cell doesn't match dart.cell<" + std::to_string(k) + ">()");
      }
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

template <size_t D>
struct OrbitNeighborhoodIterator {
  Dart<D> startDart;
  size_t iMap;
  size_t k;

  // for k=0
  size_t jMap;
  Dart<D> iPartner;

  OrbitNeighborhoodIterator(Dart<D> d, size_t k_) : startDart(d), iMap(0), k(k_), jMap(0) {
    if (k == 0) {
      iMap = 1;
      iPartner = startDart.partner(iMap);
    }
    while (!isValid() && !finished()) advance();
  }

  void advance() {
    if (k > 0) {
      iMap++;
    } else {
      jMap++;
      if (jMap >= iMap) {
        jMap = 0;
        iMap++;
        if (iMap < D) iPartner = startDart.partner(iMap);
      }
    }
  }

  bool isValid() const { return (k == 0) || (iMap + 1 != k); }

  bool finished() const { return iMap >= D; }

  const OrbitNeighborhoodIterator<D>& operator++() {
    do {
      advance();
    } while (!isValid() && !finished());
    return *this;
  }

  // any two finished iterators are equal, otherwise compare internals
  bool operator==(const OrbitNeighborhoodIterator<D>& other) const {
    return (finished() && other.finished()) || (startDart == other.startDart && iMap == other.iMap && k == other.k);
  }

  bool operator!=(const OrbitNeighborhoodIterator<D>& other) const { return !(*this == other); }

  std::pair<Dart<D>, bool> operator*() const {
    if (k > 0) {
      bool orientationPreserving = (iMap + 1) < k;

      Dart<D> currE;
      if (iMap < D) {
        currE = startDart.partner(iMap);
      } else {
        currE = startDart;
      }

      return std::make_pair(currE, orientationPreserving);
    } else {
      Dart<D> currE;
      if (iPartner != startDart) { // not INVALID_IND
        currE = iPartner.partner(jMap);
      } else {
        currE = startDart;
      }
      return std::make_pair(currE, true); // orientation always preserving for k=0
    }
  }
};


template <size_t D>
class OrbitNeighborhood {
public:
  OrbitNeighborhood(Dart<D> d, size_t k_) : k(k_), dStart(d), cachedEnd(d, k) { cachedEnd.iMap = D + 1; }

  OrbitNeighborhoodIterator<D> begin() const { return OrbitNeighborhoodIterator<D>(dStart, k); }

  // since cachedEnd.finished() == true, checking equality with cachedEnd checks if an iterator is finished()
  OrbitNeighborhoodIterator<D> end() const { return cachedEnd; }

private:
  size_t k;
  Dart<D> dStart;
  OrbitNeighborhoodIterator<D> cachedEnd;
};


// k-cells are orbits generated by compositions of dart maps
// 0-cells are generated by <map[0].map[1], map[0].map[2], ..., map[D-2].map[D-1]>
// 1-cells are generated by <map[1], ..., map[D-1]>
// 2-cells are generated by <map[0], map[2], ..., map[D-1]>
// see e.g. https://doc.cgal.org/latest/Combinatorial_map/index.html#title3
// so incidences are intersections of these orbits. Use union-find to index these orbit intersections

template <size_t D>
std::vector<Dart<D>> incidenceNeighboringDarts(Dart<D> d, size_t k1, size_t k2, bool verbose) {
  if (verbose) std::cout << "... computing neighboring darts for " << d << std::endl;
  std::vector<Dart<D>> result;
  if (k1 == 0) { // k1 = 0, k2 > 0. Take 0-cell compositions, skipping k2-1
    for (size_t iMap = 1; iMap < D; ++iMap) {
      if (iMap + 1 == k2) continue;
      for (size_t jMap = 0; jMap < iMap; ++jMap) {
        if (jMap + 1 == k2) continue;
        Dart<D> iPartner = d.partner(iMap);
        if (iPartner == d) continue; // skip (INVALID_IND)
        Dart<D> next = iPartner.partner(jMap);
        if (next == iPartner) continue; // skip (INVALID_IND)
        result.push_back(next);
      }
    }
  } else { // k1 > 0, k2 > 0: skip maps k1 - 1, k2 - 1
    for (size_t iMap = 0; iMap < D; ++iMap) {
      if (iMap + 1 != k1 && iMap + 1 != k2) { // iMap != k-1
        Dart<D> next = d.partner(iMap);
        if (next == d) continue; // skip (INVALID_IND)
        result.push_back(next);
      }
    }
  }
  return result;
}

namespace unionfind {
// find root, and update all nodes in path to point to root, updating their `sharesParentSign` fields as necessary
inline size_t findRoot(size_t x, std::vector<size_t>& parent, std::vector<bool>& sharesParentSign) {
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
}

// join together x and y, updating their `sharesParentSign` fields as necessary
inline void unite(size_t x, size_t y, bool samesign, std::vector<size_t>& parent, std::vector<bool>& sharesParentSign,
                  std::vector<size_t>& rank) {
  size_t rootX = findRoot(x, parent, sharesParentSign), rootY = findRoot(y, parent, sharesParentSign);
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
}

} // namespace unionfind

} // namespace combinatorial_map
} // namespace geometrycentral
