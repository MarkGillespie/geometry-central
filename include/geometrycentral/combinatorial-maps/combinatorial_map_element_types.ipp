// Implementations for combinatorial_map_mesh_types.h

namespace std {
template <size_t D>
struct hash<geometrycentral::combinatorial_map::Dart<D>> {
  std::size_t operator()(const geometrycentral::combinatorial_map::Dart<D>& e) const {
    return std::hash<size_t>{}(e.getIndex());
  }
};

template <size_t k, size_t D>
struct hash<geometrycentral::combinatorial_map::Cell<k, D>> {
  std::size_t operator()(const geometrycentral::combinatorial_map::Cell<k, D>& c) const {
    return std::hash<size_t>{}(c.getIndex());
  }
};

template <size_t k1, size_t k2, size_t D>
struct hash<geometrycentral::combinatorial_map::Incidence<k1, k2, D>> {
  std::size_t operator()(const geometrycentral::combinatorial_map::Incidence<k1, k2, D>& e) const {
    return std::hash<size_t>{}(e.getIndex());
  }
};
} // namespace std

namespace geometrycentral {
namespace combinatorial_map {

// ==========================================================
// ================        Dart        ==================
// ==========================================================

// Constructors
template <size_t D>
inline Dart<D>::Dart() {}

template <size_t D>
inline Dart<D>::Dart(CombinatorialMap<D>* mesh_, size_t ind_) : Element<Dart<D>, CombinatorialMap<D>>(mesh_, ind_) {}

// Navigators

template <size_t D>
template <size_t iPartner>
inline Dart<D> Dart<D>::partner() const {
  size_t partnerInd = this->mesh->dartPartner(this->ind, iPartner);
  return (partnerInd == INVALID_IND) ? *this : Dart<D>(this->mesh, partnerInd);
}

template <size_t D>
inline Dart<D> Dart<D>::partner(size_t d) const {
  size_t partnerInd = this->mesh->dartPartner(this->ind, d);
  return (partnerInd == INVALID_IND) ? *this : Dart<D>(this->mesh, partnerInd);
}

template <size_t D>
inline Dart<D> Dart<D>::next() const {
  return partner(0);
}

template <size_t D>
template <size_t k>
inline Cell<k, D> Dart<D>::cell() const { // k-cell
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return Cell<k, D>(this->mesh, this->mesh->dCellArr[k][this->ind], this->mesh->dCellSgn[k][this->ind]);
}

template <size_t D>
inline Vertex<D> Dart<D>::vertex() const {
  return cell<0>();
}

template <size_t D>
inline Edge<D> Dart<D>::edge() const {
  return cell<1>();
}

template <size_t D>
inline Face<D> Dart<D>::face() const {
  return cell<2>();
}

template <size_t D>
inline Cell<3, D> Dart<D>::cell() const { // 3-cell
  return cell<3>();
}

template <size_t D>
inline Vertex<D> Dart<D>::tailVertex() const {
  return vertex();
}
template <size_t D>
inline Vertex<D> Dart<D>::tipVertex() const {
  return next().vertex();
}

template <size_t D>
template <size_t k1, size_t k2>
inline Incidence<k1, k2, D> Dart<D>::incidence() const {
  static_assert(k1 < k2, "an incidence must have cell dimensions k1 < k2");
  static_assert(k2 <= D, "cell dimension k2 must be less than or equal to complex dimension D");
  this->mesh->ensureHaveIncidences(k1, k2);
  return Incidence<k1, k2, D>(this->mesh, this->mesh->dIncidenceArr[std::make_pair(k1, k2)][this->ind]);
}

//== Aliases for some common incidences
template <size_t D>
inline Incidence<0, D, D> Dart<D>::vertexCorner() const {
  return incidence<0, D>();
}
template <size_t D>
inline Incidence<1, D, D> Dart<D>::edgeCorner() const {
  return incidence<1, D>();
}
template <size_t D>
inline Incidence<0, 2, D> Dart<D>::faceCorner() const {
  return incidence<0, 2>();
}

template <size_t D>
inline bool Dart<D>::isDead() const {
  return this->mesh->dartIsDead(this->ind);
}

// Range iterators
template <size_t D>
inline bool DartRangeF<D>::elementOkay(const CombinatorialMap<D>& mesh, size_t ind) {
  return !mesh.dartIsDead(ind);
}

// ==========================================================
// ================        k-Cell        ==================
// ==========================================================

// Constructors
template <size_t k, size_t D>
inline Cell<k, D>::Cell() {}

template <size_t k, size_t D>
inline Cell<k, D>::Cell(CombinatorialMap<D>* mesh_, size_t ind_)
    : Element<Cell<k, D>, CombinatorialMap<D>>(mesh_, ind_) {}

template <size_t k, size_t D>
inline Cell<k, D>::Cell(CombinatorialMap<D>* mesh_, size_t ind_, bool orientation)
    : Element<Cell<k, D>, CombinatorialMap<D>>(mesh_, ind_), mOrientation(orientation) {}

// Navigators
template <size_t k, size_t D>
inline Dart<D> Cell<k, D>::dart() const {
  return Dart<D>(this->mesh, this->mesh->cDartArr[k][this->ind]);
}

template <size_t k, size_t D>
template <size_t k2>
inline Dart<D> Cell<k, D>::dartInCell(Cell<k2, D> cell) const {
  return this->mesh->adjacentDartInCell(*this, cell);
}

template <size_t k1, size_t D>
inline std::vector<Dart<D>> Cell<k1, D>::adjacentDarts() const {
  return this->mesh->adjacentDarts(*this);
}

template <size_t k1, size_t D>
template <size_t k2>
inline std::vector<Cell<k2, D>> Cell<k1, D>::adjacentCells() const {
  return this->mesh->template adjacentCells<k1, k2>(*this);
}

template <size_t k, size_t D>
inline std::vector<Vertex<D>> Cell<k, D>::adjacentVertices() const {
  return this->mesh->template adjacentVertices<k>(*this);
}

template <size_t k, size_t D>
inline std::vector<Edge<D>> Cell<k, D>::adjacentEdges() const {
  return this->mesh->template adjacentEdges<k>(*this);
}

template <size_t k, size_t D>
inline std::vector<Face<D>> Cell<k, D>::adjacentFaces() const {
  return this->mesh->template adjacentFaces<k>(*this);
}

template <size_t k, size_t D>
inline std::vector<Cell<3, D>> Cell<k, D>::adjacentCells() const { // adjacentCells gives 3-cells
  return this->mesh->template adjacentCells<k, 3>(*this);
}

// function for iterating over adjacent incidences with higher-dimensional cells
// OrderedIncidence<a, b> is just an ordinary incidence, but with a and b ordered properly, i.e. Incidence<a,b> if a < b
// and Incidence<b,a> otherwise
template <size_t k, size_t D>
template <size_t k2>
std::vector<OrderedIncidence<k, k2, D>> Cell<k, D>::adjacentIncidences() const {
  return this->mesh->template adjacentIncidences<k, k2>(*this);
}

template <size_t k, size_t D>
std::vector<Incidence<0, D, D>> Cell<k, D>::adjacentVertexCorners() const {
  static_assert(k == 0 || k == D,
                "adjacentVertexCorners() is only defined for 0-cells and D-cells on a D-dimensional cell complex");
  return adjacentIncidences<k == 0 ? D : 0>();
}

template <size_t k, size_t D>
std::vector<Incidence<1, D, D>> Cell<k, D>::adjacentEdgeCorners() const {
  static_assert(k == 1 || k == D,
                "adjacentEdgeCorners() is only defined for 1-cells and D-cells on a D-dimensional cell complex");
  return adjacentIncidences<k == 1 ? D : 1>();
}

template <size_t k, size_t D>
std::vector<Incidence<0, 2, D>> Cell<k, D>::adjacentFaceCorners() const {
  static_assert(k == 0 || k == 2,
                "adjacentFaceCorners() is only defined for 0-cells and 2-cells on a D-dimensional cell complex");
  return adjacentIncidences<k == 0 ? 2 : 0>();
}

template <size_t k, size_t D>
inline bool Cell<k, D>::isDead() const {
  return this->mesh->cellIsDead<k>(this->ind);
}

template <size_t k, size_t D>
inline bool Cell<k, D>::isBoundary() const {
  if (k == D - 1) {
    return dart().partner(D - 1) == dart();
  } else {
    for (Cell<D - 1, D> facet : adjacentCells<D - 1>()) {
      if (facet.isBoundary()) return true;
    }
    return false;
  }
}

template <size_t k, size_t D>
inline bool Cell<k, D>::orientation() const {
  return mOrientation;
}

template <size_t k, size_t D>
inline void Cell<k, D>::flipOrientation() {
  mOrientation = !mOrientation;
}

template <size_t k, size_t D>
inline void Cell<k, D>::setOrientation(bool orientation) {
  mOrientation = orientation;
}

template <size_t k, size_t D>
inline int Cell<k, D>::sign() const {
  return mOrientation ? 1 : -1;
}

template <size_t k, size_t D>
inline bool Cell<k, D>::orientationInCell(Cell<k + 1, D> c) const {
  static_assert(k + 1 <= D, "cannot construct a (D+1)-cell");
  // get the orientation of shared dart, and then flip if cells are oppositely oriented
  // (note that == on bools is XOR)
  size_t iDart = this->mesh->adjacentDartInCell(*this, c).getIndex();
  return (this->mesh->dCellSgn[k + 1][iDart] == this->mesh->dCellSgn[k][iDart]) == (c.orientation() == orientation());
}

template <size_t k, size_t D>
inline bool Cell<k, D>::orientationInCell(Cell<k - 1, D> c) const {
  static_assert(k > 0, "cannot construct a (-1)-cell");
  return c.orientationInCell(*this);
}

template <size_t k, size_t D>
inline int Cell<k, D>::signInCell(Cell<k + 1, D> c) const {
  return orientationInCell(c) ? 1 : -1;
}

template <size_t k, size_t D>
inline int Cell<k, D>::signInCell(Cell<k - 1, D> c) const {
  return orientationInCell(c) ? 1 : -1;
}

// Range iterators
template <size_t k, size_t D>
inline bool CellRangeF<k, D>::elementOkay(const CombinatorialMap<D>& mesh, size_t ind) {
  return !(mesh.template cellIsDead<k>(ind));
}

// ==========================================================
// ==============       k-Cell Incidence       ==============
// ==========================================================

// Constructors
template <size_t k1, size_t k2, size_t D>
inline Incidence<k1, k2, D>::Incidence() {}

template <size_t k1, size_t k2, size_t D>
inline Incidence<k1, k2, D>::Incidence(CombinatorialMap<D>* mesh_, size_t ind_)
    : Element<Incidence<k1, k2, D>, CombinatorialMap<D>>(mesh_, ind_) {}

// Navigators
template <size_t k1, size_t k2, size_t D>
inline Dart<D> Incidence<k1, k2, D>::dart() const {
  return Dart<D>(this->mesh, this->mesh->iDartArr[std::make_pair(k1, k2)][this->ind]);
}

template <size_t k1, size_t k2, size_t D>
template <size_t k>
inline Cell<k, D> Incidence<k1, k2, D>::cell() const {
  static_assert(k == k1 || k == k2, "Incidence<k1, k2>::cell<k>() is only defined for k = k1 or k = k2");
  return dart().template cell<k>();
}

template <size_t k1, size_t k2, size_t D>
inline Vertex<D> Incidence<k1, k2, D>::vertex() const {
  return cell<0>();
}
template <size_t k1, size_t k2, size_t D>
inline Edge<D> Incidence<k1, k2, D>::edge() const {
  return cell<1>();
}
template <size_t k1, size_t k2, size_t D>
inline Face<D> Incidence<k1, k2, D>::face() const {
  return cell<2>();
}
template <size_t k1, size_t k2, size_t D>
inline Cell<3, D> Incidence<k1, k2, D>::cell() const {
  return cell<3>();
}

template <size_t k1, size_t k2, size_t D>
inline std::vector<Dart<D>> Incidence<k1, k2, D>::adjacentDarts() const {
  return this->mesh->adjacentDarts(*this);
}

template <size_t k1, size_t k2, size_t D>
template <size_t k>
inline std::vector<Cell<k, D>> Incidence<k1, k2, D>::adjacentCells() const {
  return this->mesh->template adjacentCells<k1, k2, k>(*this);
}

template <size_t k1, size_t k2, size_t D>
inline std::vector<Vertex<D>> Incidence<k1, k2, D>::adjacentVertices() const {
  return adjacentCells<0>();
}
template <size_t k1, size_t k2, size_t D>
inline std::vector<Edge<D>> Incidence<k1, k2, D>::adjacentEdges() const {
  return adjacentCells<1>();
}
template <size_t k1, size_t k2, size_t D>
inline std::vector<Face<D>> Incidence<k1, k2, D>::adjacentFaces() const {
  return adjacentCells<2>();
}
template <size_t k1, size_t k2, size_t D>
inline std::vector<Cell<3, D>> Incidence<k1, k2, D>::adjacentCells() const {
  return adjacentCells<3>();
}

template <size_t k1, size_t k2, size_t D>
inline bool Incidence<k1, k2, D>::isDead() const {
  return this->mesh->incidenceIsDead<k1, k2>(this->ind);
}

template <size_t k1, size_t k2, size_t D>
inline bool Incidence<k1, k2, D>::isBoundary() const {
  if (k2 == D - 1) {
    return dart().partner(D - 1) == dart();
  } else {
    for (Cell<D - 1, D> facet : adjacentCells<D - 1>()) {
      if (facet.isBoundary()) return true;
    }
    return false;
  }
}

// == Range iterators

// All vertices
template <size_t k1, size_t k2, size_t D>
inline bool IncidenceRangeF<k1, k2, D>::elementOkay(const CombinatorialMap<D>& mesh, size_t ind) {
  return !(mesh.template incidenceIsDead<k1, k2>(ind));
}

} // namespace combinatorial_map

// ==========================================================
// === Customization points for MeshData<Cell<k,D>, T>, MeshData<Dart<D>, T>,
// === and MeshData<Incidence<k1,k2,D>, T>.
//
// element.h declares nElements<E>, elementCapacity<E>, dataIndexOfElement<E>,
// ElementSetType<E>, iterateElements<E>, getExpandCallbackList<E>, and
// getPermuteCallbackList<E> as the customization points MeshData<E,T> needs;
// every element type must fully specialize all of them for every concrete
// type it wants MeshData to support. element.h is shared with every other
// module in geometry-central and is deliberately left untouched here, and
// C++ does not allow *partial* specialization of a function template, so
// there is no way to write "one specialization covers Cell<k, D> for every
// k and D" directly from this file either -- the three macros below
// generate one full specialization set per concrete (k, D) / D / (k1, k2,
// D) instead (so a human only has to write and read one short invocation
// per case, not a ~10-line hand-written body), and are invoked below for
// every combination up to the bounds noted at each list. Extend a bound by
// adding more invocation lines in the same pattern if a higher-dimensional
// map needs it.
// ==========================================================

#define GC_CM_SPECIALIZE_CELL(k, D) \
  template <> \
  struct ElementSetType<combinatorial_map::Cell<k, D>> { \
    typedef combinatorial_map::CellSet<k, D> type; \
  }; \
  template <> \
  inline size_t nElements<combinatorial_map::Cell<k, D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->template nCells<k>(); \
  } \
  template <> \
  inline size_t elementCapacity<combinatorial_map::Cell<k, D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->template nCellsCapacity<k>(); \
  } \
  template <> \
  inline size_t dataIndexOfElement<combinatorial_map::Cell<k, D>>(combinatorial_map::CombinatorialMap<D>* mesh, \
                                                                  combinatorial_map::Cell<k, D> e) { \
    (void)mesh; \
    return e.getIndex(); \
  } \
  template <> \
  inline combinatorial_map::CellSet<k, D> iterateElements<combinatorial_map::Cell<k, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->template cells<k>(); \
  } \
  template <> \
  inline std::list<std::function<void(size_t)>>& getExpandCallbackList<combinatorial_map::Cell<k, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->cellExpandCallbackList[k]; \
  } \
  template <> \
  inline std::list<std::function<void(const std::vector<size_t>&)>>& \
  getPermuteCallbackList<combinatorial_map::Cell<k, D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->cellPermuteCallbackList[k]; \
  }

#define GC_CM_SPECIALIZE_DART(D) \
  template <> \
  struct ElementSetType<combinatorial_map::Dart<D>> { \
    typedef combinatorial_map::DartSet<D> type; \
  }; \
  template <> \
  inline size_t nElements<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->nDarts(); \
  } \
  template <> \
  inline size_t elementCapacity<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->nDartsCapacity(); \
  } \
  template <> \
  inline size_t dataIndexOfElement<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh, \
                                                               combinatorial_map::Dart<D> e) { \
    (void)mesh; \
    return e.getIndex(); \
  } \
  template <> \
  inline combinatorial_map::DartSet<D> iterateElements<combinatorial_map::Dart<D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->darts(); \
  } \
  template <> \
  inline std::list<std::function<void(size_t)>>& getExpandCallbackList<combinatorial_map::Dart<D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->dartExpandCallbackList; \
  } \
  template <> \
  inline std::list<std::function<void(const std::vector<size_t>&)>>& \
  getPermuteCallbackList<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    return mesh->dartPermuteCallbackList; \
  }

#define GC_CM_SPECIALIZE_INCIDENCE(k1, k2, D) \
  template <> \
  struct ElementSetType<combinatorial_map::Incidence<k1, k2, D>> { \
    typedef combinatorial_map::IncidenceSet<k1, k2, D> type; \
  }; \
  template <> \
  inline size_t nElements<combinatorial_map::Incidence<k1, k2, D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    mesh->ensureHaveIncidences(k1, k2); \
    return mesh->template nIncidences<k1, k2>(); \
  } \
  template <> \
  inline size_t elementCapacity<combinatorial_map::Incidence<k1, k2, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    mesh->ensureHaveIncidences(k1, k2); \
    return mesh->template nIncidencesCapacity<k1, k2>(); \
  } \
  template <> \
  inline size_t dataIndexOfElement<combinatorial_map::Incidence<k1, k2, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh, combinatorial_map::Incidence<k1, k2, D> e) { \
    (void)mesh; \
    return e.getIndex(); \
  } \
  template <> \
  inline combinatorial_map::IncidenceSet<k1, k2, D> iterateElements<combinatorial_map::Incidence<k1, k2, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    mesh->ensureHaveIncidences(k1, k2); \
    return mesh->template incidences<k1, k2>(); \
  } \
  template <> \
  inline std::list<std::function<void(size_t)>>& getExpandCallbackList<combinatorial_map::Incidence<k1, k2, D>>( \
      combinatorial_map::CombinatorialMap<D>* mesh) { \
    mesh->ensureHaveIncidences(k1, k2); \
    return mesh->incidenceExpandCallbackList[std::make_pair((size_t)k1, (size_t)k2)]; \
  } \
  template <> \
  inline std::list<std::function<void(const std::vector<size_t>&)>>& \
  getPermuteCallbackList<combinatorial_map::Incidence<k1, k2, D>>(combinatorial_map::CombinatorialMap<D>* mesh) { \
    mesh->ensureHaveIncidences(k1, k2); \
    return mesh->incidencePermuteCallbackList[std::make_pair((size_t)k1, (size_t)k2)]; \
  }

// --- Cell<k, D> for 0 <= k <= D <= 10 ---
GC_CM_SPECIALIZE_CELL(0, 0)
GC_CM_SPECIALIZE_CELL(0, 1) GC_CM_SPECIALIZE_CELL(1, 1)
GC_CM_SPECIALIZE_CELL(0, 2) GC_CM_SPECIALIZE_CELL(1, 2) GC_CM_SPECIALIZE_CELL(2, 2)
GC_CM_SPECIALIZE_CELL(0, 3) GC_CM_SPECIALIZE_CELL(1, 3) GC_CM_SPECIALIZE_CELL(2, 3) GC_CM_SPECIALIZE_CELL(3, 3)
GC_CM_SPECIALIZE_CELL(0, 4) GC_CM_SPECIALIZE_CELL(1, 4) GC_CM_SPECIALIZE_CELL(2, 4) GC_CM_SPECIALIZE_CELL(3, 4) GC_CM_SPECIALIZE_CELL(4, 4)
GC_CM_SPECIALIZE_CELL(0, 5) GC_CM_SPECIALIZE_CELL(1, 5) GC_CM_SPECIALIZE_CELL(2, 5) GC_CM_SPECIALIZE_CELL(3, 5) GC_CM_SPECIALIZE_CELL(4, 5) GC_CM_SPECIALIZE_CELL(5, 5)
GC_CM_SPECIALIZE_CELL(0, 6) GC_CM_SPECIALIZE_CELL(1, 6) GC_CM_SPECIALIZE_CELL(2, 6) GC_CM_SPECIALIZE_CELL(3, 6) GC_CM_SPECIALIZE_CELL(4, 6) GC_CM_SPECIALIZE_CELL(5, 6) GC_CM_SPECIALIZE_CELL(6, 6)
GC_CM_SPECIALIZE_CELL(0, 7) GC_CM_SPECIALIZE_CELL(1, 7) GC_CM_SPECIALIZE_CELL(2, 7) GC_CM_SPECIALIZE_CELL(3, 7) GC_CM_SPECIALIZE_CELL(4, 7) GC_CM_SPECIALIZE_CELL(5, 7) GC_CM_SPECIALIZE_CELL(6, 7) GC_CM_SPECIALIZE_CELL(7, 7)
GC_CM_SPECIALIZE_CELL(0, 8) GC_CM_SPECIALIZE_CELL(1, 8) GC_CM_SPECIALIZE_CELL(2, 8) GC_CM_SPECIALIZE_CELL(3, 8) GC_CM_SPECIALIZE_CELL(4, 8) GC_CM_SPECIALIZE_CELL(5, 8) GC_CM_SPECIALIZE_CELL(6, 8) GC_CM_SPECIALIZE_CELL(7, 8) GC_CM_SPECIALIZE_CELL(8, 8)
GC_CM_SPECIALIZE_CELL(0, 9) GC_CM_SPECIALIZE_CELL(1, 9) GC_CM_SPECIALIZE_CELL(2, 9) GC_CM_SPECIALIZE_CELL(3, 9) GC_CM_SPECIALIZE_CELL(4, 9) GC_CM_SPECIALIZE_CELL(5, 9) GC_CM_SPECIALIZE_CELL(6, 9) GC_CM_SPECIALIZE_CELL(7, 9) GC_CM_SPECIALIZE_CELL(8, 9) GC_CM_SPECIALIZE_CELL(9, 9)
GC_CM_SPECIALIZE_CELL(0, 10) GC_CM_SPECIALIZE_CELL(1, 10) GC_CM_SPECIALIZE_CELL(2, 10) GC_CM_SPECIALIZE_CELL(3, 10) GC_CM_SPECIALIZE_CELL(4, 10) GC_CM_SPECIALIZE_CELL(5, 10) GC_CM_SPECIALIZE_CELL(6, 10) GC_CM_SPECIALIZE_CELL(7, 10) GC_CM_SPECIALIZE_CELL(8, 10) GC_CM_SPECIALIZE_CELL(9, 10) GC_CM_SPECIALIZE_CELL(10, 10)

// --- Dart<D> for 0 <= D <= 10 ---
GC_CM_SPECIALIZE_DART(0) GC_CM_SPECIALIZE_DART(1) GC_CM_SPECIALIZE_DART(2) GC_CM_SPECIALIZE_DART(3) GC_CM_SPECIALIZE_DART(4) GC_CM_SPECIALIZE_DART(5) GC_CM_SPECIALIZE_DART(6) GC_CM_SPECIALIZE_DART(7) GC_CM_SPECIALIZE_DART(8) GC_CM_SPECIALIZE_DART(9) GC_CM_SPECIALIZE_DART(10)

// --- Incidence<k1, k2, D> for 0 <= k1 < k2 <= D <= 6 ---
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 1)
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 2)
GC_CM_SPECIALIZE_INCIDENCE(0, 2, 2) GC_CM_SPECIALIZE_INCIDENCE(1, 2, 2)
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 3)
GC_CM_SPECIALIZE_INCIDENCE(0, 2, 3) GC_CM_SPECIALIZE_INCIDENCE(1, 2, 3)
GC_CM_SPECIALIZE_INCIDENCE(0, 3, 3) GC_CM_SPECIALIZE_INCIDENCE(1, 3, 3) GC_CM_SPECIALIZE_INCIDENCE(2, 3, 3)
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 4)
GC_CM_SPECIALIZE_INCIDENCE(0, 2, 4) GC_CM_SPECIALIZE_INCIDENCE(1, 2, 4)
GC_CM_SPECIALIZE_INCIDENCE(0, 3, 4) GC_CM_SPECIALIZE_INCIDENCE(1, 3, 4) GC_CM_SPECIALIZE_INCIDENCE(2, 3, 4)
GC_CM_SPECIALIZE_INCIDENCE(0, 4, 4) GC_CM_SPECIALIZE_INCIDENCE(1, 4, 4) GC_CM_SPECIALIZE_INCIDENCE(2, 4, 4) GC_CM_SPECIALIZE_INCIDENCE(3, 4, 4)
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 5)
GC_CM_SPECIALIZE_INCIDENCE(0, 2, 5) GC_CM_SPECIALIZE_INCIDENCE(1, 2, 5)
GC_CM_SPECIALIZE_INCIDENCE(0, 3, 5) GC_CM_SPECIALIZE_INCIDENCE(1, 3, 5) GC_CM_SPECIALIZE_INCIDENCE(2, 3, 5)
GC_CM_SPECIALIZE_INCIDENCE(0, 4, 5) GC_CM_SPECIALIZE_INCIDENCE(1, 4, 5) GC_CM_SPECIALIZE_INCIDENCE(2, 4, 5) GC_CM_SPECIALIZE_INCIDENCE(3, 4, 5)
GC_CM_SPECIALIZE_INCIDENCE(0, 5, 5) GC_CM_SPECIALIZE_INCIDENCE(1, 5, 5) GC_CM_SPECIALIZE_INCIDENCE(2, 5, 5) GC_CM_SPECIALIZE_INCIDENCE(3, 5, 5) GC_CM_SPECIALIZE_INCIDENCE(4, 5, 5)
GC_CM_SPECIALIZE_INCIDENCE(0, 1, 6)
GC_CM_SPECIALIZE_INCIDENCE(0, 2, 6) GC_CM_SPECIALIZE_INCIDENCE(1, 2, 6)
GC_CM_SPECIALIZE_INCIDENCE(0, 3, 6) GC_CM_SPECIALIZE_INCIDENCE(1, 3, 6) GC_CM_SPECIALIZE_INCIDENCE(2, 3, 6)
GC_CM_SPECIALIZE_INCIDENCE(0, 4, 6) GC_CM_SPECIALIZE_INCIDENCE(1, 4, 6) GC_CM_SPECIALIZE_INCIDENCE(2, 4, 6) GC_CM_SPECIALIZE_INCIDENCE(3, 4, 6)
GC_CM_SPECIALIZE_INCIDENCE(0, 5, 6) GC_CM_SPECIALIZE_INCIDENCE(1, 5, 6) GC_CM_SPECIALIZE_INCIDENCE(2, 5, 6) GC_CM_SPECIALIZE_INCIDENCE(3, 5, 6) GC_CM_SPECIALIZE_INCIDENCE(4, 5, 6)
GC_CM_SPECIALIZE_INCIDENCE(0, 6, 6) GC_CM_SPECIALIZE_INCIDENCE(1, 6, 6) GC_CM_SPECIALIZE_INCIDENCE(2, 6, 6) GC_CM_SPECIALIZE_INCIDENCE(3, 6, 6) GC_CM_SPECIALIZE_INCIDENCE(4, 6, 6) GC_CM_SPECIALIZE_INCIDENCE(5, 6, 6)

#undef GC_CM_SPECIALIZE_CELL
#undef GC_CM_SPECIALIZE_DART
#undef GC_CM_SPECIALIZE_INCIDENCE

} // namespace geometrycentral
