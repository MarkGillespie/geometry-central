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
inline Dart<D> Dart<D>::partner(size_t d) const {
  size_t partnerInd = this->mesh->dartPartner(this->ind, d);
  return (partnerInd == INVALID_IND) ? *this : Dart<D>(this->mesh, partnerInd);
};

template <size_t D>
inline Dart<D> Dart<D>::next() const {
  return partner(0);
};

template <size_t D>
template <size_t k>
inline Cell<k, D> Dart<D>::cell() const {
  static_assert(k <= D, "cell dimension k must be less than or equal to complex dimension D");
  return Cell<k, D>(this->mesh, this->mesh->dCellArr[k][this->ind], this->mesh->dCellSgn[k][this->ind]);
};

template <size_t D>
inline Vertex<D> Dart<D>::vertex() const {
  return cell<0>();
};

template <size_t D>
inline Edge<D> Dart<D>::edge() const {
  return cell<1>();
};

template <size_t D>
inline Face<D> Dart<D>::face() const {
  return cell<2>();
};

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
};

template <size_t k1, size_t D>
inline std::vector<Dart<D>> Cell<k1, D>::adjacentDarts() const {
  return this->mesh->adjacentDarts(*this);
}

template <size_t k1, size_t D>
template <size_t k2>
inline std::vector<Cell<k2, D>> Cell<k1, D>::adjacentCells() const {
  return this->mesh->template adjacentCells<k1, k2>(*this);
}

template <size_t k1, size_t D>
inline std::vector<Vertex<D>> Cell<k1, D>::adjacentVertices() const {
  return adjacentCells<0>();
}

template <size_t k1, size_t D>
inline std::vector<Edge<D>> Cell<k1, D>::adjacentEdges() const {
  return adjacentCells<1>();
}

template <size_t k1, size_t D>
inline std::vector<Face<D>> Cell<k1, D>::adjacentFaces() const {
  return adjacentCells<2>();
}

template <size_t k, size_t D>
inline bool Cell<k, D>::isDead() const {
  return this->mesh->cellIsDead<k>(this->ind);
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

// Range iterators
template <size_t k, size_t D>
inline bool CellRangeF<k, D>::elementOkay(const CombinatorialMap<D>& mesh, size_t ind) {
  return !(mesh.template cellIsDead<k>(ind));
}


} // namespace combinatorial_map


template <>
struct ElementSetType<combinatorial_map::Vertex<2>> {
  typedef combinatorial_map::VertexSet<2> type;
};
template <>
struct ElementSetType<combinatorial_map::Vertex<3>> {
  typedef combinatorial_map::VertexSet<3> type;
};
template <>
struct ElementSetType<combinatorial_map::Dart<2>> {
  typedef combinatorial_map::DartSet<2> type;
};
template <>
struct ElementSetType<combinatorial_map::Dart<3>> {
  typedef combinatorial_map::DartSet<3> type;
};

template <>
inline size_t nElements<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nVertices();
}
template <>
inline size_t nElements<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nVertices();
}
template <>
inline size_t nElements<combinatorial_map::Edge<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nEdges();
}
template <>
inline size_t nElements<combinatorial_map::Edge<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nEdges();
}
template <>
inline size_t nElements<combinatorial_map::Face<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nFaces();
}
template <>
inline size_t nElements<combinatorial_map::Face<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nFaces();
}
template <>
inline size_t nElements<combinatorial_map::Cell<3, 3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nCells<3>();
}
template <>
inline size_t nElements<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nDarts();
}
template <>
inline size_t nElements<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nDarts();
}

template <>
inline size_t elementCapacity<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nVerticesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nVerticesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Edge<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nEdgesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Edge<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nEdgesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Face<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nFacesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Face<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nFacesCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Cell<3, 3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nCellsCapacity<3>();
}
template <>
inline size_t elementCapacity<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->nDartsCapacity();
}
template <>
inline size_t elementCapacity<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->nDartsCapacity();
}

template <>
inline size_t dataIndexOfElement<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh,
                                                               combinatorial_map::Vertex<2> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh,
                                                               combinatorial_map::Vertex<3> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Edge<2>>(combinatorial_map::CombinatorialMap<2>* mesh,
                                                             combinatorial_map::Edge<2> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Edge<3>>(combinatorial_map::CombinatorialMap<3>* mesh,
                                                             combinatorial_map::Edge<3> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Face<2>>(combinatorial_map::CombinatorialMap<2>* mesh,
                                                             combinatorial_map::Face<2> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Face<3>>(combinatorial_map::CombinatorialMap<3>* mesh,
                                                             combinatorial_map::Face<3> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Cell<3, 3>>(combinatorial_map::CombinatorialMap<3>* mesh,
                                                                combinatorial_map::Cell<3, 3> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh,
                                                             combinatorial_map::Dart<2> e) {
  return e.getIndex();
}
template <>
inline size_t dataIndexOfElement<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh,
                                                             combinatorial_map::Dart<3> e) {
  return e.getIndex();
}


template <>
inline combinatorial_map::VertexSet<2>
iterateElements<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->vertices();
}
template <>
inline combinatorial_map::VertexSet<3>
iterateElements<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->vertices();
}
template <>
inline combinatorial_map::DartSet<2>
iterateElements<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->darts();
}
template <>
inline combinatorial_map::DartSet<3>
iterateElements<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->darts();
}

template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellExpandCallbackList[0];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellExpandCallbackList[0];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Edge<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellExpandCallbackList[1];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Edge<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellExpandCallbackList[1];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Face<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellExpandCallbackList[2];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Face<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellExpandCallbackList[2];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Cell<3, 3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellExpandCallbackList[3];
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->dartExpandCallbackList;
}
template <>
inline std::list<std::function<void(size_t)>>&
getExpandCallbackList<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->dartExpandCallbackList;
}

template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Vertex<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellPermuteCallbackList[0];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Vertex<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellPermuteCallbackList[0];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Edge<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellPermuteCallbackList[1];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Edge<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellPermuteCallbackList[1];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Face<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->cellPermuteCallbackList[2];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Face<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellPermuteCallbackList[2];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Cell<3, 3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->cellPermuteCallbackList[3];
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Dart<2>>(combinatorial_map::CombinatorialMap<2>* mesh) {
  return mesh->dartPermuteCallbackList;
}
template <>
inline std::list<std::function<void(const std::vector<size_t>&)>>&
getPermuteCallbackList<combinatorial_map::Dart<3>>(combinatorial_map::CombinatorialMap<3>* mesh) {
  return mesh->dartPermuteCallbackList;
}

// template <size_t D>
// inline size_t elementCapacity<combinatorial_map::Vertex<D>>(combinatorial_map::CombinatorialMap<D>* mesh) {
//   return mesh->nVerticesCapacity();
// }

// template <size_t D>
// inline size_t elementCapacity<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh) {
//   return mesh->nDartsCapacity();
// }

} // namespace geometrycentral
