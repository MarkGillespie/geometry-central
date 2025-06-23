#pragma once

#include "geometrycentral/utilities/element.h"
#include "geometrycentral/utilities/element_iterators.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include <cstddef>
#include <iostream>
#include <list>
#include <set>
#include <typeindex>
#include <unordered_set>

namespace geometrycentral {
namespace combinatorial_map {

// === Types and inline methods for the dart mesh pointer and datatypes
template <size_t D>
class CombinatorialMap;

template <size_t D>
class Dart;

template <size_t k, size_t D>
class Cell;

template <size_t D>
using Vertex = Cell<0, D>;

template <size_t D>
using Edge = Cell<1, D>;

template <size_t D>
using Face = Cell<2, D>;

template <size_t D, size_t E>
struct DartOrbitNavigator;

// ==========================================================
// ================        Dart        ==================
// ==========================================================

template <size_t D>
class Dart : public Element<Dart<D>, CombinatorialMap<D>> {
public:
  // Constructors
  Dart();                                      // construct an empty (null) element
  Dart(CombinatorialMap<D>* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  // Dart(const Dynamic Element<Dart>& e); // construct from a dynamic element of matching type

  // Navigators
  Vertex<D> vertex() const;
  Edge<D> edge() const;
  Face<D> face() const;

  template <size_t k>
  Cell<k, D> cell() const;

  Dart<D> partner(size_t d) const;
  Dart<D> next() const;

  bool isDead() const;
};
// == Range iterators

// All darts
template <size_t D>
struct DartRangeF {
  static bool elementOkay(const CombinatorialMap<D>& mesh, size_t ind);
  typedef Dart<D> Etype;
  typedef CombinatorialMap<D> ParentMeshT;
};

template <size_t D>
using DartSet = RangeSetBase<DartRangeF<D>>;

// ==========================================================
// ================        k-Cell        ==================
// ==========================================================

template <size_t k, size_t D>
class Cell : public Element<Cell<k, D>, CombinatorialMap<D>> {
public:
  // Constructors
  Cell();                                      // construct an empty (null) element
  Cell(CombinatorialMap<D>* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.

  // Navigators
  Dart<D> dart() const;

  std::vector<Dart<D>> adjacentDarts() const;
  template <size_t k2>
  std::vector<Cell<k2, D>> adjacentCells() const;
  std::vector<Vertex<D>> adjacentVertices() const;
  std::vector<Edge<D>> adjacentEdges() const;
  std::vector<Face<D>> adjacentFaces() const;

  bool isDead() const;
};

// == Range iterators

// All vertices
template <size_t k, size_t D>
struct CellRangeF {
  static bool elementOkay(const CombinatorialMap<D>& mesh, size_t ind);
  typedef Cell<k, D> Etype;
  typedef CombinatorialMap<D> ParentMeshT;
};
template <size_t k, size_t D>
using CellSet = RangeSetBase<CellRangeF<k, D>>;

template <size_t D>
using VertexSet = CellSet<0, D>;

template <size_t D>
using EdgeSet = CellSet<1, D>;

template <size_t D>
using FaceSet = CellSet<2, D>;
} // namespace combinatorial_map

// Declare specializations of the logic templates. This is important, because these need to be declared before any of
// the templates using them are instantiated.

// template<size_t D> inline size_t nElements<combinatorial_map::Dart<D>>(combinatorial_map::CombinatorialMap<D>* mesh);

// template<size_t D> inline size_t dataIndexOfElement<combinatorial_map::Dart<D>
// >(combinatorial_map::CombinatorialMap<D>* mesh, combinatorial_map::Dart<D> e         ); template<size_t D> struct
// ElementSetType<combinatorial_map::Dart<D>      >   { typedef combinatorial_map::DartSet<D>     type; };
// template<size_t D> inline combinatorial_map::DartSet<D>       iterateElements<combinatorial_map::Dart<D>
// >(combinatorial_map::CombinatorialMap<D>* mesh); template<size_t D> inline std::list<std::function<void(size_t)>>&
// getExpandCallbackList<combinatorial_map::Dart<D>    >(combinatorial_map::CombinatorialMap<D>* mesh); template<size_t
// D> inline std::list<std::function<void(const std::vector<size_t>&)>>&
// getPermuteCallbackList<combinatorial_map::Dart<D>     >(combinatorial_map::CombinatorialMap<D>* mesh);
// template<size_t D> inline std::string typeShortName<combinatorial_map::Dart<D>>();

template <>
inline std::string typeShortName<combinatorial_map::Dart<2>>() {
  return "d";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<0, 2>>() {
  return "v";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<1, 2>>() {
  return "e";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<2, 2>>() {
  return "f";
}

template <>
inline std::string typeShortName<combinatorial_map::Dart<3>>() {
  return "d";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<0, 3>>() {
  return "v";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<1, 3>>() {
  return "e";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<2, 3>>() {
  return "f";
}
template <>
inline std::string typeShortName<combinatorial_map::Cell<3, 3>>() {
  return "c";
}
} // namespace geometrycentral
