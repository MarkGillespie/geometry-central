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

// Represents the incidence of a k1-cell inside a k2-cell for k1 < k2
// E.g. an Incidence<0, 2, D> is a corner of a face
template <size_t k1, size_t k2, size_t D>
class Incidence;

// OrderedIncidence<a, b> creates an Incidence with a and b in the right order, i.e. Incidence<a,b> if a < b, and
// Incidence<b, a> otherwise
template <size_t a, size_t b, size_t D>
using OrderedIncidence = typename std::conditional<(a < b), Incidence<a, b, D>, Incidence<b, a, D>>::type;

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

  //=== Navigators
  template <size_t k>
  Cell<k, D> cell() const; // k-cell
  //== Aliases for some common k-cells
  Vertex<D> vertex() const; // 0-cell
  Edge<D> edge() const;     // 1-cell
  Face<D> face() const;     // 2-cell
  Cell<3, D> cell() const;  // 3-cell

  Vertex<D> tailVertex() const; // same as .vertex()
  Vertex<D> tipVertex() const;  // same as .next().vertex()

  template <size_t k1, size_t k2>
  Incidence<k1, k2, D> incidence() const;
  //== Aliases for some common incidences
  Incidence<0, D, D> vertexCorner() const; // (0, D)-incidence
  Incidence<1, D, D> edgeCorner() const;   // (1, D)-incidence
  Incidence<0, 2, D> faceCorner() const;   // (0, 2)-incidence

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
// ================        k-Cell        ====================
// ==========================================================

template <size_t k, size_t D>
class Cell : public Element<Cell<k, D>, CombinatorialMap<D>> {
public:
  // Constructors
  Cell();                                      // construct an empty (null) element
  Cell(CombinatorialMap<D>* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.
  Cell(CombinatorialMap<D>* mesh, size_t ind,
       bool orientation); // construct pointing to the i'th element of that type on a mesh with specified orientation

  //=== Navigators
  Dart<D> dart() const;
  // finds a dart in this cell which is also in the input cell, or Dart<D>() if no such dart is found
  template <size_t k2>
  Dart<D> dartInCell(Cell<k2, D> cell) const;

  std::vector<Dart<D>> adjacentDarts() const;
  template <size_t k2>
  std::vector<Cell<k2, D>> adjacentCells() const; // adjacentCells<k> gives k-cells
  //== Aliases for some common k-cells
  std::vector<Vertex<D>> adjacentVertices() const; // 0-cells
  std::vector<Edge<D>> adjacentEdges() const;      // 1-cells
  std::vector<Face<D>> adjacentFaces() const;      // 2-cells
  std::vector<Cell<3, D>> adjacentCells() const;   //  3-cells

  template <size_t k2> // function for iterating over adjacent incidences with higher-dimensional cells
  std::vector<OrderedIncidence<k, k2, D>> adjacentIncidences() const;
  //== Aliases for some common incidences
  // Warning: only defined for k-cells which are part of the incidence. e.g., you can call adjacentEdgeCorners() on and
  // edge or a top-dimensional cell, but not on a 0-cell
  std::vector<Incidence<0, D, D>> adjacentVertexCorners() const; // (0, D)-incidences
  std::vector<Incidence<1, D, D>> adjacentEdgeCorners() const;   // (1, D)-incidences
  std::vector<Incidence<0, 2, D>> adjacentFaceCorners() const;   // (0, 2)-incidences

  bool isDead() const;
  bool isBoundary() const; // returns true if the cell is totally contained in the mesh boundary

  bool orientation() const;
  void flipOrientation();
  void setOrientation(bool orientation);

  int sign() const;

  bool orientationInCell(Cell<k + 1, D> c) const;
  bool orientationInCell(Cell<k - 1, D> c) const;
  int signInCell(Cell<k + 1, D> c) const; // entry in boundary_{k+1} matrix
  int signInCell(Cell<k - 1, D> c) const; // entry in boundary_k matrix

protected:
  bool mOrientation = true;
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


// ==========================================================
// ==============       k-Cell Incidence       ==============
// ==========================================================

// Class to represent ``incidence'' between a k1 cell and a k2 cell. For instance, a corners of a triangle is a
// (0,2)-incidence, and a ``hinges'' of a tetrahedron is a (1,3)-incidence

// Warning: these can be weirdly-behaved on boundary edges of triangle meshes, where there is only one dart on each
// edge. For example, a single triangle only has three (0, 1)-incidences instead of 6, since this implementation of
// (0,1) incidences identifies an incidence with the set of darts common to both k-cells.

template <size_t k1, size_t k2, size_t D>
class Incidence : public Element<Incidence<k1, k2, D>, CombinatorialMap<D>> {
  static_assert(k1 < k2, "Incidence requires k1 < k2");
  static_assert(k2 <= D, "Incidence requires k2 <= D");

public:
  // Constructors
  Incidence();                                      // construct an empty (null) element
  Incidence(CombinatorialMap<D>* mesh, size_t ind); // construct pointing to the i'th element of that type on a mesh.

  // Navigators
  Dart<D> dart() const;

  template <size_t k>
  Cell<k, D> cell() const;  // adjacent k-cell. Only defined for k=k1 or k=k2
  Vertex<D> vertex() const; // convenience alias for cell<0>()
  Edge<D> edge() const;     // convenience alias for cell<1>()
  Face<D> face() const;     // convenience alias for cell<2>()
  Cell<3, D> cell() const;  // convenience alias for cell<3>()

  std::vector<Dart<D>> adjacentDarts() const;
  template <size_t k>
  std::vector<Cell<k, D>> adjacentCells() const; // adjacentCells<k> gives k-cells
  //== Aliases for some common k-cells
  std::vector<Vertex<D>> adjacentVertices() const; // 0-cells
  std::vector<Edge<D>> adjacentEdges() const;      // 1-cells
  std::vector<Face<D>> adjacentFaces() const;      // 2-cells
  std::vector<Cell<3, D>> adjacentCells() const;   //  3-cells

  bool isDead() const;
  bool isBoundary() const;
};

// == Range iterators

// All vertices
template <size_t k1, size_t k2, size_t D>
struct IncidenceRangeF {
  static bool elementOkay(const CombinatorialMap<D>& mesh, size_t ind);
  typedef Incidence<k1, k2, D> Etype;
  typedef CombinatorialMap<D> ParentMeshT;
};
template <size_t k1, size_t k2, size_t D>
using IncidenceSet = RangeSetBase<IncidenceRangeF<k1, k2, D>>;
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
