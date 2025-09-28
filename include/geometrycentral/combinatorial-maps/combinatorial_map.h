#pragma once

#include "geometrycentral/combinatorial-maps/combinatorial_map_element_types.h"
#include "geometrycentral/utilities/mesh_data.h"
#include "geometrycentral/utilities/utilities.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <array>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <vector>

// NOTE: ipp includes at bottom of file

namespace geometrycentral {
namespace combinatorial_map {

// Typedefs and forward declarations
template <size_t k, size_t D, typename T>
using CellData = MeshData<Cell<k, D>, T>;

template <size_t D, typename T>
using VertexData = MeshData<Vertex<D>, T>;

template <size_t D, typename T>
using EdgeData = MeshData<Edge<D>, T>;

template <size_t D, typename T>
using FaceData = MeshData<Face<D>, T>;

template <size_t D, typename T>
using DartData = MeshData<Dart<D>, T>;

template <size_t k1, size_t k2, size_t D, typename T>
using IncidenceData = MeshData<Incidence<k1, k2, D>, T>;

template <size_t D>
class CombinatorialMap;

template <size_t D>
struct OrbitNeighborhoodIterator;

template <size_t D>
class OrbitNeighborhood;

template <std::size_t D, typename T>
struct NestedVectorImpl;

// = std::vector<std::vector<...<T>> nested to depth D
template <std::size_t D, typename T>
using NestedVector = typename NestedVectorImpl<D, T>::type;

// ==========================================================
// ================    Combinatorial Map   ==================
// ==========================================================

template <size_t D>
class CombinatorialMap {

public:
  // Construct a simplicial complex from a list of oriented simplices
  CombinatorialMap(const std::vector<std::array<size_t, D + 1>>& simplices);

  // Construct a cell complex given as a list of (D-1)-complexes
  CombinatorialMap(const NestedVector<D, size_t>& cells);

  // Construct a cell complex given as an array of boundary matrices
  // CombinatorialMap(const std::array<SparseMatrix<int>, D>& boundaryMatrices);

  // boundaryMap[k] is the boundary map on k+1-cells. boundaryMap[k][i] is a list of k-faces incident on (k+1)-face i,
  // with their relative orientations
  // boundaryMap[1][i] must list the edges in face i in counterclockwise order
  CombinatorialMap(const std::array<std::vector<std::vector<std::pair<size_t, bool>>>, D>& boundaryMaps);

  ~CombinatorialMap();


  // Number of mesh elements of each type
  size_t nDarts() const;
  template <size_t k> // nCells<k>() counts k-cells
  size_t nCells() const;
  //== Aliases for some common k-cells
  size_t nVertices() const; // counts 0-cells
  size_t nEdges() const;    // counts 1-cells
  size_t nFaces() const;    // counts 2-cells
  size_t nCells() const;    // counts 3-cells

  template <size_t k1, size_t k2>
  size_t nIncidences() const; // WARNING: if incidences have not been used, returns 0
  //== Aliases for some common incidences
  size_t nVertexCorners() const; // counts (0, D)-incidences
  size_t nEdgeCorners() const;   // counts (1, D)-incidences
  size_t nFaceCorners() const;   // counts (0, 2)-incidences

  // Methods for range-based for loops
  // Example: for(Vertex v : mesh.vertices()) { ... }
  DartSet<D> darts();
  template <size_t k> // call cells<k>() for k-cells
  CellSet<k, D> cells();
  //== Aliases for some common k-cells
  VertexSet<D> vertices(); // 0-cells
  EdgeSet<D> edges();      // 1-cells
  FaceSet<D> faces();      // 2-cells
  CellSet<3, D> cells();   // 3-cells

  template <size_t k1, size_t k2>
  IncidenceSet<k1, k2, D> incidences();
  //== Aliases for some common incidences
  IncidenceSet<0, D, D> vertexCorners(); // (0, D)-incidences
  IncidenceSet<1, D, D> edgeCorners();   // (1, D)-incidences
  IncidenceSet<0, 2, D> faceCorners();   // (0, 2)-incidences

  template <size_t k>
  std::vector<Dart<D>> adjacentDarts(Cell<k, D> cell) const;
  template <size_t k1, size_t k2>
  std::vector<Cell<k2, D>> adjacentCells(Cell<k1, D> cell) const;
  template <size_t k>
  std::vector<Vertex<D>> adjacentVertices(Cell<k, D> cell) const;
  template <size_t k>
  std::vector<Edge<D>> adjacentEdges(Cell<k, D> cell) const;
  template <size_t k>
  std::vector<Face<D>> adjacentFaces(Cell<k, D> cell) const;
  // Returns a dart in c1 which is also in c2, or Dart<D>() if no such dart can be found
  template <size_t k1, size_t k2>
  Dart<D> adjacentDartInCell(Cell<k1, D> c1, Cell<k2, D> c2) const;

  // OrderedIncidence<a, b> is just an ordinary incidence, but with a and b ordered properly, i.e. Incidence<a,b>
  // if a < b and Incidence<b,a> otherwise
  template <size_t k1, size_t k2>
  std::vector<OrderedIncidence<k1, k2, D>> adjacentIncidences(Cell<k1, D> cell);

  template <size_t k1, size_t k2>
  std::vector<Dart<D>> adjacentDarts(Incidence<k1, k2, D> incidence) const;
  template <size_t k1, size_t k2, size_t k>
  std::vector<Cell<k, D>> adjacentCells(Incidence<k1, k2, D> cell) const;


  // Methods for accessing elements by index
  // only valid when the  mesh is compressed
  Dart<D> dart(size_t index);
  template <size_t k>
  Cell<k, D> cell(size_t index);
  //== Aliases for some common k-cells
  Vertex<D> vertex(size_t index); // 0-cells
  Edge<D> edge(size_t index);     // 1-cells
  Face<D> face(size_t index);     // 2-cells
  Cell<3, D> cell(size_t index);  // 3-cells

  template <size_t k1, size_t k2>
  Incidence<k1, k2, D> incidence(size_t index);
  //== Aliases for some common incidences
  Incidence<0, D, D> vertexCorner(); // (0, D)-incidences
  Incidence<1, D, D> edgeCorner();   // (1, D)-incidences
  Incidence<0, 2, D> faceCorner();   // (0, 2)-incidences

  DartData<D, size_t> getDartIndices();

  VertexData<D, size_t> getVertexIndices();
  EdgeData<D, size_t> getEdgeIndices();
  FaceData<D, size_t> getFaceIndices();
  CellData<3, D, size_t> getCellIndices(); // getCellIndices() indexes 3-cells
  template <size_t k>                      // getCellIndices<k>() indexes k-cells
  CellData<k, D, size_t> getCellIndices();
  template <size_t k1, size_t k2>
  IncidenceData<k1, k2, D, size_t> getIncidenceIndices();

  template <size_t k>
  SparseMatrix<int> getBoundaryMatrix(); // boundary of k-cell as a sum of (k-1)-cells

  size_t nConnectedComponents() const; // compute number of connected components [O(n)]
  // virtual bool isManifold(); // Combinatorial maps must be manifold
  // virtual bool isEdgeManifold();
  // virtual bool isOriented(); // Combinatorial maps must be oriented
  void printStatistics() const; // print info about element counts to std::cout

  // std::vector<std::vector<size_t>> getFaceVertexList();

  template <size_t E>
  std::vector<std::vector<size_t>> getCellVertexList();

  std::unique_ptr<CombinatorialMap> copy() const;
  std::unique_ptr<CombinatorialMap<D>> dual() const;
  // std::unique_ptr<ManifoldCombinatorialMap> toManifoldMesh();

  // Compress the mesh
  bool isCompressed() const;
  void compress();

  // == Mutation routines

  // Flip an edge. Edge is rotated clockwise. Return true if the edge was actually flipped (one can only flip
  // manifold, interior, triangular edges which are not incident on degree-1 vertices). Does _not_ create any new
  // elements, or cause the mesh to become decompressed.
  bool flip(Dart<D> d);

  // == Callbacks that will be invoked on mutation to keep containers/iterators/etc valid.

  // Expansion callbacks
  // Argument is the new size of the element list. Elements up to this index may now be used (but _might_ not be
  // in use immediately).
  std::list<std::function<void(size_t)>> dartExpandCallbackList;
  std::array<std::list<std::function<void(size_t)>>, D + 1> cellExpandCallbackList;
  std::map<std::pair<size_t, size_t>, std::list<std::function<void(size_t)>>> incidenceExpandCallbackList;

  // Compression callbacks
  // Argument is a permutation to a apply, such that d_new[i] = d_old[p[i]]. THe length of the permutation is hte size
  // of the new index space. Any elements with p[i] == INVALID_IND are unused in the new index space.
  std::list<std::function<void(const std::vector<size_t>&)>> dartPermuteCallbackList;
  std::array<std::list<std::function<void(const std::vector<size_t>&)>>, D + 1> cellPermuteCallbackList;
  std::map<std::pair<size_t, size_t>, std::list<std::function<void(const std::vector<size_t>&)>>>
      incidencePermuteCallbackList;

  // Mesh delete callbacks
  // (this unfortunately seems to be necessary; objects which have registered their callbacks above
  // need to know not to try to de-register them if the mesh has been deleted)
  std::list<std::function<void()>> meshDeleteCallbackList;

  // Check capacity. Needed when implementing expandable containers for mutable meshes to ensure the contain can
  // hold a sufficient number of elements before the next resize event.
  size_t nDartsCapacity() const;

  size_t nVerticesCapacity() const;
  size_t nEdgesCapacity() const;
  size_t nFacesCapacity() const;
  template <size_t k>
  size_t nCellsCapacity() const;
  template <size_t k1, size_t k2>
  size_t nIncidencesCapacity() const;

  // Return the size corresponding to the largest raw index in the mesh (except corners, see below). That is, the
  // maximum value of he.getIndex()+1 for all halfedges, etc. This may differ from `nHalfedges()` or
  // `nHalfedgesCapacity()` for non-compressed meshes. It also may differ for corners even in the case of a compressed
  // mesh. These values may change after any mutation, not just on resize events.
  size_t dartIndexSize() const;
  size_t vertexIndexSize() const;
  size_t edgeIndexSize() const;
  size_t faceIndexSize() const;
  template <size_t k>
  size_t cellIndexSize() const;

  // Loop over the neighbors of dart `d` in the orbit defining `d`s k-cell
  // can be used as `for (std::pair<Dart<D>, bool> neighbor : mesh.orbitNeighbors<k>(dart)) { ... }`,
  // where the first component is the neighboring dart, and the second component is the relative orientation of the dart
  // compared to `d`
  OrbitNeighborhood<D> orbitNeighbors(Dart<D> d, size_t k) const;

  // == Debugging, etc

  // Performs a sanity checks on dart structure; throws on fail
  // If allowDeadDarts is false, also throws if any darts are dead
  void validateConnectivity(bool allowDeadDarts = false);

  // index k-cells and fill cDartArr[k] and dCellArr[k] based off of dartMap
  void indexCells(size_t k);

  void indexIncidences(size_t k1, size_t k2);
  void ensureHaveIncidences(size_t k1, size_t k2); // helper to populate incidence arrays lazily as needed

protected:
  // Constructor used by subclasses
  CombinatorialMap();

  // Construct directly from internal arrays
  CombinatorialMap(const std::array<std::vector<size_t>, D>& dartMap);

  void constructFromBoundaryMaps(const std::array<std::vector<std::vector<std::pair<size_t, bool>>>, D>& boundaryMaps);

  // = Core arrays which hold the connectivity
  // Note: it should always be true that heFace.size() == nDartsCapacityCount, but any elements after
  // nDartsFillCount will be valid indices (in the std::vector sense), but contain uninitialized data. Similarly,
  // any std::vector<> indices corresponding to deleted elements will hold meaningless values.
  std::array<std::vector<size_t>, D> dartMap;

  size_t dartPartner(size_t iD, size_t dim) const;
  // std::vector<size_t> dVertexArr;                  // dart.vertex()
  std::array<std::vector<size_t>, D + 1> cDartArr; // cell[k].dart()
  std::array<std::vector<size_t>, D + 1> dCellArr; // dart.cell<k>().getIndex()
  std::array<std::vector<bool>, D + 1> dCellSgn;   // dart.cell<k>().orientation()

  std::map<std::pair<size_t, size_t>, std::vector<size_t>> iDartArr;      // incidence[(k1, k2)].dart()
  std::map<std::pair<size_t, size_t>, std::vector<size_t>> dIncidenceArr; // dart.incidence<k1, k2>.getIndex()

  // Auxilliary arrays which cache other useful information

  // Track element counts (can't rely on rawVertices.size() after deletions have made the list sparse). These are the
  // actual number of valid elements, not the size of the buffer that holds them.
  size_t nDartsCount = 0;
  std::array<size_t, D + 1> nCellsCount{};
  std::map<std::pair<size_t, size_t>, size_t> nIncidencesCount;

  // == Track the capacity and fill size of our buffers.
  // These give the capacity of the currently allocated buffer.
  // Note that this is _not_ defined to be std::vector::capacity(), it's the largest size such that arr[i] is legal (aka
  // arr.size()).
  size_t nDartsCapacityCount = 0;                  // will always be even if implicit twin
  std::array<size_t, D + 1> nCellsCapacityCount{}; // will always be even if implicit twin
  std::map<std::pair<size_t, size_t>, size_t> nIncidencesCapacityCount;

  // These give the number of filled elements in the currently allocated buffer. This will also be the maximal index of
  // any element (except the weirdness of boundary loop faces). As elements get marked dead, nVerticesCount decreases
  // but nVertexFillCount does not (etc), so it denotes the end of the region in the buffer where elements have been
  // stored.
  size_t nDartsFillCount = 0; // must always be even if implicit twin
  std::array<size_t, D + 1> nCellsFillCount{};
  std::map<std::pair<size_t, size_t>, size_t> nIncidencesFillCount;

  // The mesh is _compressed_ if all of the index spaces are dense. E.g. if thare are |V| vertices, then the vertices
  // are densely indexed from 0 ... |V|-1 (and likewise for the other elements). The mesh can become not-compressed as
  // deletions mark elements with tombstones--this is how we support constant time deletion.
  // Call compress() to re-index and return to usual dense indexing.
  bool isCompressedFlag = true;

  uint64_t modificationTick = 1; // Increment every time the mesh is mutated in any way. Used to track staleness.

  // Hide copy and move constructors, we don't wanna mess with that
  CombinatorialMap(const CombinatorialMap& other) = delete;
  CombinatorialMap& operator=(const CombinatorialMap& other) = delete;
  CombinatorialMap(CombinatorialMap&& other) = delete;
  CombinatorialMap& operator=(CombinatorialMap&& other) = delete;

  // Used to resize the halfedge mesh. Expands and shifts vectors as necessary.
  Dart<D> getNewDart();
  void allocateDarts(size_t n); // ensure we have space for n more darts

  template <size_t k>
  Cell<k, D> getNewCell();

  size_t getNewCellIndex(size_t k); // equal to getNewCell<k>().index()
  size_t getNewIncidenceIndex(std::pair<size_t, size_t> k1k2);

  void allocateCells(size_t k, size_t n); // ensure we have space for n more k-cells
  template <size_t k>
  void allocateCells(size_t n); // ensure we have space for n more k-cells

  // Detect dead elements
  bool dartIsDead(size_t iD) const;
  template <size_t k>
  bool cellIsDead(size_t iC) const;
  bool cellIsDead(size_t k, size_t iC) const;
  template <size_t k1, size_t k2>
  bool incidenceIsDead(size_t iI) const;

  // Deletes leave tombstones, which can be cleaned up with compress().
  // Note that these routines merely mark the element as dead. The caller should hook up connectivity to exclude these
  // elements before invoking.
  void deleteElement(Dart<D> d); // can't use for implicit twin

  // Compression helpers
  void compressDarts();

  // = =Helpers for mutation methods and similar things

  void initializeDartNeighbors();
  void copyInternalFields(CombinatorialMap& target) const;

  // replace values of i in arr with oldToNew[i] (skipping INVALID_IND)
  void updateValues(std::vector<size_t>& arr, const std::vector<size_t>& oldToNew);

  // Elements need direct access in to members to traverse
  friend class Dart<D>;
  friend struct DartRangeF<D>;

  template <size_t k, size_t D1>
  friend class Cell;
  template <size_t k, size_t D1>
  friend struct CellRangeF;

  template <size_t k1, size_t k2, size_t D1>
  friend class Incidence;
  template <size_t k1, size_t k2, size_t D1>
  friend struct IncidenceRangeF;
};

template <size_t D>
std::vector<Dart<D>> incidenceNeighboringDarts(Dart<D> d, size_t k1, size_t k2, bool verbose = false);

// helpers
namespace unionfind {
size_t findRoot(size_t x, std::vector<size_t>& parent, std::vector<bool>& sharesParentSign);
void unite(size_t x, size_t y, bool samesign, std::vector<size_t>& parent, std::vector<bool>& sharesParentSign,
           std::vector<size_t>& rank);
} // namespace unionfind

} // namespace combinatorial_map
} // namespace geometrycentral

// clang-format off
// preserve ordering
// #include "geometrycentral/combinatorial-maps/dart_logic_templates.ipp"
#include "geometrycentral/combinatorial-maps/nested_vector.ipp"
#include "geometrycentral/combinatorial-maps/combinatorial_map.ipp"
#include "geometrycentral/combinatorial-maps/combinatorial_map_element_types.ipp"
// clang-format on
