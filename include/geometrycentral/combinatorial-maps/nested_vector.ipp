namespace geometrycentral {
namespace combinatorial_map {
//==== Helper struct for defining nested vectors
template <std::size_t D, typename T> // Recursive template struct for nested lists
struct NestedVectorImpl {
  using type = std::vector<typename NestedVectorImpl<D - 1, T>::type>;
};

template <typename T> // Base case specialization for D = 0
struct NestedVectorImpl<0, T> {
  using type = T;
};
} // namespace combinatorial_map
} // namespace geometrycentral
