/**
 * Current a2d way to create the layout and multiarray of shape (n, N1, N2)
 * is:
 * Layout<N1, N2> layout(n); MultiArray<T, Layout<N1, N2>> array(layout);
 *
 * Kokkos' way to create the array with same layout is:
 *
 * View<T*[N1][N2], LayoutType, MemoryType> view("label", n);
 */
#ifndef A2D_MULTI_ARRAY_H
#define A2D_MULTI_ARRAY_H

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <type_traits>

#include "a2dlayout.h"
#include "a2dmemory.h"
#include "a2dobjs.h"

namespace A2D {

/*
  Unpack the argument list to get the extent of indices
*/
template <int r, index_t... dims>
struct ___get_extent;

template <index_t dim0, index_t... dims>
struct ___get_extent<0, dim0, dims...> {
  static const index_t extent = dim0;
};

template <int r, index_t dim0, index_t... dims>
struct ___get_extent<r, dim0, dims...> {
  static const index_t extent = ___get_extent<r - 1, dims...>::extent;
};

template <index_t... dims>
struct __get_size;

template <>
struct __get_size<> {
  static const index_t size = 1;  // This is the multiplier of the runtime
                                  // dimension (dim0) - so it's 1 rather than 0
};

template <index_t dim0, index_t... dims>
struct __get_size<dim0, dims...> {
  static const index_t size = dim0 * __get_size<dims...>::size;
};

template <index_t dim0>
struct __get_size<dim0> {
  static const index_t size = dim0;
};

/**
 * Parse the template parameter pack to get the type that kokkos requires:
 * <T, N1, N2, N3>   ->   T*[N1][N2][N3]
 *
 * Usage:
 *   ParsePack<T, N1, N2, N3, ...>::kokkos_type = T*[N1][N2][N3]...
 *
 */
template <typename T, index_t... dims>
struct ParsePack {
  using kokkos_type = T*;
};

template <typename T, index_t... dims, index_t dim>
struct ParsePack<T, dim, dims...> {
  using nested = ParsePack<T, dims...>;
  using kokkos_type = typename nested::kokkos_type[dim];
};

/*
  Fortran ordering

  Given an entry (i, j, k) in a multi-dimensional array with shape (dim0, dim1,
  dim2), the fortran ordering stores entry A(i, j, k) in the location i + dim0 *
  (j + dim1 * k).
*/
template <index_t... dims>
class FLayout {
 public:
#ifdef A2D_USE_KOKKOS
  template <typename T>
  using KokkosTypeTemplate = typename ParsePack<T, dims...>::kokkos_type;
  using KokkosLayout = Kokkos::LayoutLeft;
#endif

  index_t dim0;  // the leading dimension (set at runtime)
  static const index_t static_size =
      __get_size<dims...>::size;  // total number of entries in the compile-time
                                  // dimensions
  static const index_t rank = sizeof...(dims) + 1;  // total number of
                                                    // dimensions (fixed and
                                                    // variable)

  /**
   * @brief Default constructor
   */
  A2D_INLINE_FUNCTION FLayout() {}

  /**
   * @brief Regular constructor
   */
  FLayout(const index_t dim0) : dim0(dim0), static_extents{dims...} {
    size = dim0 * static_size;
  }

  /**
   * @brief Copy constructor (light-weight shallow copy)
   */
  A2D_INLINE_FUNCTION FLayout(const FLayout<dims...>& src)
      : dim0(src.dim0), static_extents{dims...} {
    size = dim0 * static_size;
  }

  /**
   * @brief Destructor
   */
  A2D_INLINE_FUNCTION ~FLayout() {}

  /**
   * @brief Assignment operator
   */
  A2D_INLINE_FUNCTION FLayout<dims...>& operator=(const FLayout<dims...>& src) {
    dim0 = src.dim0;
    size = src.size;
    for (index_t i = 0; i != rank - 1; i++) {
      static_extents[i] = src.static_extents[i];
    }
    return *this;
  }

  /**
   * @brief Get number of dimensions of the multiarray layout
   */
  A2D_INLINE_FUNCTION static index_t get_rank() { return rank; }

  /**
   * @brief compute dim * dim1 * dim2 * dim3 ...
   */
  A2D_INLINE_FUNCTION static index_t get_size(index_t dim) {
    return dim * __get_size<dims...>::size;
  }

  /**
   * @brief Get number of entries of the multiarray
   */
  A2D_INLINE_FUNCTION index_t get_size() { return size; }

  /**
   * @brief Get the extent given the dimension index
   */
  A2D_INLINE_FUNCTION index_t get_extent(index_t index) const {
    if (index == 0) {
      return dim0;
    }
    return static_extents[index - 1];
  }

  template <class Idx>
  A2D_INLINE_FUNCTION index_t compute_index(Idx i1) const {
    return i1;
  }

  /**
   * @brief Find memory location index given the multidimensional coordinates
   * (i1, i2, i3, ...)
   */
  template <class Idx, class... IdxType>
  A2D_INLINE_FUNCTION index_t compute_index(Idx i1, IdxType... idx) const {
    return i1 + dim0 * __compute_index<0>(idx...);
  }

  /**
   * @brief Find the memory location offset between (i1, i2, i3, ...) and
   * (i1, 0, 0, ...) divided by the first (runtime) dimension.
   * This is a static method used in MultiArraySlice.
   */
  template <class... IdxType>
  A2D_INLINE_FUNCTION static index_t compute_slice_offset(IdxType... idx) {
    return __compute_index<0>(idx...);
  }

 private:
  index_t size;
  index_t static_extents[rank - 1];

  template <int r, class Idx, class... IdxType>
  A2D_INLINE_FUNCTION static index_t __compute_index(Idx i, IdxType... idx) {
    return i +
           ___get_extent<r, dims...>::extent * __compute_index<r + 1>(idx...);
  }

  template <int r, class Idx>
  A2D_INLINE_FUNCTION static index_t __compute_index(Idx i) {
    return i;
  }
};

/*
  C ordering

  Given an entry (i, j, k) in a multi-dimensional array with shape (dim0, dim1,
  dim2), the c ordering stores entry A(i, j, k) in the location (i * dim1 + j) *
  dim2 + k.
*/
template <index_t... dims>
class CLayout {
 public:
#ifdef A2D_USE_KOKKOS
  template <typename T>
  using KokkosTypeTemplate = typename ParsePack<T, dims...>::kokkos_type;
  using KokkosLayout = Kokkos::LayoutRight;
#endif

  index_t dim0;  // the leading dimension (set at runtime)
  static const index_t static_size =
      __get_size<dims...>::size;  // total number of entries in the compile-time
                                  // dimensions
  static const index_t rank = sizeof...(dims) + 1;  // total number of
                                                    // dimensions (fixed and
                                                    // variable)

  /**
   * @brief Default constructor
   */
  A2D_INLINE_FUNCTION CLayout() {}

  /**
   * @brief Regular constructor
   */
  CLayout(const index_t dim0) : dim0(dim0), static_extents{dims...} {
    size = dim0 * static_size;
  }

  /**
   * @brief Copy constructor (light-weight shallow copy)
   */
  A2D_INLINE_FUNCTION CLayout(const CLayout<dims...>& src)
      : dim0(src.dim0), static_extents{dims...} {
    size = dim0 * static_size;
  }

  /**
   * @brief Destructor
   */
  A2D_INLINE_FUNCTION ~CLayout() {}

  /**
   * @brief Assignment operator
   */
  A2D_INLINE_FUNCTION CLayout<dims...>& operator=(const CLayout<dims...>& src) {
    dim0 = src.dim0;
    size = src.size;
    for (index_t i = 0; i != rank - 1; i++) {
      static_extents[i] = src.static_extents[i];
    }
    return *this;
  }

  /**
   * @brief Get number of dimensions of the multiarray layout
   */
  A2D_INLINE_FUNCTION static index_t get_rank() { return rank; }

  /**
   * @brief compute dim * dim1 * dim2 * dim3 ...
   */
  A2D_INLINE_FUNCTION static index_t get_size(index_t dim) {
    return dim * __get_size<dims...>::size;
  }

  /**
   * @brief Get number of entries of the multiarray
   */
  A2D_INLINE_FUNCTION index_t get_size() { return size; }

  /**
   * @brief Get the extent given the dimension index
   */
  A2D_INLINE_FUNCTION index_t get_extent(index_t index) const {
    if (index == 0) {
      return dim0;
    }
    return static_extents[index - 1];
  }

  template <class Idx>
  A2D_INLINE_FUNCTION index_t compute_index(Idx i1) const {
    return i1;
  }

  /**
   * @brief Find memory location index given the multidimensional coordinates
   * (i1, i2, i3, ...)
   */
  template <class Idx, class... IdxType>
  A2D_INLINE_FUNCTION static index_t compute_index(Idx i1, IdxType... idx) {
    return __compute_index<0>(i1, idx...);
  }

  /**
   * @brief Find the memory location offset between (i1, i2, i3, ...) and
   * (i1, 0, 0, ...). This is a static method used in MultiArraySlice.
   */
  template <class... IdxType>
  A2D_INLINE_FUNCTION static index_t compute_slice_offset(IdxType... idx) {
    return __compute_index<0>(0, idx...);
  }

 private:
  index_t size;
  index_t static_extents[rank - 1];

  template <int r, class Idx, class... IdxType>
  A2D_INLINE_FUNCTION static index_t __compute_index(const index_t index, Idx i,
                                                     IdxType... idx) {
    return __compute_index<r + 1>(index * ___get_extent<r, dims...>::extent + i,
                                  idx...);
  }

  template <int r, class Idx>
  A2D_INLINE_FUNCTION static index_t __compute_index(const index_t index,
                                                     Idx i) {
    return index * ___get_extent<r, dims...>::extent + i;
  }
};

#ifdef A2D_USE_KOKKOS
/**
 * @brief A multi-dimensional array wrapper with Kokkos View backend
 *
 * @tparam T data type
 * @tparam Layout CLayout or FLayout
 */
template <typename T, class Layout, class MemSpace = Kokkos::HostSpace>
class MultiArray {
 public:
  using type = T;
  using ViewType = Kokkos::View<typename Layout::template KokkosTypeTemplate<T>,
                                typename Layout::KokkosLayout, MemSpace>;

  /**
   * @brief Default constructor
   *
   * This creates skeleton only, no memory allocation
   */
  MultiArray() { data = nullptr; }

  /**
   * @brief Constructor, allocate the array managed by kokkos
   *
   * @param layout CLayout or FLayout instance
   * @param data_ pointer to raw data
   */
  MultiArray(Layout layout) : layout(layout) {
    index_t dim0 = layout.get_extent(0);
    view = ViewType("default_label", dim0);
    data = view.data();
  }

  /**
   * @brief Copy constructor: shallow copy is performed
   */
  A2D_INLINE_FUNCTION MultiArray(const MultiArray<T, Layout>& src) {
    layout = src.layout;
    view = src.view;
    data = src.data;
  }

  /**
   * @brief Assignment operator: shallow copy is performed
   */
  A2D_INLINE_FUNCTION MultiArray& operator=(const MultiArray<T, Layout>& src) {
    layout = src.layout;
    view = src.view;
    data = src.data;
    return *this;
  }

  /**
   * @brief Destructor - memory is managed by Kokkos, hence no explicit
   * delete/free is needed
   */
  A2D_INLINE_FUNCTION ~MultiArray() {}

  /*
    Constant access to array elements, note: defining member function as const
    gives it better flexibility
  */
  template <class... IdxType>
  A2D_INLINE_FUNCTION T& operator()(IdxType... idx) const {
    return view(idx...);
  }

  template <class... IdxType>
  A2D_INLINE_FUNCTION T& operator[](IdxType... idx) const {
    return view(idx...);
  }

  /**
   * @brief Iterator: begin
   */
  T* begin() const { return data; }

  /**
   * @brief Iterator: end
   */
  T* end() const { return data + layout.get_size(); }

  /*
    Get the rank of the the multi-dimensional array
  */
  A2D_INLINE_FUNCTION index_t get_rank() const { return layout.get_rank(); }

  /*
    Get the extent of one of the dimensions
  */
  A2D_INLINE_FUNCTION index_t extent(index_t index) const {
    return layout.get_extent(index);
  }

  /*
    Zero all elements in the array
  */
  A2D_INLINE_FUNCTION void zero() {
    const index_t len = layout.get_size();
    for (index_t i = 0; i != len; i++) {
      data[i] = T(0);
    }
  }

  /*
    Fill all the values in the array with the specified value
  */
  A2D_INLINE_FUNCTION void fill(T value) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i != len; i++) {
      data[i] = value;
    }
  }

  /*
    Copy elements from the source to this vector
  */
  A2D_INLINE_FUNCTION void copy(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i != len; i++) {
      data[i] = src.data[i];
    }
  }

  /*
    Copy elements from the source to this vector
  */
  A2D_INLINE_FUNCTION void scale(T alpha) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] *= alpha;
    }
  }

  /*
    Duplicate the array
  */
  MultiArray<T, Layout>* duplicate() {
    const index_t len = layout.get_size();
    MultiArray<T, Layout>* array = new MultiArray<T, Layout>(layout);
    std::copy(data, data + len, array->data);
    return array;
  }

  /*
    Set a random seed for the data
  */
  void random(T lower = -1.0, T upper = 1.0) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] = lower + ((upper - lower) * (1.0 * rand())) / (1.0 * RAND_MAX);
    }
  }

  /*
    Take the dot product with the source vector data
  */
  A2D_INLINE_FUNCTION T dot(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    T result = T(0);
    for (index_t i = 0; i < len; i++) {
      result += data[i] * src.data[i];
    }
    return result;
  }

  /*
    Norm of the array
  */
  A2D_INLINE_FUNCTION T norm() {
    const index_t len = layout.get_size();
    T result = T(0);
    for (index_t i = 0; i < len; i++) {
      result += data[i] * data[i];
    }
    return sqrt(result);
  }

  /*
    Axpy: this = alpha * x + this
  */
  A2D_INLINE_FUNCTION void axpy(T alpha, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] += alpha * x.data[i];
    }
  }

  /*
    Axpby: this = alpha * x + beta * this
  */
  A2D_INLINE_FUNCTION void axpby(T alpha, T beta, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] = alpha * x.data[i] + beta * data[i];
    }
  }

  Layout layout;
  ViewType view;
  T* data;
};
#else

/**
 * @brief A multi-dimensional array that behaves like the shared_ptr
 */
template <typename T, class Layout>
class MultiArray {
 public:
  using type = T;
  using SharedPtrType = std::shared_ptr<T[]>;

  Layout layout;
  T* data;  // raw pointer

  /**
   * @brief Default constructor
   *
   * This creates skeleton only, no memory allocation
   */
  MultiArray() { data = nullptr; }

  /**
   * @brief Constructor, allocate the array managed by shared_ptr
   */
  MultiArray(Layout layout) : layout(layout) {
    data_sp = SharedPtrType(new T[layout.get_size()]);
    data = data_sp.get();
    zero();
  }

  /**
   * @brief Constructor, create a MultiArray with unmanaged data
   * values
   */
  MultiArray(Layout layout, T* vals) : layout(layout) {
    printf(
        "[MultiArray] creating multiarray with external memory (not managed by "
        "smart pointer)\n");
    data = vals;
    data_sp = nullptr;
  }

  /**
   * @brief Copy constructor: shallow copy is performed
   */
  MultiArray(const MultiArray<T, Layout>& src) {
    layout = src.layout;
    data = src.data;
    data_sp = src.data_sp;
  }

  /**
   * @brief Assignment operator: shallow copy is performed
   */
  MultiArray& operator=(const MultiArray<T, Layout>& src) {
    layout = src.layout;
    data_sp = src.data_sp;
    data = src.data;
    return *this;
  }

  /**
   * @brief Destructor - memory is managed by the smart pointer (shared_ptr),
   * hence no explicit delete/free is needed
   */
  ~MultiArray() {}

  /*
    Constant access to array elements, note: defining member function as const
    gives it better flexibility
  */
  template <class... IdxType>
  T& operator()(IdxType... idx) const {
    return data[layout.compute_index(idx...)];
  }

  template <class... IdxType>
  T& operator[](IdxType... idx) const {
    return data[layout.compute_index(idx...)];
  }

  /**
   * @brief Iterator: begin
   */
  T* begin() const { return data; }

  /**
   * @brief Iterator: end
   */
  T* end() const { return data + layout.get_size(); }

  /*
    Get the rank of the the multi-dimensional array
  */
  index_t get_rank() const { return layout.get_rank(); }

  /*
    Get the extent of one of the dimensions
  */
  index_t extent(index_t index) const { return layout.get_extent(index); }

  /*
    Zero all elements in the array
  */
  void zero() {
    const index_t len = layout.get_size();
    std::fill(data, data + len, T(0));
  }

  /*
    Fill all the values in the array with the specified value
  */
  void fill(T value) {
    const index_t len = layout.get_size();
    std::fill(data, data + len, value);
  }

  /*
    Copy elements from the source to this vector
  */
  void copy(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    std::copy(src.data, src.data + len, data);
  }

  /*
    Copy elements from the source to this vector
  */
  void scale(T alpha) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] *= alpha;
    }
  }

  /*
    Duplicate the array
  */
  MultiArray<T, Layout>* duplicate() {
    const index_t len = layout.get_size();
    MultiArray<T, Layout> array = MultiArray<T, Layout>(layout);
    std::copy(data, data + len, array->data);
    return array;
  }

  /*
    Set a random seed for the data
  */
  void random(T lower = -1.0, T upper = 1.0) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] =
          lower + ((upper - lower) * (1.0 * std::rand())) / (1.0 * RAND_MAX);
    }
  }

  /*
    Take the dot product with the source vector data
  */
  T dot(MultiArray<T, Layout>& src) {
    const index_t len = layout.get_size();
    return std::inner_product(data, data + len, src.data, T(0));
  }

  /*
    Norm of the array
  */
  T norm() {
    const index_t len = layout.get_size();
    return std::sqrt(std::inner_product(data, data + len, data, T(0)));
  }

  /*
    Axpy: this = alpha * x + this
  */
  void axpy(T alpha, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] += alpha * x.data[i];
    }
  }

  /*
    Axpby: this = alpha * x + beta * this
  */
  void axpby(T alpha, T beta, MultiArray<T, Layout>& x) {
    const index_t len = layout.get_size();
    for (index_t i = 0; i < len; i++) {
      data[i] = alpha * x.data[i] + beta * data[i];
    }
  }

 private:
  SharedPtrType data_sp;  // smart pointer
};
#endif

/**
 * @brief The primary template that is never used directly, one of the two
 * partial specializations below will be used instead.
 */
template <typename T, class Layout>
class MultiArraySlice;

/**
 * @brief A partially-specialized MultiArraySlice for FLayout
 */
template <typename T, index_t... dims>
class MultiArraySlice<T, FLayout<dims...>> {
 public:
  using Layout = FLayout<dims...>;

  template <class IdxType>
  A2D_INLINE_FUNCTION MultiArraySlice(
      const MultiArray<T, FLayout<dims...>>& array, const IdxType idx)
      : dim0(array.layout.dim0) {
    data = &array.data[idx];
  }

  /*
    Constant access to array elements, note: defining member function as const
    gives it better flexibility
  */
  template <class... IdxType>
  A2D_INLINE_FUNCTION T& operator()(IdxType... idx) const {
    return data[dim0 * FLayout<dims...>::template compute_slice_offset(idx...)];
  }

  /*
    Zero the entries of the slice
  */
  A2D_INLINE_FUNCTION void zero() {
    for (index_t i = 0; i < FLayout<dims...>::static_size; i++) {
      data[dim0 * i] = 0.0;
    }
  }

  /**
   * Get the underlying data pointer
   */
  A2D_INLINE_FUNCTION T* get_data() { return data; }

 private:
  T* data;
  index_t dim0;
};

/**
 * @brief A partially-specialized MultiArraySlice for CLayout
 */
template <typename T, index_t... dims>
class MultiArraySlice<T, CLayout<dims...>> {
 public:
  using Layout = CLayout<dims...>;

  template <class IdxType>
  A2D_INLINE_FUNCTION MultiArraySlice(
      const MultiArray<T, CLayout<dims...>>& array, const IdxType idx)
      : dim0(array.layout.dim0) {
    data = &array.data[array.layout.get_size(idx)];
  }

  /*
    Constant access to array elements, note: defining member function as const
    gives it better flexibility
  */
  template <class... IdxType>
  A2D_INLINE_FUNCTION T& operator()(IdxType... idx) const {
    return data[CLayout<dims...>::template compute_slice_offset(idx...)];
  }

  /*
    Zero the entries of the slice
  */
  A2D_INLINE_FUNCTION void zero() {
    for (index_t i = 0; i < CLayout<dims...>::static_size; i++) {
      data[i] = 0.0;
    }
  }

  /**
   * Get the underlying data pointer
   */
  A2D_INLINE_FUNCTION T* get_data() { return data; }

 private:
  T* data;
  index_t dim0;
};

/**
 * @brief Create a MultiArraySlice
 */
template <typename T, typename IdxType, class Layout>
A2D_INLINE_FUNCTION MultiArraySlice<T, Layout> MakeSlice(
    const MultiArray<T, Layout>& array, IdxType idx) {
  return MultiArraySlice<T, Layout>(array, idx);
}

}  // namespace A2D

#endif  // A2D_MULTI_ARRAY_H
