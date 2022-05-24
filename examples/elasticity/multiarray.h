#ifndef MULTI_ARRAY_H
#define MULTI_ARRAY_H

#include <type_traits>
#include <cstddef>

namespace A2D {

  /*
    Unpack the argument list to get the extent of indices
  */
  template <int r, std::size_t... dims>
  struct ___get_extent;

  template <std::size_t dim0, std::size_t... dims>
  struct ___get_extent<0, dim0, dims...> {
    static const std::size_t extent = dim0;
  };

  template <int r, std::size_t dim0, std::size_t... dims>
  struct ___get_extent<r, dim0, dims...> {
    static const std::size_t extent = ___get_extent<r-1, dims...>::extent;
  };

  /*
    Fortran ordering

    Given an entry (i, j, k) in a multi-dimensional array, the fortran ordering
    stores entry A(i, j, k) in the location i + dim1 * (j + dim2 * k).
  */
  template <std::size_t... dims>
  class FLayout {
  public:
    FLayout( const std::size_t dim1 ) : dim1(dim1), static_extents{dims...} {}
    const std::size_t dim1;
    static const std::size_t rank = sizeof... (dims) + 1;

    static const std::size_t get_rank(){
      return rank;
    }

    const std::size_t get_extent( std::size_t index ) const {
      if (index == 0){
        return dim1;
      }
      return static_extents[index-1];
    }

    template <class Idx, class... IdxType>
    const std::size_t compute_index( Idx i1, IdxType... idx ) const {
      return i1 + dim1 * __compute_index<0>(idx...);
    }

  private:
    std::size_t static_extents[rank-1];

    template<int r, class Idx, class... IdxType>
    static const std::size_t __compute_index( Idx i, IdxType... idx ){
      return i + ___get_extent<r, dims...>::extent * __compute_index<r+1>(idx...);
    }

    template<int r, class Idx>
    static const std::size_t __compute_index( Idx i ){
      return i;
    }
  };

  /*
    C ordering

    Given an entry (i, j, k) in a multi-dimensional array, the c ordering
    stores entry A(i, j, k) in the location (i * dim2 + j) * dim3 + k.
  */
  template <std::size_t... dims>
  class CLayout {
  public:
    CLayout( const std::size_t dim1 ) : dim1(dim1), static_extents{dims...} {}
    const std::size_t dim1;
    static const std::size_t rank = sizeof... (dims) + 1;

    static const std::size_t get_rank(){
      return rank;
    }

    const std::size_t get_extent( std::size_t index ) const {
      if (index == 0){
        return dim1;
      }
      return static_extents[index-1];
    }

    template <class Idx, class... IdxType>
    const std::size_t compute_index( Idx i1, IdxType... idx ) const {
      return __compute_index<0>(i1, idx...);
    }

  private:
    std::size_t static_extents[rank-1];

    template<int r, class Idx, class... IdxType>
    static const std::size_t __compute_index( const std::size_t index, Idx i, IdxType... idx ){
      return __compute_index<r+1>(index * ___get_extent<r, dims...>::extent + i, idx...);
    }

    template<int r, class Idx>
    static const std::size_t __compute_index( const std::size_t index, Idx i ){
      return index * ___get_extent<r, dims...>::extent + i;
    }
  };

  /*
    A multi-dimensional array
  */
  template <typename T, class Layout>
  class MultiArray {
  public:
    MultiArray( Layout& layout, T* data ) : layout(layout), data(data) {}

    Layout& layout;
    T* data;

    template<class... IdxType>
    T& operator()( IdxType... idx ){
      return data[layout.compute_index(idx...)];
    }

    const std::size_t rank() const {
      return layout.get_rank();
    }

    const std::size_t extent( std::size_t index ) const {
      return layout.get_extent(index);
    }
  };

} // namepsace A2D

#endif // MULTI_ARRAY_H
