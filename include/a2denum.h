#ifndef A2D_ENUM_H
#define A2D_ENUM_H
namespace A2D {

/**
 * @brief Whether a variable should be automatically differentiated or not
 */
enum class ADiffType { PASSIVE, ACTIVE };

/**
 * @brief Whether we compute first or second order derivatives
 */
enum class ADorder { ZERO, FIRST, SECOND };

/**
 * @brief The class of object
 */
enum class ADObjType { SCALAR, VECTOR, MATRIX, SYMMAT };

/**
 * @brief Is the matrix normal (not transposed) or transposed
 */
enum class MatOp { NORMAL, TRANSPOSE };

/**
 * @brief The symmetry type of the matrix
 */
enum class MatSymType { NORMAL, SYMMETRIC };

/**
 * @brief Specify bvalue, pvalue or hvalue
 */
enum class ADseed { b, p, h };

/**
 * @brief compile-time conditional for type values
 */
template <bool B, class TrueType, class FalseType>
struct conditional_type {
  using type = TrueType;
};

template <class TrueType, class FalseType>
struct conditional_type<false, TrueType, FalseType> {
  using type = FalseType;
};

/**
 * @brief compile-time conditional for non-type values, works like
 * std::conditional, but user gets non-type values instead.
 *
 * @tparam T type of the non-type value
 * @tparam B condition, true or false
 * @tparam TrueVal the value user gets if B is true
 * @tparam FalseVal the value user gets if B is false
 */
template <class T, bool B, T TrueVal, T FalseVal>
struct conditional_value {
  static constexpr T value = FalseVal;
};

template <class T, T TrueVal, T FalseVal>
struct conditional_value<T, true, TrueVal, FalseVal> {
  static constexpr T value = TrueVal;
};

}  // namespace A2D
#endif  // A2D_ENUM_H