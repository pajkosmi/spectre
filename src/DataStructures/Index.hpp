// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class template Index.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <ostream>

#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/IsInteger.hpp"

namespace PUP {
class er;
}  // namespace PUP

/// \ingroup DataStructuresGroup
/// An integer multi-index.
///
/// \tparam Dim the number of integers in the Index.
template <size_t Dim>
class Index {
 public:
  /// Construct with each element set to the same value.
  explicit Index(const size_t i0 = std::numeric_limits<size_t>::max())
      : indices_(make_array<Dim>(i0)) {}

  /// Construct specifying value in each dimension
  template <typename... I, Requires<(sizeof...(I) > 1)> = nullptr>
  explicit Index(I... i) : indices_(make_array(static_cast<size_t>(i)...)) {
    static_assert(std::conjunction_v<tt::is_integer<I>...>,
                  "You must pass in a set of size_t's to Index.");
    static_assert(Dim == sizeof...(I),
                  "The number of indices given to Index must be the same as "
                  "the dimensionality of the Index.");
  }

  explicit Index(std::array<size_t, Dim> i) : indices_(std::move(i)) {}

  size_t operator[](const size_t d) const { return gsl::at(indices_, d); }
  size_t& operator[](const size_t d) { return gsl::at(indices_, d); }

  typename std::array<size_t, Dim>::iterator begin() {
    return indices_.begin();
  }
  typename std::array<size_t, Dim>::const_iterator begin() const {
    return indices_.begin();
  }

  typename std::array<size_t, Dim>::iterator end() { return indices_.end(); }
  typename std::array<size_t, Dim>::const_iterator end() const {
    return indices_.end();
  }

  size_t size() const { return Dim; }

  /// The product of the indices.
  /// If Dim = 0, the product is defined as 1.
  template <int N = Dim, Requires<(N > 0)> = nullptr>
  constexpr size_t product() const {
    return indices_[N - 1] * product<N - 1>();
  }
  /// \cond
  // Specialization for N = 0 to stop recursion
  template <int N = Dim, Requires<(N == 0)> = nullptr>
  constexpr size_t product() const {
    return 1;
  }
  /// \endcond

  /// Return a smaller Index with the d-th element removed.
  ///
  /// \param d the element to remove.
  template <size_t N = Dim, Requires<(N > 0)> = nullptr>
  Index<Dim - 1> slice_away(const size_t d) const {
    ASSERT(d < Dim,
           "Can't slice dimension " << d << " from an Index<" << Dim << ">");
    std::array<size_t, Dim - 1> t{};
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 13
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds="
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 13
    for (size_t i = 0; i < d; ++i) {
      gsl::at(t, i) = gsl::at(indices_, i);
    }
    for (size_t i = d + 1; i < Dim; ++i) {
      gsl::at(t, i - 1) = gsl::at(indices_, i);
    }
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 13
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 13
    return Index<Dim - 1>(t);
  }

  /// \cond
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
  /// \endcond

  template <size_t N>
  friend std::ostream& operator<<(std::ostream& os,  // NOLINT
                                  const Index<N>& i);

  const size_t* data() const { return indices_.data(); }
  size_t* data() { return indices_.data(); }

  const std::array<size_t, Dim>& indices() const { return indices_; }

 private:
  std::array<size_t, Dim> indices_;
};

/// \ingroup DataStructuresGroup
/// Get the collapsed index into a 1D array of the data corresponding to this
/// Index. Note that the first dimension of the Index varies fastest when
/// computing the collapsed index.
template <size_t N>
size_t collapsed_index(const Index<N>& index, const Index<N>& extents);

template <size_t N>
std::ostream& operator<<(std::ostream& os, const Index<N>& i);

/// \cond HIDDEN_SYMBOLS
#ifdef SPECTRE_DEBUG
namespace Index_detail {
template <size_t Dim>
void collapsed_index_check(const Index<Dim>& index, const Index<Dim>& extents) {
  for (size_t d = 0; d < 1; ++d) {
    ASSERT(index[d] < extents[d], "The requested index in the dimension "
                                      << d << " with value " << index[d]
                                      << " exceeds the number of grid "
                                         "points "
                                      << extents[d]);
  }
}
}  // namespace Index_detail
#endif

// the specializations are in the header file so they can be inlined. We use
// specializations to avoid having loops since this computation is very
// straightforward.
template <>
SPECTRE_ALWAYS_INLINE size_t collapsed_index(const Index<0>& /*index*/,
                                             const Index<0>& /*extents*/) {
  return 0;
}

template <>
SPECTRE_ALWAYS_INLINE size_t collapsed_index(const Index<1>& index,
                                             const Index<1>& extents) {
  (void)extents;
#ifdef SPECTRE_DEBUG
  Index_detail::collapsed_index_check(index, extents);
#endif
  return index[0];
}

template <>
SPECTRE_ALWAYS_INLINE size_t collapsed_index(const Index<2>& index,
                                             const Index<2>& extents) {
#ifdef SPECTRE_DEBUG
  Index_detail::collapsed_index_check(index, extents);
#endif
  return index[0] + extents[0] * index[1];
}

template <>
SPECTRE_ALWAYS_INLINE size_t collapsed_index(const Index<3>& index,
                                             const Index<3>& extents) {
#ifdef SPECTRE_DEBUG
  Index_detail::collapsed_index_check(index, extents);
#endif
  return index[0] + extents[0] * (index[1] + extents[1] * index[2]);
}

template <>
SPECTRE_ALWAYS_INLINE size_t collapsed_index(const Index<4>& index,
                                             const Index<4>& extents) {
#ifdef SPECTRE_DEBUG
  Index_detail::collapsed_index_check(index, extents);
#endif
  return index[0] +
         extents[0] *
             (index[1] + extents[1] * (index[2] + extents[2] * index[3]));
}

template <size_t N>
bool operator==(const Index<N>& lhs, const Index<N>& rhs);

template <size_t N>
bool operator!=(const Index<N>& lhs, const Index<N>& rhs);
/// \endcond
