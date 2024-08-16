// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \ingroup GeneralRelativityGroup
/// Holds functions related to general relativity.
namespace gr {

namespace Tags {
/*!
 * \brief Compute item for spatial metric determinant \f$\gamma\f$ and inverse
 * \f$\gamma^{ij}\f$ in terms of the spatial metric \f$\gamma_{ij}\f$.
 *
 * \details Can be retrieved using `gr::Tags::DetSpatialMetric` and
 * `gr::Tags::InverseSpatialMetric`.
 */
template <typename DataType, size_t SpatialDim, typename Frame>
struct DetAndInverseSpatialMetricCompute
    : ::Tags::Variables<
          tmpl::list<DetSpatialMetric<DataType>,
                     InverseSpatialMetric<DataType, SpatialDim, Frame>>>,
      db::ComputeTag {
  using argument_tags = tmpl::list<SpatialMetric<DataType, SpatialDim, Frame>>;
  using base = ::Tags::Variables<
      tmpl::list<DetSpatialMetric<DataType>,
                 InverseSpatialMetric<DataType, SpatialDim, Frame>>>;
  using return_type = typename base::type;
  static constexpr auto function = static_cast<void (*)(
      const gsl::not_null<return_type*>,
      const Tensor<DataType, tmpl::integral_list<std::int32_t, 1, 1>,
                   tmpl::list<SpatialIndex<SpatialDim, UpLo::Lo, Frame>,
                              SpatialIndex<SpatialDim, UpLo::Lo, Frame>>>&)>(
      &determinant_and_inverse);
};

/*!
 * \brief Compute item to get the square root of the determinant of the spatial
 * metric \f$\sqrt{\gamma}\f$ via `gr::Tags::DetAndInverseSpatialMetric`.
 *
 * \details Can be retrieved using `gr::Tags::SqrtDetSpatialMetric`.
 */
template <typename DataType, size_t SpatialDim, typename Frame>
struct SqrtDetSpatialMetricCompute : SqrtDetSpatialMetric<DataType>,
                                     db::ComputeTag {
  using argument_tags = tmpl::list<DetSpatialMetric<DataType>>;

  using return_type = Scalar<DataType>;

  static void function(const gsl::not_null<Scalar<DataType>*> result,
                       const Scalar<DataType>& det_spatial_metric) {
    get(*result) = sqrt(get(det_spatial_metric));
  }

using base = SqrtDetSpatialMetric<DataType>;
};
}  // namespace Tags
}  // namespace gr
