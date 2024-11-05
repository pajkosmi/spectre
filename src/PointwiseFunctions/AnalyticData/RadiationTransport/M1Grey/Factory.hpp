// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/HomogeneousSphere.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"
#include "Utilities/TMPL.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {
/*!
 * \brief Typelist of all analytic data of M1Grey evolution system
 */
template <size_t Dim>
using all_data = tmpl::list<HomogeneousSphereImpl, SphericalGaussian<Dim>>;
}  // namespace RadiationTransport::M1Grey::AnalyticData
