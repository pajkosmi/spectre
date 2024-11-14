// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/HomogeneousSphere.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"
#include "PointwiseFunctions/AnalyticSolutions/RadiationTransport/M1Grey/ConstantM1.hpp"
#include "Utilities/TMPL.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {
/*!
 * \brief Typelist of all analytic data of M1Grey evolution system
 */

template <size_t Dim>
using all_data = tmpl::list<HomogeneousSphereImpl, SphericalGaussian<Dim>,
                            Solutions::ConstantM1>;
}  // namespace RadiationTransport::M1Grey::AnalyticData
