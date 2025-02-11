// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Index;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

namespace evolution::dg::subcell {
/// @{
/*!
 * \brief Compute and add the flux divergences on a Cartesian mesh for the
 * spherically symmetric (1D) cartoon method.
 *
 * For more details, see the Cartoon method documentation in the SpECTRE
 * book.
 */
void add_cartesian_flux_divergence(gsl::not_null<DataVector*> dt_var,
                                   double one_over_delta,
                                   const DataVector& inv_jacobian,
                                   const DataVector& boundary_correction,
                                   const Index<1>& subcell_extents,
                                   size_t dimension);

void add_cartesian_flux_divergence(gsl::not_null<DataVector*> dt_var,
                                   double one_over_delta,
                                   const DataVector& inv_jacobian,
                                   const DataVector& boundary_correction,
                                   const Index<2>& subcell_extents,
                                   size_t dimension);

void add_cartesian_flux_divergence(
    gsl::not_null<DataVector*> dt_var, double one_over_delta,
    const DataVector& inv_jacobian, const DataVector& boundary_correction,
    const Index<3>& subcell_extents, size_t dimension,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords);
/// @}
}  // namespace evolution::dg::subcell
