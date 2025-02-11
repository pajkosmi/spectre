// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/DgSubcell/CartoonFluxDivergence.hpp"

#include <cstddef>

#include <iostream>
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

namespace evolution::dg::subcell {
void add_cartesian_flux_divergence(const gsl::not_null<DataVector*> dt_var,
                                   const double one_over_delta,
                                   const DataVector& inv_jacobian,
                                   const DataVector& boundary_correction,
                                   const Index<1>& subcell_extents,
                                   const size_t dimension) {
  (void)dimension;
  ASSERT(dimension == 0, "dimension must be 0 but is " << dimension);
  for (size_t i = 0; i < subcell_extents[0]; ++i) {
    (*dt_var)[i] += one_over_delta * inv_jacobian[i] *
                    (boundary_correction[i + 1] - boundary_correction[i]);
  }
}

void add_cartesian_flux_divergence(const gsl::not_null<DataVector*> dt_var,
                                   const double one_over_delta,
                                   const DataVector& inv_jacobian,
                                   const DataVector& boundary_correction,
                                   const Index<2>& subcell_extents,
                                   const size_t dimension) {
  ASSERT(dimension == 0 or dimension == 1,
         "dimension must be 0 or 1 but is " << dimension);
  Index<2> subcell_face_extents = subcell_extents;
  ++subcell_face_extents[dimension];
  for (size_t j = 0; j < subcell_extents[1]; ++j) {
    for (size_t i = 0; i < subcell_extents[0]; ++i) {
      Index<2> index(i, j);
      const size_t volume_index = collapsed_index(index, subcell_extents);
      const size_t boundary_correction_lower_index =
          collapsed_index(index, subcell_face_extents);
      ++index[dimension];
      const size_t boundary_correction_upper_index =
          collapsed_index(index, subcell_face_extents);
      (*dt_var)[volume_index] +=
          one_over_delta * inv_jacobian[volume_index] *
          (boundary_correction[boundary_correction_upper_index] -
           boundary_correction[boundary_correction_lower_index]);
    }
  }
}

void add_cartesian_flux_divergence(
    const gsl::not_null<DataVector*> dt_var, const double one_over_delta,
    const DataVector& inv_jacobian, const DataVector& boundary_correction,
    const Index<3>& subcell_extents, const size_t dimension,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords) {
  ASSERT(dimension == 0 or dimension == 1 or dimension == 2,
         "dimension must be 0, 1, or 2 but is " << dimension);
  Index<3> subcell_face_extents = subcell_extents;
  ++subcell_face_extents[dimension];

  double dfdx = 0.0;

  for (size_t k = 0; k < subcell_extents[2]; ++k) {
    for (size_t j = 0; j < subcell_extents[1]; ++j) {
      for (size_t i = 0; i < subcell_extents[0]; ++i) {
        Index<3> index(i, j, k);

        // spherical symmetry change
        const size_t volume_index = collapsed_index(index, subcell_extents);
        const size_t boundary_correction_lower_index =
            collapsed_index(index, subcell_face_extents);
        ++index[dimension];
        const size_t boundary_correction_upper_index =
            collapsed_index(index, subcell_face_extents);
        // We are evaluating div F from d_t U + div F = S.
        // For 1D cartoon, div F = dF/dx + 2 * F / x (see SpECTRE book).
        // Because at the origin, the flux is zero, we apply L'Hospital's rule
        // to the 2 * F / x term.  Differentiating top and bottom terms and
        // collecting: div F = 3 dF/dx.
        if (inertial_coords.get(0)[volume_index] == 0.0) {
          dfdx = 3.0 * one_over_delta * inv_jacobian[volume_index] *
                 (boundary_correction[boundary_correction_upper_index] -
                  boundary_correction[boundary_correction_lower_index]);
        } else {
          double inertial_coord_lower_face =
              inertial_coords.get(0)[volume_index] -
              0.5 / one_over_delta / inv_jacobian[volume_index];

          double inertial_coord_upper_face =
              inertial_coords.get(0)[volume_index] +
              0.5 / one_over_delta / inv_jacobian[volume_index];

          // (x_low / x_center)^2
          double lower_face_weight = inertial_coord_lower_face *
                                     inertial_coord_lower_face /
                                     (inertial_coords.get(0)[volume_index] *
                                      inertial_coords.get(0)[volume_index]);
          // (x_hi / x_center)^2
          double upper_face_weight = inertial_coord_upper_face *
                                     inertial_coord_upper_face /
                                     (inertial_coords.get(0)[volume_index] *
                                      inertial_coords.get(0)[volume_index]);

          // If not at the origin div F is simply = df/dx  = 1 / (delta x) *
          // ((x_hi / x_center)^2 * Flux_upper - (x_low / x_center)^2 *
          // Flux_lower).
          dfdx = one_over_delta * inv_jacobian[volume_index] *
                 (upper_face_weight *
                      boundary_correction[boundary_correction_upper_index] -
                  lower_face_weight *
                      boundary_correction[boundary_correction_lower_index]);
        }
        (*dt_var)[volume_index] += dfdx;
      }
    }
  }
}
}  // namespace evolution::dg::subcell
