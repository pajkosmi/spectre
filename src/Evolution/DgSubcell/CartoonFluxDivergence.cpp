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
        if (inertial_coords.get(0)[volume_index] == 0.0) {
          // if (abs(inertial_coords.get(0)[volume_index] -
          //         0.5 / one_over_delta / inv_jacobian[volume_index])
          //         < 1.0e-14) {
          dfdx = 3.0 * one_over_delta * inv_jacobian[volume_index] *
                 (boundary_correction[boundary_correction_upper_index] -
                  boundary_correction[boundary_correction_lower_index]);
        } else {
          // if spherically symmetric inner boundary, no flux should leave or
          // enter
          double rescale_delta = 1.0;
          double fac = 1.0;
          // if (abs(inertial_coords.get(0)[volume_index] -
          //         0.5 / one_over_delta / inv_jacobian[volume_index]) <
          //     1.0e-14) {
          //   fac = 0.0;
          //   rescale_delta =
          //       boundary_correction[boundary_correction_upper_index] /
          //       boundary_correction[boundary_correction_lower_index];
          // }

          // std::cout << "boundary_correction " << boundary_correction << "\n";
          // std::cout << "boundary_correction_upper_index " <<
          // boundary_correction_upper_index << "
          // boundary_correction_lower_index " <<
          // boundary_correction_lower_index << "\n"; std::cout <<
          // "inertial_coords " << inertial_coords.get(0) << "\n"; std::cout <<
          // "volume index " << volume_index << "\n"; std::cout <<
          // "inertial_coords left face " << inertial_coords.get(0) -
          //       0.5 / one_over_delta / inv_jacobian[volume_index] << "\n";

          dfdx =
              2.0 *
                  (0.5 *
                   (boundary_correction[boundary_correction_upper_index] +
                    rescale_delta *
                        boundary_correction[boundary_correction_lower_index])) /
                  abs(inertial_coords.get(0)[volume_index]) +
              one_over_delta * inv_jacobian[volume_index] *
                  (boundary_correction[boundary_correction_upper_index] -
                   fac * boundary_correction[boundary_correction_lower_index]);

          // product rule
          double inertial_coord_lower_face =
              inertial_coords.get(0)[volume_index] -
              0.5 / one_over_delta / inv_jacobian[volume_index];

          double inertial_coord_upper_face =
              inertial_coords.get(0)[volume_index] +
              0.5 / one_over_delta / inv_jacobian[volume_index];

          double lower_face_weight = inertial_coord_lower_face *
                                     inertial_coord_lower_face /
                                     (inertial_coords.get(0)[volume_index] *
                                      inertial_coords.get(0)[volume_index]);
          double upper_face_weight = inertial_coord_upper_face *
                                     inertial_coord_upper_face /
                                     (inertial_coords.get(0)[volume_index] *
                                      inertial_coords.get(0)[volume_index]);

          dfdx = one_over_delta * inv_jacobian[volume_index] *
                 (upper_face_weight *
                      boundary_correction[boundary_correction_upper_index] -
                  lower_face_weight *
                      boundary_correction[boundary_correction_lower_index]);
        }
        (*dt_var)[volume_index] += dfdx;
        // (*dt_var)[volume_index] +=
        //     one_over_delta * inv_jacobian[volume_index] *
        //     (boundary_correction[boundary_correction_upper_index] -
        //      boundary_correction[boundary_correction_lower_index]);
      }
    }
  }
}
}  // namespace evolution::dg::subcell
