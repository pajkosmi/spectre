// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Transpose.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/Tags.hpp"
#include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.hpp"
// #include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/StdArrayHelpers.hpp"

#include <iostream>

// namespace cartoon_namespace {
// template <typename ResultTags, size_t Dim, typename DerivativeFrame>
// void cartoon_tensor_derivatives() {
//   std::cout << "test function call \n";
// }

// }  // namespace cartoon_namespace

namespace partial_derivatives_detail {
template <size_t Dim, typename VariableTags, typename DerivativeTags>
struct LogicalImpl;

// This routine has been optimized to perform really well. The following
// describes what optimizations were made.
//
// - The `partial_derivatives` functions below have an overload where the
//   logical derivatives may be passed in instead of being computed. In the
//   overloads where the logical derivatives are not passed in they must be
//   computed. However, it is more efficient to allocate the memory for the
//   logical partial derivatives with respect to each coordinate at once. This
//   requires the `partial_derivatives_impl` to accept raw pointers to doubles
//   for the logical derivatives so it can be used for all overloads.
//
// - The resultant Variables `du` is a not_null pointer so that mutating compute
//   items can be supported.
//
// - The storage indices into the inverse Jacobian are precomputed to avoid
//   having to recompute them for each tensor component of `u`.
//
// - The DataVectors lhs and logical_du are non-owning DataVectors to be able to
//   plug into the optimized expression templates. This requires a `const_cast`
//   even though we will never change the `double*`.
//
// - Loop over every Tensor component in the variables by incrementing a raw
//   pointer to the contiguous data (vs. looping over each Tensor in the
//   variables with a tmpl::for_each then iterating over each component of this
//   Tensor).
//
// - We factor out the `logical_deriv_index == 0` case so that we do not need to
//   zero the memory in `du` before the computation.

// spacetime_metric component index = 0 tensor index = (0,0)
// spacetime_metric component index = 1 tensor index = (1,0)
// spacetime_metric component index = 2 tensor index = (2,0)
// spacetime_metric component index = 3 tensor index = (3,0)
// spacetime_metric component index = 4 tensor index = (1,1)
// spacetime_metric component index = 5 tensor index = (2,1)
// spacetime_metric component index = 6 tensor index = (3,1)
// spacetime_metric component index = 7 tensor index = (2,2)
// spacetime_metric component index = 8 tensor index = (3,2)
// spacetime_metric component index = 9 tensor index = (3,3)

template <typename ResultTags, size_t Dim, typename DerivativeFrame>
void lie_drag_covariant_rank_2_tensor(
    double& dfdx,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const gsl::span<const double>& volume_vars, size_t volume_index,
    const size_t component_index, const size_t num_grid_points,
    const size_t deriv_index) {
  // becomes 10 if component_index >= 10, b/c of floor division
  int shift_index_by_10 = 10 * (component_index / 10);

  if (UNLIKELY(shift_index_by_10 > 10)) {
    ERROR("shift_index_by_10 factor is too high (" << shift_index_by_10
                                                   << "). Should be <= 10.");
  }

  size_t component_index_shift = component_index - shift_index_by_10;

  auto index_shift = [&num_grid_points, &volume_index,
                      &volume_vars](size_t shifted_index) {
    return volume_vars[shifted_index * num_grid_points + volume_index];
  };

  auto var = [&volume_vars, &num_grid_points, &volume_index, &index_shift,
              &shift_index_by_10](size_t shifted_index) {
    return index_shift(shifted_index + shift_index_by_10);
  };

  if (deriv_index == 1) {
    // y deriv
    if (component_index_shift == 0) {
      // dy_gtt = 0
      dfdx = 0.0;
    } else if (component_index_shift == 1) {
      // d_y g_tx = -g_02 = -g(2)
      dfdx = -var(2);
    } else if (component_index_shift == 2) {
      // d_y g_ty = g_01 = g(1)
      dfdx = var(1);
    } else if (component_index_shift == 3) {
      // dy g_tz = 0
      dfdx = 0.0;
    } else if (component_index_shift == 4) {
      // d_y g_xx = -2 g_12 = -2 g(5)
      dfdx = -2.0 * var(5);
    } else if (component_index_shift == 5) {
      // d_y g_xy = g_11 - g_22 = g(4) - g(7)
      dfdx = var(4) - var(7);
    } else if (component_index_shift == 6) {
      // d_y g_xz = -g_23 = -g(8)
      dfdx = -var(8);
    } else if (component_index_shift == 7) {
      // d_y g_yy = 2 * g_12 = 2 g(5)
      dfdx = 2.0 * var(5);
    } else if (component_index_shift == 8) {
      // d_y g_yz = g_13 = g(6)
      dfdx = var(6);
    } else if (component_index_shift == 9) {
      // d_y g_zz = 0
      dfdx = 0;
    } else {
      ERROR("NOT 0-9 COMPONENT INDEX. COMPONENT = " << component_index_shift
                                                    << ".");
    }
  } else {
    // z deriv
    if (component_index_shift == 0) {
      // dz g_tt = 0
      dfdx = 0.0;
    } else if (component_index_shift == 1) {
      // d_z g_tx = -g_03 = -g(3)
      dfdx = -var(3);
    } else if (component_index_shift == 2) {
      // d_z g_ty = 0
      dfdx = 0.0;
    } else if (component_index_shift == 3) {
      // dz g_tz = g01 = g(1)
      dfdx = var(1);
    } else if (component_index_shift == 4) {
      // d_z g_xx = -2 g_13 = -2 * g(6)
      dfdx = -2.0 * var(6);
    } else if (component_index_shift == 5) {
      // d_z g_xy =  -g_23 = -g(8)
      dfdx = -var(8);
    } else if (component_index_shift == 6) {
      // d_z g_xz = g_11 - g_33 = g(4) - g(9)
      dfdx = var(4) - var(9);
    } else if (component_index_shift == 7) {
      // d_z g_yy = 0
      dfdx = 0.0;
    } else if (component_index_shift == 8) {
      // d_z g_yz = g_12 = g(5)
      dfdx = var(5);
    } else if (component_index_shift == 9) {
      // d_z g_zz = 2 * g_13 = 2 * g(6)
      dfdx = 2.0 * var(6);
    } else {
      ERROR("NOT 0-9 COMPONENT INDEX. COMPONENT = " << component_index_shift
                                                    << ".");
    }
  }
  // scale by -1
  dfdx *= -1.0;  /// abs(inertial_coords.get(0)[volume_index]);
}

// Phi component index = 0 tensor index = (0,0,0)
// Phi component index = 1 tensor index = (1,0,0)
// Phi component index = 2 tensor index = (2,0,0)
// Phi component index = 3 tensor index = (0,1,0)
// Phi component index = 4 tensor index = (1,1,0)
// Phi component index = 5 tensor index = (2,1,0)
// Phi component index = 6 tensor index = (0,2,0)
// Phi component index = 7 tensor index = (1,2,0)
// Phi component index = 8 tensor index = (2,2,0)
// Phi component index = 9 tensor index = (0,3,0)
// Phi component index = 10 tensor index = (1,3,0)
// Phi component index = 11 tensor index = (2,3,0)
// Phi component index = 12 tensor index = (0,1,1)
// Phi component index = 13 tensor index = (1,1,1)
// Phi component index = 14 tensor index = (2,1,1)
// Phi component index = 15 tensor index = (0,2,1)
// Phi component index = 16 tensor index = (1,2,1)
// Phi component index = 17 tensor index = (2,2,1)
// Phi component index = 18 tensor index = (0,3,1)
// Phi component index = 19 tensor index = (1,3,1)
// Phi component index = 20 tensor index = (2,3,1)
// Phi component index = 21 tensor index = (0,2,2)
// Phi component index = 22 tensor index = (1,2,2)
// Phi component index = 23 tensor index = (2,2,2)
// Phi component index = 24 tensor index = (0,3,2)
// Phi component index = 25 tensor index = (1,3,2)
// Phi component index = 26 tensor index = (2,3,2)
// Phi component index = 27 tensor index = (0,3,3)
// Phi component index = 28 tensor index = (1,3,3)
// Phi component index = 29 tensor index = (2,3,3)
template <typename ResultTags, size_t Dim, typename DerivativeFrame>
void lie_drag_covariant_rank_3_tensor(
    double& dfdx,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const gsl::span<const double>& volume_vars, size_t volume_index,
    const size_t component_index, const size_t num_grid_points,
    const size_t deriv_index) {
  auto index_shift = [&num_grid_points, &volume_index,
                      &volume_vars](size_t shifted_index) {
    return volume_vars[shifted_index * num_grid_points + volume_index];
  };

  int shift_index_by_20 = 20;

  auto var = [&volume_vars, &num_grid_points, &volume_index, &index_shift,
              &shift_index_by_20](size_t shifted_index) {
    return index_shift(shifted_index + shift_index_by_20);
  };

  size_t component_index_shift = component_index - shift_index_by_20;

  if (UNLIKELY(component_index_shift > 49)) {
    ERROR("component_index_shift factor is too high ("
          << component_index_shift << "). Should be <= 50.");
  }

  if (UNLIKELY(component_index_shift < 0)) {
    ERROR("component_index_shift factor is too low (" << component_index_shift
                                                      << "). Should be > 0.");
  }

  if (deriv_index == 1) {
    // y deriv
    // if (component_index % 3 == 0) {
    //   // All Phi_xab components. Eqn (33)
    //   // dfdx =
    // } else if ((component_index - 1) % 3 == 0) {
    //   // All Phi_yab components. Eqn (38)
    //   // dfdx =

    // } else if ((component_index - 2) % 3 == 0) {
    //   // All Phi_zab components. Eqn (43)
    //   // dfdx =
    // }
    if (component_index_shift == 0) {
      // dy Phi_xtt = Phi_y00 = p(1) [Eqn (33)]
      dfdx = var(1);
    } else if (component_index_shift == 1) {
      // dy Phi_ytt = -Phi_x00 = -p(0) [Eqn (38)]
      dfdx = -var(0);
    } else if (component_index_shift == 2) {
      // dy Phi_ztt = 0 [Eqn(43)]
      dfdx = 0.0;
    } else if (component_index_shift == 3) {
      // dy Phi_xxt = Phi_y01 + Phi_x02 = p(4) + p(6)
      dfdx = var(4) + var(6);
    } else if (component_index_shift == 4) {
      // dy Phi_yxt = -Phi_x01 + Phi_y02 = -p(3) + p(7)
      dfdx = -var(3) + var(7);
    } else if (component_index_shift == 5) {
      // dy Phi_zxt = Phi_z02 = p(8)
      dfdx = var(8);
    } else if (component_index_shift == 6) {
      // dy Phi_xyt = Phi_y02 - Phi_x01 = p(7) - p(3)
      dfdx = var(7) - var(3);
    } else if (component_index_shift == 7) {
      // dy Phi_yyt = -Phi_x02 - Phi_y01 = -p(6) - p(4)
      dfdx = -var(6) - var(4);
    } else if (component_index_shift == 8) {
      // dy Phi_zyt = -Phi_z01 = -p(5)
      dfdx = -var(5);
    } else if (component_index_shift == 9) {
      // dy Phi_xzt = Phi_y03 = p(10)
      dfdx = var(10);
    } else if (component_index_shift == 10) {
      // dy Phi_yzt = -Phi_x03 = -p(9)
      dfdx = -var(9);
    } else if (component_index_shift == 11) {
      // dy Phi_zzt = 0
      dfdx = 0.0;
    } else if (component_index_shift == 12) {
      // dy Phi_xxx = Phi_y11 + 2 * Phi_x21 = p(13) + 2 * p(15)
      dfdx = var(13) + 2.0 * var(15);
    } else if (component_index_shift == 13) {
      // dy Phi_yxx = -Phi_x11 + 2 * Phi_y21 = -p(12) + 2 * p(16)
      dfdx = -var(12) + 2.0 * var(16);
    } else if (component_index_shift == 14) {
      // dy Phi_zxx = 2 * Phi_z21 = 2 * p(17)
      dfdx = 2.0 * var(17);
    } else if (component_index_shift == 15) {
      // dy Phi_xyx = Phi_y12 + Phi_x22 - Phi_x11 = p(16) + p(21) - p(12)
      dfdx = var(16) + var(21) - var(12);
    } else if (component_index_shift == 16) {
      // dy Phi_yyx = -Phi_x12 + Phi_y22 - Phi_y11 = -p(15) + p(22) - p(13)
      dfdx = -var(15) + var(22) - var(13);
    } else if (component_index_shift == 17) {
      // dy Phi_zyx = Phi_z22 - Phi_z11 = p(23) - p(14)
      dfdx = var(23) - var(14);
    } else if (component_index_shift == 18) {
      // dy Phi_xzx = Phi_y13 + Phi_x23 = p(19) + p(24)
      dfdx = var(19) + var(24);
    } else if (component_index_shift == 19) {
      // dy Phi_yzx = -Phi_x13 + Phi_y23 = -p(18) + p(25)
      dfdx = -var(18) + var(25);
    } else if (component_index_shift == 20) {
      // dy Phi_zzx = Phi_z23 = p(26)
      dfdx = var(26);
    } else if (component_index_shift == 21) {
      // dy Phi_xyy = Phi_y22 - 2Phi_x12 = p(22) - 2 p(15)
      dfdx = var(22) - 2.0 * var(15);
    } else if (component_index_shift == 22) {
      // dy Phi_yyy = -Phi_x22 - 2Phi_y12 = -p(21) - 2 p(16)
      dfdx = -var(21) - 2.0 * var(16);
    } else if (component_index_shift == 23) {
      // dy Phi_zyy = -2Phi_z12 = -2 p(17)
      dfdx = -2.0 * var(17);
    } else if (component_index_shift == 24) {
      // dy Phi_xzy = Phi_y23 - Phi_x13 = p(25) - p(18)
      dfdx = var(25) - var(18);
    } else if (component_index_shift == 25) {
      // dy Phi_yzy = -Phi_x23 - Phi_y13 = -p(24) - p(19)
      dfdx = -var(24) - var(19);
    } else if (component_index_shift == 26) {
      // dy Phi_zzy = -Phi_z13 = -p(20)
      dfdx = -var(20);
    } else if (component_index_shift == 27) {
      // dy Phi_xzz = Phi_y33 = p(28)
      dfdx = var(28);
    } else if (component_index_shift == 28) {
      // dy Phi_yzz = -Phi_x33 = -p(27)
      dfdx = -var(27);
    } else if (component_index_shift == 29) {
      // dy Phi_zzz = 0
      dfdx = 0.0;
    } else {
      // One of the above should be selected
      ERROR("D/DY: NOT 0-29 COMPONENT INDEX. COMPONENT = "
            << component_index_shift << ".");
    }
  } else {
    // z deriv
    // if (component_index % 3 == 0) {
    //   // All Phi_xab components. Eqn (49)
    //   // dfdx =
    // } else if ((component_index - 1) % 3 == 0) {
    //   // All Phi_yab components. Eqn (54)
    //   // dfdx =

    // } else if ((component_index - 2) % 3 == 0) {
    //   // All Phi_zab components. Eqn (59)
    //   // dfdx =
    // }
    if (component_index_shift == 0) {
      // dz Phi_xtt = Phi_z00 = p(2) [Eqn (49)]
      dfdx = var(2);
    } else if (component_index_shift == 1) {
      // dz Phi_ytt = 0 [Eqn (38)]
      dfdx = 0.0;
    } else if (component_index_shift == 2) {
      // dz Phi_ztt = -Phi_x00 = -p(0) [Eqn(43)]
      dfdx = -var(0);
    } else if (component_index_shift == 3) {
      // dz Phi_xxt = Phi_z01 + Phi_x03 = p(5) + p(9)
      dfdx = var(5) + var(9);
    } else if (component_index_shift == 4) {
      // dz Phi_yxt = Phi_y03 = p(10)
      dfdx = var(10);
    } else if (component_index_shift == 5) {
      // dz Phi_zxt = -Phi_x01 + Phi_z03 = -p(3) + p(11)
      dfdx = -var(3) + var(11);
    } else if (component_index_shift == 6) {
      // dz Phi_xyt = Phi_z02 = p(8)
      dfdx = var(8);
    } else if (component_index_shift == 7) {
      // dz Phi_yyt = 0
      dfdx = 0.0;
    } else if (component_index_shift == 8) {
      // dz Phi_zyt = -Phi_x02 = -p(6)
      dfdx = -var(6);
    } else if (component_index_shift == 9) {
      // dz Phi_xzt = Phi_z03 - Phi_x01 = p(11) - p(3)
      dfdx = var(11) - var(3);
    } else if (component_index_shift == 10) {
      // dz Phi_yzt = -Phi_y01 = -p(4)
      dfdx = -var(4);
    } else if (component_index_shift == 11) {
      // dz Phi_zzt = -Phi_x03 - Phi_z01 = -p(9) - p(5)
      dfdx = -var(9) - var(5);
    } else if (component_index_shift == 12) {
      // dz Phi_xxx = Phi_z11 + 2 * Phi_x31 = p(14) + 2 * p(18)
      dfdx = var(14) + 2.0 * var(18);
    } else if (component_index_shift == 13) {
      // dz Phi_yxx = 2 * Phi_y31 = 2 * p(19)
      dfdx = 2.0 * var(19);
    } else if (component_index_shift == 14) {
      // dz Phi_zxx = -Phi_x11 + 2 * Phi_z31 = -p(12) + 2 * p(20)
      dfdx = -var(12) + 2.0 * var(20);
    } else if (component_index_shift == 15) {
      // dz Phi_xyx = Phi_z12 + Phi_x32 = p(17) + p(24)
      dfdx = var(17) + var(24);
    } else if (component_index_shift == 16) {
      // dz Phi_yyx = Phi_y32 = p(25)
      dfdx = var(25);
    } else if (component_index_shift == 17) {
      // dz Phi_zyx = -Phi_x12 + Phi_z32 = -p(15) + p(26)
      dfdx = -var(15) + var(26);
    } else if (component_index_shift == 18) {
      // dz Phi_xzx = Phi_z13 + Phi_x33 - Phi_x11 = p(20) + p(27) - p(12)
      dfdx = var(20) + var(27) - var(12);
    } else if (component_index_shift == 19) {
      // dz Phi_yzx = Phi_y33 - Phi_y11 = p(28) - p(13)
      dfdx = var(28) - var(13);
    } else if (component_index_shift == 20) {
      // dz Phi_zzx = -x13 + z33 - z11 = -p(18) + p(29) - p(14)
      dfdx = -var(18) + var(29) - var(14);
    } else if (component_index_shift == 21) {
      // dz Phi_xyy = z22 = p(23)
      dfdx = var(23);
    } else if (component_index_shift == 22) {
      // dz Phi_yyy = 0
      dfdx = 0;
    } else if (component_index_shift == 23) {
      // dz Phi_zyy = -x22 = -p(21)
      dfdx = -var(21);
    } else if (component_index_shift == 24) {
      // dz Phi_xzy = z23 - x21 = p(26) - p(15)
      dfdx = var(26) - var(15);
    } else if (component_index_shift == 25) {
      // dz Phi_yzy = -y21 = -p(16)
      dfdx = -var(16);
    } else if (component_index_shift == 26) {
      // dz Phi_zzy = -x23 - z21 = -p(24) - p(17)
      dfdx = -var(24) - var(17);
    } else if (component_index_shift == 27) {
      // dz Phi_xzz = z33 - 2 * x13 = p(29) - 2.0 * var(18)
      dfdx = var(29) - 2.0 * var(18);
    } else if (component_index_shift == 28) {
      // dz Phi_yzz = -2 * y13 = -2 * p(19)
      dfdx = -2.0 * var(19);
    } else if (component_index_shift == 29) {
      // dz Phi_zzz = -x33 - 2.0 * z13 = -p(27) - 2.0 * p(20)
      dfdx = -var(27) - 2.0 * var(20);
    } else {
      // One of the above should be selected
      ERROR("D/DZ: NOT 0-29 COMPONENT INDEX. COMPONENT = "
            << component_index_shift << ".");
    }
  }
  // scale by -1/x
  // dfdx *= 1.0 / abs(inertial_coords.get(0)[volume_index]);
}

// need template parameters or else duplicate definitions
//
template <typename ResultTags, size_t Dim, typename DerivativeFrame>
void cartoon_tensor_derivatives(
    double& dfdx, const size_t deriv_index, const size_t volume_index,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian,
    const DataVector logical_du,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const gsl::span<const double>& volume_vars, const size_t component_index,
    const size_t num_grid_points) {
  if (deriv_index == 0) {
    // x derivatives
    // scale lhs (logical_du) by inverse jacobian--giving pdu
    // df/dx (numerical derivative)
    dfdx = inverse_jacobian.get(0, 0)[volume_index] * logical_du[volume_index];
  } else {
    // y & z derivatives

    // Overleaf definition for Lie
    // dragging a covariant rank 2 tensor (g_ab and Pi_ab) or a
    // covariant rank 3 tensor (Phi_iab)

    // Mike TODO:
    // i) load

    if (component_index < 20) {
      // g_ab & Pi_ab
      lie_drag_covariant_rank_2_tensor<ResultTags, Dim, DerivativeFrame>(
          dfdx, inertial_coords, volume_vars, volume_index, component_index,
          num_grid_points, deriv_index);
    } else {
      // Phi_iab
      lie_drag_covariant_rank_3_tensor<ResultTags, Dim, DerivativeFrame>(
          dfdx, inertial_coords, volume_vars, volume_index, component_index,
          num_grid_points, deriv_index);
    }

    // scale by -1 / x
    dfdx *= -1.0 / inertial_coords.get(0)[volume_index];
    // dfdx = -1.0 / abs(inertial_coords.get(0)[volume_index]) *
    //        (volume_vars[component_index * num_grid_points + volume_index]);
  }
}

template <typename ResultTags, size_t Dim, typename DerivativeFrame>
void partial_derivatives_cartoon(
    const gsl::not_null<Variables<ResultTags>*> du,
    const std::array<const double*, Dim>& logical_partial_derivatives_of_u,
    const size_t number_of_independent_components,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const gsl::span<const double>& volume_vars,
    const Index<3>& subcell_extents) {
  // pdu points to du
  double* pdu = du->data();

  const size_t num_grid_points = du->number_of_grid_points();
  DataVector lhs{};
  DataVector logical_du{};

  using tag_list =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>,
                 gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>;

  const Variables<tag_list> metric_quantities{
      const_cast<double*>(volume_vars.data()), volume_vars.size()};

  tnsr::iaa<DataVector, Dim, Frame::Inertial> deriv_spacetime_metric;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> deriv_pi_evolution;
  tnsr::ijaa<DataVector, Dim, Frame::Inertial> deriv_phi_evolution;

  ::fd::general_cartoon_deriv(
      deriv_spacetime_metric,
      get<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>>(
          metric_quantities),
      inertial_coords);
  ::fd::general_cartoon_deriv(
      deriv_pi_evolution, get<gh::Tags::Pi<DataVector, 3>>(metric_quantities),
      inertial_coords);
  ::fd::general_cartoon_deriv(
      deriv_phi_evolution, get<gh::Tags::Phi<DataVector, 3>>(metric_quantities),
      inertial_coords);

  size_t shifted_index = 0;
  // Loop over different variables stored in u
  for (size_t component_index = 0;
       component_index < number_of_independent_components; ++component_index) {
    // loop over derivative directions
    for (size_t deriv_index = 0; deriv_index < Dim; ++deriv_index) {
      // lhs points to pdu (shifts by num grid points below)
      lhs.set_data_ref(pdu, num_grid_points);

      // clang-tidy: const cast is fine since we won't modify the data and we
      // need it to easily hook into the expression templates.

      // logical_du now points to logical_partial_derivatives_of_u
      logical_du.set_data_ref(
          const_cast<double*>(                                 // NOLINT
              gsl::at(logical_partial_derivatives_of_u, 0)) +  // NOLINT
              component_index * num_grid_points,
          num_grid_points);

      if (deriv_index == 0) {
        // usual x derivative
        lhs = inverse_jacobian.get(0, 0) * logical_du;
      } else {
        // g & Pi
        if (component_index < 20) {
          // floor division in parenthesis for positive numbers
          shifted_index = component_index - 10 * (component_index / 10);
          const auto input_tensor_index =
              get<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>>(
                  metric_quantities)
                  .get_tensor_index(shifted_index);
          auto output_tensor_index =
              prepend(input_tensor_index, size_t{deriv_index});
          if (component_index < 10) {
            // g
            lhs = deriv_spacetime_metric.get(output_tensor_index);
          } else {
            // Pi
            lhs = deriv_pi_evolution.get(output_tensor_index);
          }
        } else {
          // Phi calculation
          shifted_index = component_index - 30 * (component_index / 30);
          const auto input_tensor_index =
              get<gh::Tags::Phi<DataVector, 3>>(metric_quantities)
                  .get_tensor_index(shifted_index);
          auto output_tensor_index =
              prepend(input_tensor_index, size_t{deriv_index});
          lhs = deriv_phi_evolution.get(output_tensor_index);
        }
      }
      pdu += num_grid_points;  // NOLINT
    }
  }
}

// MIKE:
template <typename ResultTags, size_t Dim, typename DerivativeFrame>
void partial_derivatives_impl(
    const gsl::not_null<Variables<ResultTags>*> du,
    const std::array<const double*, Dim>& logical_partial_derivatives_of_u,
    const size_t number_of_independent_components,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  // pointer to where derivatives will be stored
  double* pdu = du->data();
  const size_t num_grid_points =
      du->number_of_grid_points();  // grid points per element
  DataVector lhs{};
  DataVector logical_du{};

  std::array<std::array<size_t, Dim>, Dim> indices{};
  for (size_t deriv_index = 0; deriv_index < Dim; ++deriv_index) {
    for (size_t d = 0; d < Dim; ++d) {
      gsl::at(gsl::at(indices, d), deriv_index) =
          InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>::get_storage_index(d, deriv_index);
    }
  }
  // loop over different components (variables) of u
  for (size_t component_index = 0;
       component_index < number_of_independent_components; ++component_index) {
    // loop over direction of derivatives
    for (size_t deriv_index = 0; deriv_index < Dim; ++deriv_index) {
      // pdu now points to first "num_grid_points"
      lhs.set_data_ref(pdu, num_grid_points);
      // clang-tidy: const cast is fine since we won't modify the data and we
      // need it to easily hook into the expression templates.

      // logical_du now points to logical_partial_derivatives_of_u?
      logical_du.set_data_ref(
          const_cast<double*>(                                 // NOLINT
              gsl::at(logical_partial_derivatives_of_u, 0)) +  // NOLINT
              component_index * num_grid_points,
          num_grid_points);
      // scale lhs (logical_du) by inverse jacobian--giving pdu
      lhs = (*(inverse_jacobian.begin() + gsl::at(indices[0], deriv_index))) *
            logical_du;

      for (size_t logical_deriv_index = 1; logical_deriv_index < Dim;
           ++logical_deriv_index) {
        // clang-tidy: const cast is fine since we won't modify the data
        // and we
        // need it to easily hook into the expression templates.
        logical_du.set_data_ref(const_cast<double*>(  // NOLINT
                                    gsl::at(logical_partial_derivatives_of_u,
                                            logical_deriv_index)) +
                                    component_index * num_grid_points,
                                num_grid_points);
        lhs +=
            (*(inverse_jacobian.begin() +
               gsl::at(gsl::at(indices, logical_deriv_index), deriv_index))) *
            logical_du;
      }
      // clang-tidy: no pointer arithmetic
      // shift pdu to next variable
      pdu += num_grid_points;  // NOLINT
    }
  }
}
}  // namespace partial_derivatives_detail

template <typename DerivativeTags, typename VariableTags, size_t Dim>
void logical_partial_derivatives(
    const gsl::not_null<std::array<Variables<DerivativeTags>, Dim>*>
        logical_partial_derivatives_of_u,
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
  if (UNLIKELY((*logical_partial_derivatives_of_u)[0].number_of_grid_points() !=
               u.number_of_grid_points())) {
    for (auto& deriv : *logical_partial_derivatives_of_u) {
      deriv.initialize(u.number_of_grid_points());
    }
  }
  std::array<double*, Dim> deriv_pointers{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(deriv_pointers, i) =
        gsl::at(*logical_partial_derivatives_of_u, i).data();
  }
  if constexpr (Dim == 1) {
    Variables<DerivativeTags>* temp = nullptr;
    partial_derivatives_detail::LogicalImpl<Dim, VariableTags, DerivativeTags>::
        apply(make_not_null(&deriv_pointers), temp, temp, u, mesh);
    return;
  } else {
    auto buffer = cpp20::make_unique_for_overwrite<double[]>(
        2 * u.number_of_grid_points() *
        Variables<DerivativeTags>::number_of_independent_components);
    Variables<DerivativeTags> temp0(
        &buffer[0],
        u.number_of_grid_points() *
            Variables<DerivativeTags>::number_of_independent_components);
    Variables<DerivativeTags> temp1(
        &buffer[u.number_of_grid_points() *
                Variables<DerivativeTags>::number_of_independent_components],
        u.number_of_grid_points() *
            Variables<DerivativeTags>::number_of_independent_components);
    partial_derivatives_detail::LogicalImpl<Dim, VariableTags, DerivativeTags>::
        apply(make_not_null(&deriv_pointers), &temp0, &temp1, u, mesh);
  }
}

template <typename DerivativeTags, typename VariableTags, size_t Dim>
std::array<Variables<DerivativeTags>, Dim> logical_partial_derivatives(
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
  auto logical_partial_derivatives_of_u =
      make_array<Dim>(Variables<DerivativeTags>(u.number_of_grid_points()));
  logical_partial_derivatives<DerivativeTags>(
      make_not_null(&logical_partial_derivatives_of_u), u, mesh);
  return logical_partial_derivatives_of_u;
}

template <typename ResultTags, typename DerivativeTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    const gsl::not_null<Variables<ResultTags>*> du,
    const std::array<Variables<DerivativeTags>, Dim>&
        logical_partial_derivatives_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  auto& partial_derivatives_of_u = *du;
  // For mutating compute items we must set the size.
  if (UNLIKELY(partial_derivatives_of_u.number_of_grid_points() !=
               logical_partial_derivatives_of_u[0].number_of_grid_points())) {
    partial_derivatives_of_u.initialize(
        logical_partial_derivatives_of_u[0].number_of_grid_points());
  }

  std::array<const double*, Dim> logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivs, i) =
        gsl::at(logical_partial_derivatives_of_u, i).data();
  }
  partial_derivatives_detail::partial_derivatives_impl(
      make_not_null(&partial_derivatives_of_u), logical_derivs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian);
}

template <typename ResultTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    const gsl::not_null<Variables<ResultTags>*> du,
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  using DerivativeTags =
      tmpl::front<tmpl::split_at<VariableTags, tmpl::size<ResultTags>>>;
  static_assert(
      std::is_same_v<
          tmpl::transform<ResultTags, tmpl::bind<tmpl::type_from, tmpl::_1>>,
          tmpl::transform<db::wrap_tags_in<Tags::deriv, DerivativeTags,
                                           tmpl::size_t<Dim>, DerivativeFrame>,
                          tmpl::bind<tmpl::type_from, tmpl::_1>>>);
  auto& partial_derivatives_of_u = *du;
  // For mutating compute items we must set the size.
  if (UNLIKELY(partial_derivatives_of_u.number_of_grid_points() !=
               mesh.number_of_grid_points())) {
    partial_derivatives_of_u.initialize(mesh.number_of_grid_points());
  }

  const size_t vars_size =
      u.number_of_grid_points() *
      Variables<DerivativeTags>::number_of_independent_components;
  const auto logical_derivs_data = cpp20::make_unique_for_overwrite<double[]>(
      (Dim > 1 ? (Dim + 1) : Dim) * vars_size);
  std::array<double*, Dim> logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivs, i) = &(logical_derivs_data[i * vars_size]);
  }
  Variables<DerivativeTags> temp{};
  if constexpr (Dim > 1) {
    temp.set_data_ref(&logical_derivs_data[Dim * vars_size], vars_size);
  }
  partial_derivatives_detail::LogicalImpl<
      Dim, VariableTags, DerivativeTags>::apply(make_not_null(&logical_derivs),
                                                &partial_derivatives_of_u,
                                                &temp, u, mesh);

  std::array<const double*, Dim> const_logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(const_logical_derivs, i) = gsl::at(logical_derivs, i);
  }
  partial_derivatives_detail::partial_derivatives_impl(
      make_not_null(&partial_derivatives_of_u), const_logical_derivs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian);
}

template <typename DerivativeTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
Variables<db::wrap_tags_in<Tags::deriv, DerivativeTags, tmpl::size_t<Dim>,
                           DerivativeFrame>>
partial_derivatives(
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  Variables<db::wrap_tags_in<Tags::deriv, DerivativeTags, tmpl::size_t<Dim>,
                             DerivativeFrame>>
      partial_derivatives_of_u(u.number_of_grid_points());
  partial_derivatives(make_not_null(&partial_derivatives_of_u), u, mesh,
                      inverse_jacobian);
  return partial_derivatives_of_u;
}

namespace partial_derivatives_detail {
template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<1, VariableTags, DerivativeTags> {
  static constexpr const size_t Dim = 1;
  template <typename T>
  static void apply(const gsl::not_null<std::array<double*, Dim>*> logical_du,
                    Variables<T>* /*unused_in_1d*/,
                    Variables<DerivativeTags>* const /*unused_in_1d*/,
                    const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
    auto& logical_partial_derivatives_of_u = *logical_du;
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    dgemm_<true>('N', 'N', mesh.extents(0), deriv_size / mesh.extents(0),
                 mesh.extents(0), 1.0, differentiation_matrix_xi.data(),
                 differentiation_matrix_xi.spacing(), u.data(), mesh.extents(0),
                 0.0, logical_partial_derivatives_of_u[0], mesh.extents(0));
  }
};

template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<2, VariableTags, DerivativeTags> {
  static constexpr size_t Dim = 2;
  template <typename T>
  static void apply(const gsl::not_null<std::array<double*, Dim>*> logical_du,
                    Variables<T>* const partial_u_wrt_eta,
                    Variables<DerivativeTags>* const u_eta_fastest,
                    const Variables<VariableTags>& u, const Mesh<2>& mesh) {
    static_assert(
        Variables<DerivativeTags>::number_of_independent_components <=
            Variables<T>::number_of_independent_components,
        "Temporary buffer in logical partial derivatives is too small");
    auto& logical_partial_derivatives_of_u = *logical_du;
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t num_components_times_xi_slices = deriv_size / mesh.extents(0);
    dgemm_<true>('N', 'N', mesh.extents(0), num_components_times_xi_slices,
                 mesh.extents(0), 1.0, differentiation_matrix_xi.data(),
                 differentiation_matrix_xi.spacing(), u.data(), mesh.extents(0),
                 0.0, logical_partial_derivatives_of_u[0], mesh.extents(0));

    transpose<Variables<VariableTags>, Variables<DerivativeTags>>(
        make_not_null(u_eta_fastest), u, mesh.extents(0),
        num_components_times_xi_slices);
    const Matrix& differentiation_matrix_eta =
        Spectral::differentiation_matrix(mesh.slice_through(1));
    const size_t num_components_times_eta_slices = deriv_size / mesh.extents(1);
    dgemm_<true>('N', 'N', mesh.extents(1), num_components_times_eta_slices,
                 mesh.extents(1), 1.0, differentiation_matrix_eta.data(),
                 differentiation_matrix_eta.spacing(), u_eta_fastest->data(),
                 mesh.extents(1), 0.0, partial_u_wrt_eta->data(),
                 mesh.extents(1));
    raw_transpose(make_not_null(logical_partial_derivatives_of_u[1]),
                  partial_u_wrt_eta->data(), num_components_times_xi_slices,
                  mesh.extents(0));
  }
};

template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<3, VariableTags, DerivativeTags> {
  static constexpr size_t Dim = 3;
  template <class T>
  static void apply(const gsl::not_null<std::array<double*, Dim>*> logical_du,
                    Variables<T>* const partial_u_wrt_eta_or_zeta,
                    Variables<DerivativeTags>* const u_eta_or_zeta_fastest,
                    const Variables<VariableTags>& u, const Mesh<3>& mesh) {
    static_assert(
        Variables<DerivativeTags>::number_of_independent_components <=
            Variables<T>::number_of_independent_components,
        "Temporary buffer in logical partial derivatives is too small");
    auto& logical_partial_derivatives_of_u = *logical_du;
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    const size_t num_components_times_xi_slices = deriv_size / mesh.extents(0);
    dgemm_<true>('N', 'N', mesh.extents(0), num_components_times_xi_slices,
                 mesh.extents(0), 1.0, differentiation_matrix_xi.data(),
                 differentiation_matrix_xi.spacing(), u.data(), mesh.extents(0),
                 0.0, logical_partial_derivatives_of_u[0], mesh.extents(0));

    transpose<Variables<VariableTags>, Variables<DerivativeTags>>(
        make_not_null(u_eta_or_zeta_fastest), u, mesh.extents(0),
        num_components_times_xi_slices);
    const Matrix& differentiation_matrix_eta =
        Spectral::differentiation_matrix(mesh.slice_through(1));
    const size_t num_components_times_eta_slices = deriv_size / mesh.extents(1);
    dgemm_<true>('N', 'N', mesh.extents(1), num_components_times_eta_slices,
                 mesh.extents(1), 1.0, differentiation_matrix_eta.data(),
                 differentiation_matrix_eta.spacing(),
                 u_eta_or_zeta_fastest->data(), mesh.extents(1), 0.0,
                 partial_u_wrt_eta_or_zeta->data(), mesh.extents(1));
    raw_transpose(make_not_null(logical_partial_derivatives_of_u[1]),
                  partial_u_wrt_eta_or_zeta->data(),
                  num_components_times_xi_slices, mesh.extents(0));

    const size_t chunk_size = mesh.extents(0) * mesh.extents(1);
    const size_t number_of_chunks = deriv_size / chunk_size;
    transpose(make_not_null(u_eta_or_zeta_fastest), u, chunk_size,
              number_of_chunks);
    const Matrix& differentiation_matrix_zeta =
        Spectral::differentiation_matrix(mesh.slice_through(2));
    const size_t num_components_times_zeta_slices =
        deriv_size / mesh.extents(2);
    dgemm_<true>('N', 'N', mesh.extents(2), num_components_times_zeta_slices,
                 mesh.extents(2), 1.0, differentiation_matrix_zeta.data(),
                 differentiation_matrix_zeta.spacing(),
                 u_eta_or_zeta_fastest->data(), mesh.extents(2), 0.0,
                 partial_u_wrt_eta_or_zeta->data(), mesh.extents(2));
    raw_transpose(make_not_null(logical_partial_derivatives_of_u[2]),
                  partial_u_wrt_eta_or_zeta->data(), number_of_chunks,
                  chunk_size);
  }
};
}  // namespace partial_derivatives_detail

// Macro to explicitly instantiate partial_derivatives()
// for a given system of equations
#define INSTANTIATE_PARTIAL_DERIVATIVES_WITH_SYSTEM(SYSTEM, DIM,             \
                                                    DERIVATIVE_FRAME)        \
  template void partial_derivatives(                                         \
      gsl::not_null<Variables<                                               \
          db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables, \
                           tmpl::size_t<DIM>, DERIVATIVE_FRAME>>*>           \
          du,                                                                \
      const std::array<Variables<typename SYSTEM::gradient_variables>, DIM>& \
          logical_partial_derivatives_of_u,                                  \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);            \
  template void partial_derivatives(                                         \
      gsl::not_null<Variables<                                               \
          db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables, \
                           tmpl::size_t<DIM>, DERIVATIVE_FRAME>>*>           \
          du,                                                                \
      const Variables<typename SYSTEM::gradient_variables>& u,               \
      const Mesh<DIM>& mesh,                                                 \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);            \
  template Variables<                                                        \
      db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables,     \
                       tmpl::size_t<DIM>, DERIVATIVE_FRAME>>                 \
  partial_derivatives<typename SYSTEM::gradient_variables>(                  \
      const Variables<typename SYSTEM::gradient_variables>& u,               \
      const Mesh<DIM>& mesh,                                                 \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);
