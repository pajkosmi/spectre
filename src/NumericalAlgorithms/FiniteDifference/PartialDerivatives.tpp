// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.hpp"

#include <array>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

namespace fd {
namespace detail {
template <size_t Dim>
void logical_partial_derivatives_impl(
    const gsl::not_null<std::array<gsl::span<double>, Dim>*>
        logical_derivatives,
    gsl::span<double>* buffer, const gsl::span<const double>& volume_vars,
    const DirectionMap<Dim, gsl::span<const double>>& ghost_cell_vars,
    const Mesh<Dim>& volume_mesh, size_t number_of_variables, size_t fd_order);
}  // namespace detail

template <typename DerivativeTags, size_t Dim, typename DerivativeFrame>
void partial_derivatives(
    const gsl::not_null<Variables<db::wrap_tags_in<
        Tags::deriv, DerivativeTags, tmpl::size_t<Dim>, DerivativeFrame>>*>
        partial_derivatives,
    const gsl::span<const double>& volume_vars,
    const DirectionMap<Dim, gsl::span<const double>>& ghost_cell_vars,
    const Mesh<Dim>& volume_mesh, const size_t number_of_variables,
    const size_t fd_order,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  ASSERT(partial_derivatives->size() == Dim * volume_vars.size(),
         "The partial derivatives Variables must have size "
             << Dim * volume_vars.size()
             << " (Dim * volume_vars.size()) but has size "
             << partial_derivatives->size() << " and "
             << partial_derivatives->number_of_grid_points()
             << " grid points.");
  const size_t logical_derivs_internal_buffer_size =
      Dim == 1
          ? static_cast<size_t>(0)
          : (volume_vars.size() +
             2 * alg::max_element(ghost_cell_vars,
                                  [](const auto& a, const auto& b) {
                                    return a.second.size() < b.second.size();
                                  })
                     ->second.size() +
             volume_vars.size());
  DataVector buffer(Dim * volume_vars.size() +
                    logical_derivs_internal_buffer_size);
  std::array<gsl::span<double>, Dim> logical_partial_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_partial_derivs, i) =
        gsl::make_span(&buffer[i * volume_vars.size()], volume_vars.size());
  }
  if constexpr (Dim > 1) {
    gsl::span<double> span_buffer = gsl::make_span(
        &buffer[Dim * volume_vars.size()], logical_derivs_internal_buffer_size);
    detail::logical_partial_derivatives_impl(
        make_not_null(&logical_partial_derivs), &span_buffer, volume_vars,
        ghost_cell_vars, volume_mesh, number_of_variables, fd_order);
  } else {
    // No buffer in 1d
    logical_partial_derivatives(make_not_null(&logical_partial_derivs),
                                volume_vars, ghost_cell_vars, volume_mesh,
                                number_of_variables, fd_order);
  }

  std::array<const double*, Dim> logical_partial_derivs_ptrs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_partial_derivs_ptrs, i) =
        gsl::at(logical_partial_derivs, i).data();
  }
  ::partial_derivatives_detail::partial_derivatives_impl(
      partial_derivatives, logical_partial_derivs_ptrs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian);
}

template <typename DerivativeTags, size_t Dim, typename DerivativeFrame>
void cartoon_partial_derivatives(
    const gsl::not_null<Variables<db::wrap_tags_in<
        Tags::deriv, DerivativeTags, tmpl::size_t<Dim>, DerivativeFrame>>*>
        partial_derivatives,
    const gsl::span<const double>& volume_vars,
    const DirectionMap<Dim, gsl::span<const double>>& ghost_cell_vars,
    const Mesh<Dim>& volume_mesh, const size_t number_of_variables,
    const size_t fd_order,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords) {
  const size_t logical_derivs_internal_buffer_size =
      Dim == 1
          ? static_cast<size_t>(0)
          : (volume_vars.size() +
             2 * alg::max_element(ghost_cell_vars,
                                  [](const auto& a, const auto& b) {
                                    return a.second.size() < b.second.size();
                                  })
                     ->second.size() +
             volume_vars.size());
  DataVector buffer(Dim * volume_vars.size() +
                    logical_derivs_internal_buffer_size);
  std::array<gsl::span<double>, Dim> logical_partial_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_partial_derivs, i) =
        gsl::make_span(&buffer[i * volume_vars.size()], volume_vars.size());
  }

  logical_partial_derivatives(make_not_null(&logical_partial_derivs),
                              volume_vars, ghost_cell_vars, volume_mesh,
                              number_of_variables, fd_order);

  std::array<const double*, Dim> logical_partial_derivs_ptrs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_partial_derivs_ptrs, i) =
        gsl::at(logical_partial_derivs, i).data();
  }

  ::partial_derivatives_detail::partial_derivatives_cartoon(
      partial_derivatives, logical_partial_derivs_ptrs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian, inertial_coords, volume_vars, volume_mesh.extents());
}

// The math governing the general cartoon partial derivatives follows a
// generalization of Section A.1 of \cite BaumgarteShapiro.  Equation (A.8)
// provides a sample for a mixed rank-2 tensor.
// L_xi T^a_b = xi^c \partial_c T^a_b - T^c_b \partial_c xi^a + T^a_c \partial_b
// xi^c where xi^a is a Killing vector. Asserting L_xi T^a_b = 0 to impose a
// cartoon symmetry, and setting xi^a to a given killing vector, either in the y
// or z direction, one simplifies to e.g., \partial_y g^a_b = - 1 / x *(-g^c_b *
// \partial_c \xi^a + g^a_c * \partial_b \xi^c).  Here we notice a pattern
// emerge that can be applied to higher rank, mixed tensors.  For example,
// \partial_y H^a_b^c_d = -1 / x * (-H^dummy_b^c_d * \partial_dummy \xi^a +
// H^a_dummy^c_d * \partial_b \xi^dummy - H^a_b^dummy_d * \partial_dummy \xi^c
// + H^a_b^c_dummy * \partial_d \xi^dummy).  Note the pattern: when the dummy
// index of H is contravariant, there is a negative sign and summation over the
// direction of the partial derivative of the killing vector.  When the dummy
// index of H is covariant, there is a positive sign and summation over the
// components of the Killing vector.  The number of terms in the summand is
// equal to the number of indices of H.

// cartoon general call
template <typename InputTensorDataType, typename DataType, size_t SpatialDim,
          typename Frame>
void general_cartoon_deriv(
    TensorMetafunctions::prepend_spatial_index<InputTensorDataType, SpatialDim,
                                               UpLo::Lo, Frame>& deriv_tensor,
    const InputTensorDataType& tensor,
    const tnsr::I<DataType, SpatialDim, Frame>& inertial_coords) {
  // check output rank is 1 larger than input rank
  ASSERT(deriv_tensor.rank() - tensor.rank() == 1,
         "Deriv tensor rank must be exactly 1 higher than input tensor.");

  // read valence
  const std::array<UpLo, InputTensorDataType::rank()> valences =
      tensor.index_valences();
  auto type_index = tensor.index_types();

  // Sign of addition factor
  double sign = 1.0;

  // tensor
  const size_t valence_size = valences.size();
  auto input_tensor_array = make_array<valence_size>(0_st);

  std::array<size_t, 3> Killing_indices;

  size_t max_dummy = 3;
  size_t shift_index = 0;

  // construct array containing derivatives of Killing vectors
  // \partial_deriv_index xi_i^Killing_index for Killing vectors along direction
  tnsr::iab<DataType, SpatialDim> da_killing_vectors{};

  // i = 0 is unused
  // i = 1 & 2 corresponds to spherical symmetry b/c we have 2 Killing vectors.
  // i = 1 corresponds to Killing vector xi_1 = (0, -y, x, 0) and i = 2
  // corresponds to Killing vector .
  // Only i = 2 corresponds to axisymmetry, with the only Killing vector xi_2 =
  // (0, -z, 0, x).
  for (size_t i = 0; i < 3; i++) {
    // derivatives of Killing vectors - t, x, y, z
    for (size_t deriv_index = 0; deriv_index < 4; deriv_index++) {
      // four components of Killing vectors
      for (size_t Killing_index = 0; Killing_index < 4; Killing_index++) {
        // \partial_y xi_1^x = -1 for xi along the y-direction
        if (i == 1 and deriv_index == 2 and Killing_index == 1) {
          da_killing_vectors.get(i, deriv_index, Killing_index) =
              0.0 * tensor.get(input_tensor_array) - 1.0;
        } else if (i == 1 and deriv_index == 1 and Killing_index == 2) {
          // \partial_x xi_1^y = 1
          da_killing_vectors.get(i, deriv_index, Killing_index) =
              0.0 * tensor.get(input_tensor_array) + 1.0;
        } else if (i == 2 and deriv_index == 3 and Killing_index == 1) {
          // \partial_z xi_2^x = -1
          da_killing_vectors.get(i, deriv_index, Killing_index) =
              0.0 * tensor.get(input_tensor_array) - 1.0;
        } else if (i == 2 and deriv_index == 1 and Killing_index == 3) {
          // \partial_x xi_2^z = 1
          da_killing_vectors.get(i, deriv_index, Killing_index) =
              0.0 * tensor.get(input_tensor_array) + 1.0;
        } else {
          da_killing_vectors.get(i, deriv_index, Killing_index) =
              0.0 * tensor.get(input_tensor_array);
        }
      }
    }
  }
  std::string symmetry = "Spherical";  // Axisymmetric, Spherical

  size_t max_deriv_index = 3;
  size_t start_deriv_index = symmetry == "Spherical" ? 1 : 2;

  // loop over different Killing vector(s) that depend on symmetry
  for (size_t deriv_index = start_deriv_index; deriv_index < max_deriv_index;
       deriv_index++) {
    // loop over tensor indices e.g., a, b, c, d: H^a_b^c_d
    for (size_t i = 0; i < tensor.size(); i++) {
      const auto input_tensor_index = tensor.get_tensor_index(i);
      auto output_tensor_index =
          prepend(input_tensor_index, size_t{deriv_index});
      // initialize derivative tensor cartoon terms to 0
      deriv_tensor.get(output_tensor_index) = 0.0 * inertial_coords.get(0);
      // loop over different ranks
      for (size_t rank = 0; rank < valence_size; rank++) {
        // select sign based on contravariant or covariant index
        // e.g., (+/-)H^a_b^c_d * \partial ...
        sign = (valences[rank] == UpLo::Lo) ? 1.0 : -1.0;
        // If it's a spacetime index, we have 4 dummy indices, otherwise spatial
        // has 3
        max_dummy = type_index[rank] == IndexType::Spacetime ? 4 : 3;
        // helper index for tensor algebra
        shift_index = type_index[rank] == IndexType::Spacetime ? 0 : 1;

        // loop over dimension of a given dummy rank, summing each step
        for (size_t dummy = 0; dummy < max_dummy; dummy++) {
          for (size_t j = 0; j < input_tensor_index.size(); j++) {
            input_tensor_array[j] = input_tensor_index[j];
          }
          Killing_indices[0] = deriv_index;
          if (valences[rank] == UpLo::Lo) {
            // covariant index
            Killing_indices[1] = input_tensor_index[rank] + shift_index;
            Killing_indices[2] = dummy + shift_index;
          } else {
            // contravariant
            Killing_indices[1] = dummy + shift_index;
            Killing_indices[2] = input_tensor_index[rank] + shift_index;
          }

          input_tensor_array[rank] = dummy;

          // MIKE: L'Hospital rule check here?
          // Note negative sign needed
          // \partial_... H^a_b^c_d = -1.0 * ((+/-1) * H^a_b^c_... *
          // \partial_... \xi^a + (+/-1)*...) / x.
          deriv_tensor.get(output_tensor_index) -=
              sign * tensor.get(input_tensor_array) *
              da_killing_vectors.get(Killing_indices) / inertial_coords.get(0);
        }
      }
    }
  }
};

}  // namespace fd

// TimeDerivative
