// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <random>
#include <unordered_set>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Evolution/DgSubcell/SliceData.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.tpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"

namespace {
template <size_t Dim>
void test(const gsl::not_null<std::mt19937*> generator,
          const gsl::not_null<std::uniform_real_distribution<>*> dist,
          const size_t points_per_dimension, const size_t fd_order) {
  CAPTURE(points_per_dimension);
  CAPTURE(fd_order);
  CAPTURE(Dim);
  const size_t max_degree = fd_order - 1;
  const size_t stencil_width = fd_order + 1;
  const size_t number_of_vars = 2;  // arbitrary, 2 is "cheap but not trivial"

  const Mesh<Dim> mesh{points_per_dimension, Spectral::Basis::FiniteDifference,
                       Spectral::Quadrature::CellCentered};
  auto logical_coords = logical_coordinates(mesh);
  // Make the logical coordinates different in each direction
  for (size_t i = 1; i < Dim; ++i) {
    logical_coords.get(i) += 4.0 * i;
  }

  // Compute polynomial on cell centers in FD cluster of points
  const auto set_polynomial = [max_degree](
                                  const gsl::not_null<DataVector*> var1_ptr,
                                  const gsl::not_null<DataVector*> var2_ptr,
                                  const auto& local_logical_coords) {
    *var1_ptr = 0.0;
    *var2_ptr = 100.0;  // some constant offset to distinguish the var values
    for (size_t degree = 1; degree <= max_degree; ++degree) {
      for (size_t i = 0; i < Dim; ++i) {
        *var1_ptr += pow(local_logical_coords.get(i), degree);
        *var2_ptr += pow(local_logical_coords.get(i), degree);
      }
    }
  };
  const auto set_polynomial_derivative =
      [max_degree](const gsl::not_null<std::array<DataVector, Dim>*> d_var1_ptr,
                   const gsl::not_null<std::array<DataVector, Dim>*> d_var2_ptr,
                   const auto& local_logical_coords) {
        for (size_t deriv_dim = 0; deriv_dim < Dim; ++deriv_dim) {
          gsl::at(*d_var1_ptr, deriv_dim) = 0.0;
          // constant deriv is zero
          gsl::at(*d_var2_ptr, deriv_dim) = 0.0;
          for (size_t degree = 1; degree <= max_degree; ++degree) {
            gsl::at(*d_var1_ptr, deriv_dim) +=
                degree * pow(local_logical_coords.get(deriv_dim), degree - 1);
            gsl::at(*d_var2_ptr, deriv_dim) +=
                degree * pow(local_logical_coords.get(deriv_dim), degree - 1);
          }
        }
      };

  DataVector volume_vars{mesh.number_of_grid_points() * number_of_vars, 0.0};
  DataVector var1(volume_vars.data(), mesh.number_of_grid_points());
  DataVector var2(volume_vars.data() + mesh.number_of_grid_points(),  // NOLINT
                  mesh.number_of_grid_points());
  set_polynomial(&var1, &var2, logical_coords);

  DataVector expected_deriv{Dim * volume_vars.size()};
  std::array<DataVector, Dim> expected_d_var1{};
  std::array<DataVector, Dim> expected_d_var2{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(expected_d_var1, i)
        .set_data_ref(&expected_deriv[i * volume_vars.size()],
                      mesh.number_of_grid_points());
    gsl::at(expected_d_var2, i)
        .set_data_ref(&expected_deriv[i * volume_vars.size() +
                                      mesh.number_of_grid_points()],
                      mesh.number_of_grid_points());
  }
  set_polynomial_derivative(&expected_d_var1, &expected_d_var2, logical_coords);

  // Compute the polynomial at the cell center for the neighbor data that we
  // "received".
  //
  // We do this by computing the solution in our entire neighbor, then using
  // slice_data to get the subset of points that are needed.
  DirectionMap<Dim, DataVector> neighbor_data{};
  for (const auto& direction : Direction<Dim>::all_directions()) {
    auto neighbor_logical_coords = logical_coords;
    neighbor_logical_coords.get(direction.dimension()) +=
        direction.sign() * 2.0;
    DataVector neighbor_vars{mesh.number_of_grid_points() * number_of_vars,
                             0.0};
    DataVector neighbor_var1(neighbor_vars.data(),
                             mesh.number_of_grid_points());
    DataVector neighbor_var2(
        neighbor_vars.data() + mesh.number_of_grid_points(),  // NOLINT
        mesh.number_of_grid_points());
    set_polynomial(&neighbor_var1, &neighbor_var2, neighbor_logical_coords);

    const auto sliced_data = evolution::dg::subcell::detail::slice_data_impl(
        gsl::make_span(neighbor_vars.data(), neighbor_vars.size()),
        mesh.extents(), (stencil_width - 1) / 2 + 1,
        std::unordered_set{direction.opposite()}, 0, {});
    CAPTURE((stencil_width - 1) / 2 + 1);
    REQUIRE(sliced_data.size() == 1);
    REQUIRE(sliced_data.contains(direction.opposite()));
    neighbor_data[direction] = sliced_data.at(direction.opposite());
    REQUIRE(neighbor_data.at(direction).size() ==
            number_of_vars * (fd_order / 2 + 1) *
                mesh.slice_away(0).number_of_grid_points());
  }

  // Note: reconstructed_num_pts assumes isotropic extents
  DataVector logical_derivative_buffer{volume_vars.size() * Dim};
  std::array<gsl::span<double>, Dim> logical_derivative_view{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivative_view, i) = gsl::make_span(
        &logical_derivative_buffer[i * volume_vars.size()], volume_vars.size());
  }

  DirectionMap<Dim, gsl::span<const double>> ghost_cell_vars{};
  for (const auto& [direction, data] : neighbor_data) {
    ghost_cell_vars[direction] = gsl::make_span(data.data(), data.size());
  }

  ::fd::logical_partial_derivatives(
      make_not_null(&logical_derivative_view),
      gsl::make_span(volume_vars.data(), volume_vars.size()), ghost_cell_vars,
      mesh, number_of_vars, fd_order);

  // Scale to volume_vars since that sets the subtraction error threshold.
  Approx custom_approx = Approx::custom().epsilon(1.0e-14).scale(
      *std::max_element(volume_vars.begin(), volume_vars.end()));

  for (size_t i = 0; i < Dim; ++i) {
    CAPTURE(i);
    {
      CAPTURE(var1);
      const DataVector fd_d_var1(&gsl::at(logical_derivative_view, i)[0],
                                 mesh.number_of_grid_points());
      CHECK_ITERABLE_CUSTOM_APPROX(fd_d_var1, gsl::at(expected_d_var1, i),
                                   custom_approx);
    }
    {
      CAPTURE(var2);
      const DataVector fd_d_var2(
          &gsl::at(logical_derivative_view, i)[mesh.number_of_grid_points()],
          mesh.number_of_grid_points());
      CHECK_ITERABLE_CUSTOM_APPROX(fd_d_var2, gsl::at(expected_d_var2, i),
                                   custom_approx);
    }
  }

  // Test partial derivative with random Jacobian. We know we calculated the
  // logical partial derivatives correctly, just need to make sure we forward to
  // the other functions correctly.
  const auto inverse_jacobian = make_with_random_values<
      InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>>(
      generator, dist, DataVector{mesh.number_of_grid_points()});

  using derivative_tags = tmpl::list<Tags::TempScalar<0, DataVector>,
                                     Tags::TempScalar<1, DataVector>>;
  Variables<db::wrap_tags_in<Tags::deriv, derivative_tags, tmpl::size_t<Dim>,
                             Frame::Inertial>>
      partial_derivatives{mesh.number_of_grid_points()};
  ::fd::partial_derivatives<derivative_tags>(
      make_not_null(&partial_derivatives),
      gsl::make_span(volume_vars.data(), volume_vars.size()), ghost_cell_vars,
      mesh, number_of_vars, fd_order, inverse_jacobian);

  std::array<const double*, Dim> expected_logical_derivs_ptrs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(expected_logical_derivs_ptrs, i) =
        gsl::at(expected_d_var1, i).data();
  }
  Variables<db::wrap_tags_in<Tags::deriv, derivative_tags, tmpl::size_t<Dim>,
                             Frame::Inertial>>
      expected_partial_derivatives{mesh.number_of_grid_points()};
  ::partial_derivatives_detail::partial_derivatives_impl(
      make_not_null(&expected_partial_derivatives),
      expected_logical_derivs_ptrs,
      Variables<derivative_tags>::number_of_independent_components,
      inverse_jacobian);

  using d_var1_tag = Tags::deriv<Tags::TempScalar<0, DataVector>,
                                 tmpl::size_t<Dim>, Frame::Inertial>;
  using d_var2_tag = Tags::deriv<Tags::TempScalar<1, DataVector>,
                                 tmpl::size_t<Dim>, Frame::Inertial>;
  CHECK_ITERABLE_CUSTOM_APPROX(get<d_var1_tag>(partial_derivatives),
                               get<d_var1_tag>(expected_partial_derivatives),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<d_var2_tag>(partial_derivatives),
                               get<d_var2_tag>(expected_partial_derivatives),
                               custom_approx);
}
void test_cartoon() {
  // random number tools
  std::random_device rand;
  std::mt19937 gen(rand());
  double lower = 0.1;
  double upper = 100;
  std::uniform_real_distribution<double> unif(lower, upper);

  // sample inertial coords
  tnsr::I<DataVector, 3> inertial_coords{4};

  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    get<0>(inertial_coords)[i] = unif(gen);
    get<1>(inertial_coords)[i] = 0.0;
    get<2>(inertial_coords)[i] = 0.0;
  }
  auto rad = inertial_coords.get(0);

  // zeros for reference
  tnsr::I<DataVector, 3, Frame::Inertial> zeros{inertial_coords.size() + 1,
                                                0.0};
  // sample tensors (rank 1)
  tnsr::I<DataVector, 3, Frame::Inertial> rank1_contra{
      inertial_coords.size() + 1, 0.0};
  tnsr::iJ<DataVector, 3, Frame::Inertial> d_rank1_contra_cartoon{
      inertial_coords.size() + 1, 11.0};
  const auto d_rank1_contra_reference = d_rank1_contra_cartoon;

  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    get<0>(rank1_contra)[i] = unif(gen);
    get<1>(rank1_contra)[i] = unif(gen);
    get<2>(rank1_contra)[i] = unif(gen);
  }

  ::fd::general_cartoon_deriv(d_rank1_contra_cartoon, rank1_contra,
                              inertial_coords);

  // rank 2 contravariant derivative test tensor (deriv_index, component)
  for (size_t i = 0; i < inertial_coords.size(); i++) {
    // x deriv should remain unchanged
    CHECK(d_rank1_contra_cartoon.get(0, i) ==
          d_rank1_contra_reference.get(0, i));
    // y and z derivs should change
    CHECK(d_rank1_contra_cartoon.get(1, i) !=
          d_rank1_contra_reference.get(1, i));
    CHECK(d_rank1_contra_cartoon.get(2, i) !=
          d_rank1_contra_reference.get(2, i));
  }

  // y derivs
  // dy V^x = -V^y / x
  CHECK(d_rank1_contra_cartoon.get(1, 0) == -rank1_contra.get(1) / rad);
  // dy V^y = V^x / x
  CHECK(d_rank1_contra_cartoon.get(1, 1) == rank1_contra.get(0) / rad);
  // dy V^z = 0.0
  CHECK(d_rank1_contra_cartoon.get(1, 2) == zeros.get(0));

  // z derivs
  // dz V^x = -V^z / x
  CHECK(d_rank1_contra_cartoon.get(2, 0) == -rank1_contra.get(2) / rad);
  // dz V^y = 0.0
  CHECK(d_rank1_contra_cartoon.get(2, 1) == zeros.get(0));
  // dz V^z = V^x / x
  CHECK(d_rank1_contra_cartoon.get(2, 2) == rank1_contra.get(0) / rad);
  // end rank 1 tests

  // rank 2 tests
  // spatial metric-like quantities
  tnsr::aa<DataVector, 3, Frame::Inertial> rank2_co{inertial_coords.size() + 1,
                                                    0.0};
  tnsr::iaa<DataVector, 3, Frame::Inertial> d_rank2_co_cartoon{
      inertial_coords.size() + 1, 22.0};
  const auto d_rank2_co_cartoon_reference = d_rank2_co_cartoon;

  // random values for rank 2 tensor
  for (size_t a = 0; a < 4; a++) {
    for (size_t b = 0; b < 4; b++) {
      for (size_t i = 0; i <= inertial_coords.size(); i++) {
        rank2_co.get(a, b)[i] = unif(gen);
      }
    }
  }

  ::fd::general_cartoon_deriv(d_rank2_co_cartoon, rank2_co, inertial_coords);

  // make sure things should be unchanged remain unchanged & should change have
  // changed
  for (size_t a = 0; a < 4; a++) {
    for (size_t b = 0; b < 4; b++) {
      // x derivatives should remain constant
      CHECK(d_rank2_co_cartoon_reference.get(0, a, b) ==
            d_rank2_co_cartoon.get(0, a, b));
      // y & z derivatives should change
      CHECK(d_rank2_co_cartoon_reference.get(1, a, b) !=
            d_rank2_co_cartoon.get(1, a, b));
      CHECK(d_rank2_co_cartoon_reference.get(2, a, b) !=
            d_rank2_co_cartoon.get(2, a, b));
    }
  }

  // y derivs
  // dy g00 = 0
  CHECK(d_rank2_co_cartoon.get(1, 0, 0) == zeros.get(0));
  // dy g01 = -g02
  CHECK(d_rank2_co_cartoon.get(1, 0, 1) == -rank2_co.get(0, 2) / rad);
  // dy g02 = g01
  CHECK(d_rank2_co_cartoon.get(1, 0, 2) == rank2_co.get(0, 1) / rad);
  // dy g03 = 0
  CHECK(d_rank2_co_cartoon.get(1, 0, 3) == zeros.get(0));

  // dy g10 = -g20 (symmetry)
  CHECK(d_rank2_co_cartoon.get(1, 1, 0) == -rank2_co.get(2, 0) / rad);
  // dy g11 = -2g12
  CHECK(d_rank2_co_cartoon.get(1, 1, 1) == -2.0 * rank2_co.get(1, 2) / rad);
  // dy g12 = g11 - g22
  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    CHECK(d_rank2_co_cartoon.get(1, 1, 2)[i] ==
          approx((rank2_co.get(1, 1)[i] - rank2_co.get(2, 2)[i]) / rad[i]));
  }
  // dy g13 = -g23
  CHECK(d_rank2_co_cartoon.get(1, 1, 3) == -rank2_co.get(2, 3) / rad);

  // dy g20 = g01 (symmetric)
  CHECK(d_rank2_co_cartoon.get(1, 2, 0) == rank2_co.get(0, 1) / rad);
  // dy g21 = g11 - g22 (symmetric)
  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    CHECK(d_rank2_co_cartoon.get(1, 2, 1)[i] ==
          approx((rank2_co.get(1, 1)[i] - rank2_co.get(2, 2)[i]) / rad[i]));
  }
  // dy g22 = 2 g12
  CHECK(d_rank2_co_cartoon.get(1, 2, 2) == 2.0 * rank2_co.get(1, 2) / rad);
  // dy g23 = g13
  CHECK(d_rank2_co_cartoon.get(1, 2, 3) == rank2_co.get(1, 3) / rad);

  // dy g30 = 0
  CHECK(d_rank2_co_cartoon.get(1, 3, 0) == zeros.get(0));
  // dy g31 = -g23
  CHECK(d_rank2_co_cartoon.get(1, 3, 1) == -rank2_co.get(2, 3) / rad);
  // dy g32 = g13
  CHECK(d_rank2_co_cartoon.get(1, 3, 2) == rank2_co.get(1, 3) / rad);
  // dy g33 = 0
  CHECK(d_rank2_co_cartoon.get(1, 3, 3) == zeros.get(0));
  // end y derivs

  // begin z derivs
  // dz g00 = 0
  CHECK(d_rank2_co_cartoon.get(2, 0, 0) == zeros.get(0));
  // dz g01 = -g03
  CHECK(d_rank2_co_cartoon.get(2, 0, 1) == -rank2_co.get(0, 3) / rad);
  // dz g02 = 0
  CHECK(d_rank2_co_cartoon.get(2, 0, 2) == zeros.get(0));
  // dz g03 = g01
  CHECK(d_rank2_co_cartoon.get(2, 0, 3) == rank2_co.get(0, 1) / rad);
  // end z derivs

  // dz g10 = -g03 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 1, 0) == -rank2_co.get(0, 3) / rad);
  // dz g11 = -2g13
  CHECK(d_rank2_co_cartoon.get(2, 1, 1) == -2.0 * rank2_co.get(1, 3) / rad);
  // dz g12 = -g23
  CHECK(d_rank2_co_cartoon.get(2, 1, 2) == -rank2_co.get(2, 3) / rad);
  // dz g13 = g11 - g33
  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    CHECK(d_rank2_co_cartoon.get(2, 1, 3)[i] ==
          approx((rank2_co.get(1, 1)[i] - rank2_co.get(3, 3)[i]) / rad[i]));
  }

  // dz g20 = 0 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 2, 0) == zeros.get(0));
  // dz g21 = -g23 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 2, 1) == -rank2_co.get(2, 3) / rad);
  // dz g22 = 0
  CHECK(d_rank2_co_cartoon.get(2, 2, 2) == zeros.get(0));
  // dz g23 = g12
  CHECK(d_rank2_co_cartoon.get(2, 2, 3) == rank2_co.get(1, 2) / rad);

  // dz g20 = 0 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 2, 0) == zeros.get(0));
  // dz g21 = -g23 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 2, 1) == -rank2_co.get(2, 3) / rad);
  // dz g22 = 0
  CHECK(d_rank2_co_cartoon.get(2, 2, 2) == zeros.get(0));
  // dz g23 = g12
  CHECK(d_rank2_co_cartoon.get(2, 2, 3) == rank2_co.get(1, 2) / rad);

  // dz g30 = g01 (symmetric)
  CHECK(d_rank2_co_cartoon.get(2, 3, 0) == rank2_co.get(0, 1) / rad);
  // dz g31 = g11 - g33 (symmetric)
  for (size_t i = 0; i <= inertial_coords.size(); i++) {
    CHECK(d_rank2_co_cartoon.get(2, 3, 1)[i] ==
          approx((rank2_co.get(1, 1)[i] - rank2_co.get(3, 3)[i]) / rad[i]));
  }
  // dz g32 = g12
  CHECK(d_rank2_co_cartoon.get(2, 3, 2) == rank2_co.get(1, 2) / rad);
  // dz g33 = 2 g13
  CHECK(d_rank2_co_cartoon.get(2, 3, 3) == 2.0 * rank2_co.get(1, 3) / rad);

  // end z derivs

  // begin rank3 tests

  // end rank3 tests

  // ::fd::general_cartoon_deriv(di_gauge_h, *gauge_h, inertial_coords);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.FiniteDifference.PartialDerivatives",
                  "[Unit][NumericalAlgorithms]") {
  MAKE_GENERATOR(generator);
  test_cartoon();

  // std::uniform_real_distribution<> dist{-1.0, 1.0};
  // for (const size_t fd_order : {2_st, 4_st, 6_st, 8_st}) {
  //   test<1>(make_not_null(&generator), make_not_null(&dist), fd_order + 2,
  //           fd_order);
  //   test<2>(make_not_null(&generator), make_not_null(&dist), fd_order + 2,
  //           fd_order);
  //   test<3>(make_not_null(&generator), make_not_null(&dist), fd_order + 2,
  //           fd_order);
  // }
}
