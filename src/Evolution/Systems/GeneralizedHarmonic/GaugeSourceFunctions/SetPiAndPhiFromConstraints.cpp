// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/SetPiAndPhiFromConstraints.hpp"

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/EagerMath/Trace.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/DgSubcell/Mesh.hpp"
#include "Evolution/DgSubcell/Projection.hpp"
#include "Evolution/DgSubcell/Reconstruction.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Dispatch.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Gauges.hpp"
#include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/DerivSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/GaugeSource.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Pi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpacetimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfLapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpatialDerivOfShift.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivOfSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/TimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/GenerateInstantiations.hpp"

#include <iostream>

namespace gh::gauges {
template <size_t Dim>
void SetPiAndPhiFromConstraints<Dim>::apply(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> phi,
    const double time, const Mesh<Dim>& mesh,
    const ElementMap<Dim, Frame::Grid>& logical_to_grid_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        grid_to_inertial_map,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    const tnsr::I<DataVector, Dim, Frame::ElementLogical>& logical_coordinates,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const gauges::GaugeCondition& gauge_condition) {
  const auto grid_coords = logical_to_grid_map(logical_coordinates);
  const auto inv_jac_logical_to_grid =
      logical_to_grid_map.inv_jacobian(logical_coordinates);
  const auto [inertial_coords, inv_jac_grid_to_inertial, jac_grid_to_inertial,
              mesh_velocity] =
      grid_to_inertial_map.coords_frame_velocity_jacobians(grid_coords, time,
                                                           functions_of_time);

  InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian{};

  for (size_t logical_i = 0; logical_i < Dim; ++logical_i) {
    for (size_t inertial_i = 0; inertial_i < Dim; ++inertial_i) {
      inverse_jacobian.get(logical_i, inertial_i) =
          inv_jac_logical_to_grid.get(logical_i, 0) *
          inv_jac_grid_to_inertial.get(0, inertial_i);
      for (size_t grid_i = 1; grid_i < Dim; ++grid_i) {
        inverse_jacobian.get(logical_i, inertial_i) +=
            inv_jac_logical_to_grid.get(logical_i, grid_i) *
            inv_jac_grid_to_inertial.get(grid_i, inertial_i);
      }
    }
  }

  // TODO:
  //  i) reconstruct to DG grid
  //  ii) apply partial derivative on dg grid
  //  iii) project derivative to FD grid

  // const auto dg_spacetime_metric = spacetime_metric;
  // need dg_mesh and subcell_mesh.  Can I get from mesh or do I need to pass in
  // data box?

  // Pseudo code below
  // reconstruct spacetime metric and inverse jacobian
  // src/Evolution/DgSubcell/Mesh.hpp
  const auto dg_mesh = evolution::dg::subcell::fd::dg_mesh(
      mesh, Spectral::Basis::Legendre, Spectral::Quadrature::Gauss);

  using tag_list =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>,
                 domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                               Frame::Inertial>>;

  Variables<tag_list> metric_and_inv_jac{mesh.number_of_grid_points()};

  get<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>>(
      metric_and_inv_jac) = spacetime_metric;
  get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                    Frame::Inertial>>(metric_and_inv_jac) =
      inverse_jacobian;

  // metric and inv jacobian now on DG grid
  auto dg_metric_and_inv_jac = evolution::dg::subcell::fd::reconstruct(
      metric_and_inv_jac, dg_mesh, mesh.extents(),
      evolution::dg::subcell::fd::ReconstructionMethod::DimByDim);

  // Phi on dg grid now
  auto dg_phi = partial_derivative(
      get<gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>>(
          dg_metric_and_inv_jac),
      dg_mesh,
      get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                        Frame::Inertial>>(
          dg_metric_and_inv_jac));

  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < a; ++b) {
      phi->get(0, a, b) = evolution::dg::subcell::fd::project(
          dg_phi.get(0, a, b), dg_mesh, mesh.extents());

      // zero out lower components
      phi->get(1, a, b) = 0 * phi->get(0, a, b);
      phi->get(2, a, b) = 0 * phi->get(0, a, b);
    }
  }
  // evolution::dg::subcell::fd::project(phi, dg_phi, dg_mesh, mesh.extents());

  //
  // end pseudo code

  // original phi derivative calculation
  // calculate phi
  // partial_derivative(phi, spacetime_metric, mesh, inverse_jacobian);

  // MIKE: general cartoon
  // y derivative
  //  Phi_y12 = 1 / x * (g_11 - g_22)
  // phi->get(1, 1, 2) = 1.0 / inertial_coords.get(0) *
  //                     (spacetime_metric.get(1, 1) - spacetime_metric.get(2,
  //                     2));
  // phi->get(1, 2, 1) = phi->get(1, 1, 2);

  // // z derivative
  // // Phi_z13 = 1 / x * (g_11 - g_33)
  // phi->get(2, 1, 3) = 1.0 / inertial_coords.get(0) *
  //                     (spacetime_metric.get(1, 1) - spacetime_metric.get(3,
  //                     3));
  // phi->get(2, 3, 1) = 99*phi->get(0, 2, 2);

  ::fd::general_cartoon_deriv(*phi, spacetime_metric, inertial_coords);

  Variables<
      tmpl::list<gr::Tags::SpatialMetric<DataVector, Dim>,
                 gr::Tags::SqrtDetSpatialMetric<DataVector>,
                 gr::Tags::InverseSpatialMetric<DataVector, Dim>,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
                 gr::Tags::InverseSpacetimeMetric<DataVector, Dim>,
                 gr::Tags::SpacetimeNormalVector<DataVector, Dim>,
                 ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<Dim>,
                               Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::Shift<DataVector, Dim>,
                               tmpl::size_t<Dim>, Frame::Inertial>,
                 ::Tags::deriv<gr::Tags::SpatialMetric<DataVector, Dim>,
                               tmpl::size_t<Dim>, Frame::Inertial>,
                 gr::Tags::ExtrinsicCurvature<DataVector, Dim>,
                 gr::Tags::TraceExtrinsicCurvature<DataVector>,
                 gr::Tags::SpatialChristoffelFirstKind<DataVector, Dim>,
                 gr::Tags::TraceSpatialChristoffelFirstKind<DataVector, Dim>,
                 gh::Tags::GaugeH<DataVector, Dim>,
                 gh::Tags::SpacetimeDerivGaugeH<DataVector, Dim>,
                 ::Tags::dt<gr::Tags::Lapse<DataVector>>,
                 ::Tags::dt<gr::Tags::Shift<DataVector, Dim>>,
                 ::Tags::dt<gr::Tags::SpatialMetric<DataVector, Dim>>>>
      buffer(get<0, 0>(spacetime_metric).size());

  auto& [spatial_metric, sqrt_det_spatial_metric, inverse_spatial_metric, lapse,
         shift, inverse_spacetime_metric, spacetime_unit_normal_vector, d_lapse,
         d_shift, d_spatial_metric, ex_curvature, trace_ex_curvature,
         spatial_christoffel_first, trace_spatial_christoffel_first, gauge_h,
         d4_gauge_h, dt_lapse, dt_shift, dt_spatial_metric] = buffer;

  // spatial metric from spacetime metric
  gr::spatial_metric(make_not_null(&spatial_metric), spacetime_metric);
  determinant_and_inverse(make_not_null(&sqrt_det_spatial_metric),
                          make_not_null(&inverse_spatial_metric),
                          spatial_metric);
  // calc sqrt det g
  get(sqrt_det_spatial_metric) = sqrt(get(sqrt_det_spatial_metric));
  // shift from spacetime metric [x]
  gr::shift(make_not_null(&shift), spacetime_metric, inverse_spatial_metric);
  //" " from lapse
  gr::lapse(make_not_null(&lapse), shift, spacetime_metric);
  // inverse from lapse and shift [x]
  gr::inverse_spacetime_metric(make_not_null(&inverse_spacetime_metric), lapse,
                               shift, inverse_spatial_metric);
  // normal from lapse and shift [x]
  gr::spacetime_normal_vector(make_not_null(&spacetime_unit_normal_vector),
                              lapse, shift);
  // spatial derivative of lapse from lapse, normal vector, and phi [x]
  spatial_deriv_of_lapse(make_not_null(&d_lapse), lapse,
                         spacetime_unit_normal_vector, *phi);
  // shift derivative from inverse  [x]
  spatial_deriv_of_shift(make_not_null(&d_shift), lapse,
                         inverse_spacetime_metric, spacetime_unit_normal_vector,
                         *phi);
  deriv_spatial_metric(make_not_null(&d_spatial_metric), *phi);

  extrinsic_curvature(make_not_null(&ex_curvature),
                      spacetime_unit_normal_vector, *pi, *phi);
  trace(make_not_null(&trace_ex_curvature), ex_curvature,
        inverse_spatial_metric);
  gr::christoffel_first_kind(make_not_null(&spatial_christoffel_first),
                             d_spatial_metric);
  // Gamma_a = Gamma_ab g^bc
  trace_last_indices(make_not_null(&trace_spatial_christoffel_first),
                     spatial_christoffel_first, inverse_spatial_metric);

  // Here we use `derivatives_of_spacetime_metric` to get \f$ \partial_a
  // g_{bc}\f$ instead, and use only the derivatives of \f$ g_{bi}\f$. [x]
  tnsr::abb<DataVector, Dim, Frame::Inertial> d4_spacetime_metric{};
  gh::spacetime_derivative_of_spacetime_metric(
      make_not_null(&d4_spacetime_metric), lapse, shift, *pi, *phi);

  Scalar<DataVector> half_pi_two_normals{get(lapse).size(), 0.0};
  tnsr::i<DataVector, Dim, Frame::Inertial> half_phi_two_normals{
      get(lapse).size(), 0.0};
  for (size_t a = 0; a < Dim + 1; ++a) {
    get(half_pi_two_normals) += spacetime_unit_normal_vector.get(a) *
                                spacetime_unit_normal_vector.get(a) *
                                pi->get(a, a);
    for (size_t i = 0; i < Dim; ++i) {
      half_phi_two_normals.get(i) += 0.5 * spacetime_unit_normal_vector.get(a) *
                                     spacetime_unit_normal_vector.get(a) *
                                     phi->get(i, a, a);
    }
    for (size_t b = a + 1; b < Dim + 1; ++b) {
      get(half_pi_two_normals) += 2.0 * spacetime_unit_normal_vector.get(a) *
                                  spacetime_unit_normal_vector.get(b) *
                                  pi->get(a, b);
      for (size_t i = 0; i < Dim; ++i) {
        half_phi_two_normals.get(i) += spacetime_unit_normal_vector.get(a) *
                                       spacetime_unit_normal_vector.get(b) *
                                       phi->get(i, a, b);
      }
    }
  }
  get(half_pi_two_normals) *= 0.5;

  // Note: we pass in pi to compute d4_gauge_h, but we don't use d4_gauge_h. We
  // actually reset pi from gauge_h below.
  dispatch(make_not_null(&gauge_h), make_not_null(&d4_gauge_h), lapse, shift,
           sqrt_det_spatial_metric, inverse_spatial_metric, d4_spacetime_metric,
           half_pi_two_normals, half_phi_two_normals, spacetime_metric, *phi,
           mesh, time, inertial_coords, inverse_jacobian, gauge_condition);

  // Compute lapse and shift time derivatives
  get(dt_lapse) =
      -get(lapse) * (get<0>(gauge_h) + get(lapse) * get(trace_ex_curvature));
  for (size_t i = 0; i < Dim; ++i) {
    get(dt_lapse) +=
        shift.get(i) * (d_lapse.get(i) + get(lapse) * gauge_h.get(i + 1));
  }

  for (size_t i = 0; i < Dim; ++i) {
    dt_shift.get(i) = get<0>(shift) * d_shift.get(0, i);
    for (size_t k = 1; k < Dim; ++k) {
      dt_shift.get(i) += shift.get(k) * d_shift.get(k, i);
    }
    for (size_t j = 0; j < Dim; ++j) {
      dt_shift.get(i) +=
          get(lapse) * inverse_spatial_metric.get(i, j) *
          (get(lapse) *
               (gauge_h.get(j + 1) + trace_spatial_christoffel_first.get(j)) -
           d_lapse.get(j));
    }
  }

  time_deriv_of_spatial_metric(make_not_null(&dt_spatial_metric), lapse, shift,
                               *phi, *pi);
  // Mike: want to match gauge_h with trace_spatial_christoffel_first.  Should
  // they be the same?
  gh::pi(pi, lapse, dt_lapse, shift, dt_shift, spatial_metric,
         dt_spatial_metric, *phi);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) \
  template class SetPiAndPhiFromConstraints<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
}  // namespace gh::gauges
