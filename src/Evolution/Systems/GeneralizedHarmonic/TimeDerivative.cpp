// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/EagerMath/Trace.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/DuDtTempTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Dispatch.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Gauges.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpacetimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

#include <iostream>

namespace gh {
template <size_t Dim>
void TimeDerivative<Dim>::apply(
    const gsl::not_null<tnsr::aa<DataVector, Dim>*> dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*> dt_pi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim>*> dt_phi,
    const gsl::not_null<Scalar<DataVector>*> temp_gamma1,
    const gsl::not_null<Scalar<DataVector>*> temp_gamma2,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> gauge_function,
    const gsl::not_null<tnsr::ab<DataVector, Dim>*>
        spacetime_deriv_gauge_function,
    const gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
    const gsl::not_null<Scalar<DataVector>*> half_pi_two_normals,
    const gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
    const gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> pi_one_normal,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> gauge_constraint,
    const gsl::not_null<tnsr::i<DataVector, Dim>*> half_phi_two_normals,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*>
        shift_dot_three_index_constraint,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*>
        mesh_velocity_dot_three_index_constraint,
    const gsl::not_null<tnsr::ia<DataVector, Dim>*> phi_one_normal,
    const gsl::not_null<tnsr::aB<DataVector, Dim>*> pi_2_up,
    const gsl::not_null<tnsr::iaa<DataVector, Dim>*> three_index_constraint,
    const gsl::not_null<tnsr::Iaa<DataVector, Dim>*> phi_1_up,
    const gsl::not_null<tnsr::iaB<DataVector, Dim>*> phi_3_up,
    const gsl::not_null<tnsr::abC<DataVector, Dim>*>
        christoffel_first_kind_3_up,
    const gsl::not_null<Scalar<DataVector>*> lapse,
    const gsl::not_null<tnsr::I<DataVector, Dim>*> shift,
    const gsl::not_null<tnsr::II<DataVector, Dim>*> inverse_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
    const gsl::not_null<tnsr::AA<DataVector, Dim>*> inverse_spacetime_metric,
    const gsl::not_null<tnsr::abb<DataVector, Dim>*> christoffel_first_kind,
    const gsl::not_null<tnsr::Abb<DataVector, Dim>*> christoffel_second_kind,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> trace_christoffel,
    const gsl::not_null<tnsr::A<DataVector, Dim>*> normal_spacetime_vector,
    const tnsr::iaa<DataVector, Dim>& d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim>& d_pi,
    const tnsr::ijaa<DataVector, Dim>& d_phi,
    const tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::aa<DataVector, Dim>& pi,
    const tnsr::iaa<DataVector, Dim>& phi_test,
    const Scalar<DataVector>& gamma0, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2,
    const gauges::GaugeCondition& gauge_condition, const Mesh<Dim>& mesh,
    double time,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        mesh_velocity) {
  const size_t number_of_points = get<0, 0>(*dt_spacetime_metric).size();
  // Need constraint damping on interfaces in DG schemes
  *temp_gamma1 = gamma1;
  *temp_gamma2 = gamma2;

  const tnsr::ii<DataVector, Dim> spatial_metric{};
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      make_const_view(make_not_null(&spatial_metric.get(i, j)),
                      spacetime_metric.get(i + 1, j + 1), 0, number_of_points);
    }
  }

  // hack phi
  auto phi = phi_test;

  // y (1) derivatives
  // gtt
  phi.get(1, 0, 0) = 0.0 * spacetime_metric.get(0, 0);
  // gtx
  phi.get(1, 0, 1) = -spacetime_metric.get(0, 2) / inertial_coords.get(0);
  // Symmetry
  phi.get(1, 1, 0) = phi.get(1, 0, 1);
  // gty
  phi.get(1, 0, 2) = spacetime_metric.get(0, 1) / inertial_coords.get(0);
  // Symmetry
  phi.get(1, 2, 0) = phi.get(1, 0, 2);
  // gtz
  phi.get(1, 0, 3) = 0.0 * spacetime_metric.get(0, 2);
  // Symmetry
  phi.get(1, 3, 0) = phi.get(1, 0, 3);
  // gxx
  phi.get(1, 1, 1) = -2.0 * spatial_metric.get(0, 1) / inertial_coords.get(0);
  // gxy
  phi.get(1, 1, 2) = (spatial_metric.get(0, 0) - spatial_metric.get(1, 1)) /
                     inertial_coords.get(0);
  // gyx = gxy symmetry
  phi.get(1, 2, 1) = phi.get(1, 1, 2);
  // gxz
  phi.get(1, 1, 3) = -spatial_metric.get(1, 2) / inertial_coords.get(0);
  // gzx = gxz symmetry
  phi.get(1, 3, 1) = phi.get(1, 1, 3);
  // gyy
  phi.get(1, 2, 2) = 2.0 * spatial_metric.get(0, 1) / inertial_coords.get(0);
  // gyz
  phi.get(1, 2, 3) = spatial_metric.get(0, 2) / inertial_coords.get(0);
  // gzy = gyz symmetry
  phi.get(1, 3, 2) = phi.get(1, 2, 3);
  // gzz                                  !just zero
  phi.get(1, 3, 3) = 0.0 * phi.get(1, 2, 2);

  // z (2) derivatives
  // gtt
  phi.get(2, 0, 0) = 0.0 * spacetime_metric.get(0, 0);
  // gtx
  phi.get(2, 0, 1) = -spacetime_metric.get(0, 3) / inertial_coords.get(0);
  // Symmetry
  phi.get(2, 1, 0) = phi.get(2, 0, 1);
  // gty
  phi.get(2, 0, 2) = 0.0 * spacetime_metric.get(0, 2);
  // Symmetry
  phi.get(2, 2, 0) = phi.get(2, 0, 2);
  // gtz
  phi.get(2, 0, 3) = spacetime_metric.get(0, 1) / inertial_coords.get(0);
  // Symmetry
  phi.get(2, 3, 0) = phi.get(2, 0, 3);
  // gxx
  phi.get(2, 1, 1) = -2.0 * spatial_metric.get(0, 2) / inertial_coords.get(0);
  // gxy
  phi.get(2, 1, 2) = -spatial_metric.get(1, 2) / inertial_coords.get(0);
  // gyx = gxy symmetry
  phi.get(2, 2, 1) = phi.get(2, 1, 2);
  // gxz
  phi.get(2, 1, 3) = (spatial_metric.get(0, 0) - spatial_metric.get(2, 2)) /
                     inertial_coords.get(0);
  // gzx = gxz symmetry
  phi.get(2, 3, 1) = phi.get(2, 1, 3);
  // gyy                                  !just zero
  phi.get(2, 2, 2) = phi.get(1, 1, 1) - phi.get(1, 1, 1);
  // gyz
  phi.get(2, 2, 3) = spatial_metric.get(0, 1) / inertial_coords.get(0);
  // gzy = gyz symmetry
  phi.get(2, 3, 2) = phi.get(2, 2, 3);
  // gzz
  phi.get(2, 3, 3) = 2.0 * spatial_metric.get(0, 2) / inertial_coords.get(0);

  // end hack phi

  // Mike: somehow hack pi?  How to account for time derivative
  // Mike: shift should be zero

  determinant_and_inverse(det_spatial_metric, inverse_spatial_metric,
                          spatial_metric);
  gr::shift(shift, spacetime_metric, *inverse_spatial_metric);
  gr::lapse(lapse, *shift, spacetime_metric);
  gr::inverse_spacetime_metric(inverse_spacetime_metric, *lapse, *shift,
                               *inverse_spatial_metric);
  // Compute the part of the dt_spacetime_metric equation that doesn't involve
  // constraints so we can use it for da_spacetime_metric to compute Christoffel
  // symbols.
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      // Mike: hack pi to zero
      dt_spacetime_metric->get(mu, nu) = -get(*lapse) * pi.get(mu, nu);
      for (size_t m = 0; m < Dim; ++m) {
        dt_spacetime_metric->get(mu, nu) += shift->get(m) * phi.get(m, mu, nu);
        // std::cout << " m n " << mu << " " << nu << " " << pi.get(mu, nu) <<
        // "\n";
      }
    }
  }
  // std::cout << "shift " << shift->get(0) << " " << shift->get(1) << " " <<
  // shift->get(2) << "\n";

  // assign spatial metric derivative based on phi
  // Mike: derivative should have nonzero terms for d_y g_12 and d_z g_13
  const std::optional da_spacetime_metric{tnsr::abb<DataVector, Dim>{}};
  for (size_t a = 0; a < Dim + 1; ++a) {
    for (size_t b = a; b < Dim + 1; ++b) {
      make_const_view(make_not_null(&da_spacetime_metric.value().get(0, a, b)),
                      dt_spacetime_metric->get(a, b), 0, number_of_points);
      for (size_t i = 0; i < Dim; ++i) {
        make_const_view(
            make_not_null(&da_spacetime_metric.value().get(i + 1, a, b)),
            phi.get(i, a, b), 0, number_of_points);
      }
    }
  }

  // Christoffel calculations
  gr::christoffel_first_kind(christoffel_first_kind,
                             da_spacetime_metric.value());
  trace_last_indices(trace_christoffel, *christoffel_first_kind,
                     *inverse_spacetime_metric);
  gr::spacetime_normal_vector(normal_spacetime_vector, *lapse, *shift);

  get(*gamma1gamma2) = get(gamma1) * get(gamma2);
  const DataVector& gamma12 = get(*gamma1gamma2);

  // raise Phi
  for (size_t m = 0; m < Dim; ++m) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        phi_1_up->get(m, mu, nu) =
            inverse_spatial_metric->get(m, 0) * phi.get(0, mu, nu);
        for (size_t n = 1; n < Dim; ++n) {
          phi_1_up->get(m, mu, nu) +=
              inverse_spatial_metric->get(m, n) * phi.get(n, mu, nu);
        }
      }
    }
  }

  for (size_t m = 0; m < Dim; ++m) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
        phi_3_up->get(m, nu, alpha) =
            inverse_spacetime_metric->get(alpha, 0) * phi.get(m, nu, 0);
        for (size_t beta = 1; beta < Dim + 1; ++beta) {
          phi_3_up->get(m, nu, alpha) +=
              inverse_spacetime_metric->get(alpha, beta) * phi.get(m, nu, beta);
        }
      }
    }
  }

  // raise Pi
  for (size_t nu = 0; nu < Dim + 1; ++nu) {
    for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
      pi_2_up->get(nu, alpha) =
          inverse_spacetime_metric->get(alpha, 0) * pi.get(nu, 0);
      for (size_t beta = 1; beta < Dim + 1; ++beta) {
        pi_2_up->get(nu, alpha) +=
            inverse_spacetime_metric->get(alpha, beta) * pi.get(nu, beta);
      }
    }
  }

  // Raise Christoffel
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
        christoffel_first_kind_3_up->get(mu, nu, alpha) =
            inverse_spacetime_metric->get(alpha, 0) *
            christoffel_first_kind->get(mu, nu, 0);
        for (size_t beta = 1; beta < Dim + 1; ++beta) {
          christoffel_first_kind_3_up->get(mu, nu, alpha) +=
              inverse_spacetime_metric->get(alpha, beta) *
              christoffel_first_kind->get(mu, nu, beta);
        }
      }
    }
  }

  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    pi_one_normal->get(mu) = get<0>(*normal_spacetime_vector) * pi.get(0, mu);
    for (size_t nu = 1; nu < Dim + 1; ++nu) {
      pi_one_normal->get(mu) +=
          normal_spacetime_vector->get(nu) * pi.get(nu, mu);
    }
  }

  get(*half_pi_two_normals) =
      get<0>(*normal_spacetime_vector) * get<0>(*pi_one_normal);
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    get(*half_pi_two_normals) +=
        normal_spacetime_vector->get(mu) * pi_one_normal->get(mu);
  }
  get(*half_pi_two_normals) *= 0.5;

  for (size_t n = 0; n < Dim; ++n) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      phi_one_normal->get(n, nu) =
          get<0>(*normal_spacetime_vector) * phi.get(n, 0, nu);
      for (size_t mu = 1; mu < Dim + 1; ++mu) {
        phi_one_normal->get(n, nu) +=
            normal_spacetime_vector->get(mu) * phi.get(n, mu, nu);
      }
    }
  }

  for (size_t n = 0; n < Dim; ++n) {
    half_phi_two_normals->get(n) =
        get<0>(*normal_spacetime_vector) * phi_one_normal->get(n, 0);
    for (size_t mu = 1; mu < Dim + 1; ++mu) {
      half_phi_two_normals->get(n) +=
          normal_spacetime_vector->get(mu) * phi_one_normal->get(n, mu);
    }
    half_phi_two_normals->get(n) *= 0.5;
  }

  for (size_t n = 0; n < Dim; ++n) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        // Mike zero out: n b/c only finite difference in x
        if (n == 0) {
          three_index_constraint->get(n, mu, nu) =
              d_spacetime_metric.get(n, mu, nu) - phi.get(n, mu, nu);
        } else {
          three_index_constraint->get(n, mu, nu) = 0.0;
        }
        // std::cout << n << " " << mu << " " << nu << " " <<
        // three_index_constraint->get(n, mu, nu) << "\n"; std::cout << n << "
        // mu " << mu << " nu " << mu << " d_spacetime_metric.get(n, mu, nu) "
        // << d_spacetime_metric.get(n, mu, nu) << " phi " << phi.get(n, mu, nu)
        // << "\n";
      }
    }
  }

  get(*gamma1_plus_1) = 1.0 + gamma1.get();
  const DataVector& gamma1p1 = get(*gamma1_plus_1);

  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    gauge_constraint->get(mu) = trace_christoffel->get(mu);
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      shift_dot_three_index_constraint->get(mu, nu) =
          get<0>(*shift) * three_index_constraint->get(0, mu, nu);
      if (mesh_velocity.has_value()) {
        mesh_velocity_dot_three_index_constraint->get(mu, nu) =
            get<0>(*mesh_velocity) * three_index_constraint->get(0, mu, nu);
      }
      for (size_t m = 1; m < Dim; ++m) {
        shift_dot_three_index_constraint->get(mu, nu) +=
            shift->get(m) * three_index_constraint->get(m, mu, nu);
        if (mesh_velocity.has_value()) {
          mesh_velocity_dot_three_index_constraint->get(mu, nu) +=
              mesh_velocity->get(m) * three_index_constraint->get(m, mu, nu);
        }
      }
    }
  }

  const bool using_harmonic_gauge = gauge_condition.is_harmonic();
  if (not using_harmonic_gauge) {
    // Compute gauge condition.
    get(*sqrt_det_spatial_metric) = sqrt(get(*det_spatial_metric));
    raise_or_lower_first_index(christoffel_second_kind, *christoffel_first_kind,
                               *inverse_spacetime_metric);
  }
  gauges::dispatch<Dim>(
      gauge_function, spacetime_deriv_gauge_function, *lapse, *shift,
      *sqrt_det_spatial_metric, *inverse_spatial_metric, *da_spacetime_metric,
      *half_pi_two_normals, *half_phi_two_normals, spacetime_metric, phi, mesh,
      time, inertial_coords, inverse_jacobian, gauge_condition);

  // spacetime derivative of killing vector for Y and Z derivatives
  // auto da_killing_vectorY = inverse_spacetime_metric;
  // auto da_killing_vectorY = spacetime_metric;
  tnsr::ab<DataVector, Dim> da_killing_vectorY{};
  auto da_killing_vectorZ = da_killing_vectorY;

  for (size_t a = 0; a < Dim + 1; ++a) {
    for (size_t b = 0; b < Dim + 1; ++b) {
      if (a == 2 and b == 1) {
        da_killing_vectorY.get(a, b) = 0.0 * spacetime_metric.get(a, b) - 1.0;
        da_killing_vectorZ.get(a, b) = 0.0 * spacetime_metric.get(a, b);
      } else if (a == 1 and b == 2) {
        da_killing_vectorY.get(a, b) = 0.0 * spacetime_metric.get(a, b) + 1.0;
        da_killing_vectorZ.get(a, b) = 0.0 * spacetime_metric.get(a, b);
      } else if (a == 3 and b == 1) {
        da_killing_vectorZ.get(a, b) = 0.0 * spacetime_metric.get(a, b) - 1.0;
        da_killing_vectorY.get(a, b) = 0.0 * spacetime_metric.get(a, b);
      } else if (a == 1 and b == 3) {
        da_killing_vectorZ.get(a, b) = 0.0 * spacetime_metric.get(a, b) + 1.0;
        da_killing_vectorY.get(a, b) = 0.0 * spacetime_metric.get(a, b);
      } else {
        da_killing_vectorY.get(a, b) = 0.0 * spacetime_metric.get(a, b);
        da_killing_vectorZ.get(a, b) = 0.0 * spacetime_metric.get(a, b);
      }
    }
  }

  // Zero out Y & Z derivatives of gauge function
  for (size_t a = 0; a < Dim + 1; ++a) {
    for (size_t b = 2; b < 4; ++b) {
      spacetime_deriv_gauge_function->get(b, a) = 0.0;
    }
  }

  // for (size_t dim_index = 1; dim_index < 3; ++dim_index) {
  for (size_t dim_index = 2; dim_index < 4; ++dim_index) {
    for (size_t b = 0; b < Dim + 1; ++b) {
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t c = 0; c < Dim + 1; ++c) {
          for (size_t d = 0; d < Dim + 1; ++d) {
            for (size_t e = 0; e < Dim + 1; ++e) {
              // Gamma' =
              spacetime_deriv_gauge_function->get(dim_index, b) +=
                  christoffel_first_kind->get(d, a, b) *
                  inverse_spacetime_metric->get(d, e) *
                  inverse_spacetime_metric->get(a, c) *
                  phi.get(dim_index - 1, e, c);
            }
          }
        }
      }
    }
  }

  for (size_t b = 0; b < Dim + 1; ++b) {
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t c = 0; c < Dim + 1; ++c) {
        for (size_t d = 0; d < Dim + 1; ++d) {
          // Gamma''
          spacetime_deriv_gauge_function->get(2, b) -=
              1.0 / inertial_coords.get(0) *
              inverse_spacetime_metric->get(d, a) *
              (christoffel_first_kind->get(c, a, b) *
                   da_killing_vectorY.get(d, c) +
               christoffel_first_kind->get(d, c, b) *
                   da_killing_vectorY.get(a, c) +
               christoffel_first_kind->get(d, a, c) *
                   da_killing_vectorY.get(b, c));

          spacetime_deriv_gauge_function->get(3, b) -=
              1.0 / inertial_coords.get(0) *
              inverse_spacetime_metric->get(d, a) *
              (christoffel_first_kind->get(c, a, b) *
                   da_killing_vectorZ.get(d, c) +
               christoffel_first_kind->get(d, c, b) *
                   da_killing_vectorZ.get(a, c) +
               christoffel_first_kind->get(d, a, c) *
                   da_killing_vectorZ.get(b, c));
        }
      }
      std::cout << a << " " << b << " "
                << spacetime_deriv_gauge_function->get(2, b) << "\n";
    }
  }

  // Cartoon deriv gauge function H_a
  // begin with guessing vector transformation
  // next pass is product rule of Lie derivative

  // // dy H_x
  // spacetime_deriv_gauge_function->get(2, 1) =
  //     -gauge_function->get(2) / inertial_coords.get(0);

  // // dy H_y
  // spacetime_deriv_gauge_function->get(2, 2) =
  //     gauge_function->get(1) / inertial_coords.get(0);

  // // dy H_z
  // spacetime_deriv_gauge_function->get(2, 3) = 0.0;

  // // dz H_x
  // spacetime_deriv_gauge_function->get(3, 1) =
  //     -gauge_function->get(3) / inertial_coords.get(0);

  // // dz H_y
  // spacetime_deriv_gauge_function->get(3, 2) = 0.0;

  // // dz H_z
  // spacetime_deriv_gauge_function->get(3, 3) =
  //     gauge_function->get(1) / inertial_coords.get(0);

  if (not using_harmonic_gauge) {
    // Compute source function last so that we don't need to recompute any of
    // the other temporary tags.
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      gauge_constraint->get(nu) += gauge_function->get(nu);
    }
  }

  get(*normal_dot_gauge_constraint) =
      get<0>(*normal_spacetime_vector) * get<0>(*gauge_constraint);
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    get(*normal_dot_gauge_constraint) +=
        normal_spacetime_vector->get(mu) * gauge_constraint->get(mu);
  }

  // Invalidate da_spacetime_metric since we will be modifying some of the
  // data it points to.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  const_cast<std::optional<tnsr::abb<DataVector, Dim>>&>(da_spacetime_metric) =
      std::nullopt;

  // Here are the actual equations

  // Equation for dt_spacetime_metric
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_spacetime_metric->get(mu, nu) +=
          gamma1p1 * shift_dot_three_index_constraint->get(mu, nu);
      // dt_spacetime_metric->get(mu, nu) = 0.0;
      if (mesh_velocity.has_value()) {
        dt_spacetime_metric->get(mu, nu) +=
            get(gamma1) * mesh_velocity_dot_three_index_constraint->get(mu, nu);
      }
    }
  }
  // std::cout << "dt_spacetime_metric->get(0, nu) " <<
  // dt_spacetime_metric->get(0, 0) << " " << dt_spacetime_metric->get(0, 1) <<
  // " " << dt_spacetime_metric->get(0, 2) << " "  <<
  // dt_spacetime_metric->get(0, 3) << " \n"; std::cout <<
  // "dt_spacetime_metric->get(1, nu) " << dt_spacetime_metric->get(1, 0) << " "
  // << dt_spacetime_metric->get(1, 1) << " " << dt_spacetime_metric->get(1, 2)
  // << " "  << dt_spacetime_metric->get(1, 3) << " \n"; std::cout <<
  // "dt_spacetime_metric->get(2, nu) " << dt_spacetime_metric->get(2, 0) << " "
  // << dt_spacetime_metric->get(2, 1) << " " << dt_spacetime_metric->get(2, 2)
  // << " "  << dt_spacetime_metric->get(2, 3) << " \n"; std::cout <<
  // "dt_spacetime_metric->get(3, nu) " << dt_spacetime_metric->get(3, 0) << " "
  // << dt_spacetime_metric->get(3, 1) << " " << dt_spacetime_metric->get(3, 2)
  // << " "  << dt_spacetime_metric->get(3, 3) << " \n";

  // Equation for dt_pi

  // We first compute the n_a contributions but only for a=0 since n_i=0
  // identically. We also use dt_Pi_{00} as temporary storage to avoid any extra
  // allocations and multiply normal_dot_gauge_constraint=(n^a C_a) by gamma0
  // since it always shows up multiplied by gamma0 in the equations. This
  // reduces the number of multiplications that are needed for the RHS
  // evaluation.
  //
  // WARNING: normal_dot_gauge_constraint is rescaled by gamma0!
  get(*normal_dot_gauge_constraint) *= get(gamma0);

  // Use dt_pi_{00} as temporary storage.
  get<0, 0>(*dt_pi) = -get(gamma0) * get(*lapse);
  for (size_t i = 1; i < Dim + 1; ++i) {
    dt_pi->get(0, i) =
        get<0, 0>(*dt_pi) * gauge_constraint->get(i) -
        get(*normal_dot_gauge_constraint) * spacetime_metric.get(0, i);
  }
  get<0, 0>(*dt_pi) =
      2.0 * get<0, 0>(*dt_pi) * get<0>(*gauge_constraint) -
      get(*normal_dot_gauge_constraint) * get<0, 0>(spacetime_metric);

  // Set space-space components
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_pi->get(mu, nu) =
          -get(*normal_dot_gauge_constraint) * spacetime_metric.get(mu, nu);
    }
  }

  // Add additional pieces to dt_pi that aren't just n_a*(stuff)
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_pi->get(mu, nu) -= get(*half_pi_two_normals) * pi.get(mu, nu);

      if (not using_harmonic_gauge) {
        dt_pi->get(mu, nu) -= spacetime_deriv_gauge_function->get(mu, nu) +
                              spacetime_deriv_gauge_function->get(nu, mu);
      }
      for (size_t delta = 0; delta < Dim + 1; ++delta) {
        dt_pi->get(mu, nu) -= 2 * pi.get(mu, delta) * pi_2_up->get(nu, delta);
        if (not using_harmonic_gauge) {
          dt_pi->get(mu, nu) += 2 *
                                christoffel_second_kind->get(delta, mu, nu) *
                                gauge_function->get(delta);
        }
        for (size_t n = 0; n < Dim; ++n) {
          dt_pi->get(mu, nu) +=
              2 * phi_1_up->get(n, mu, delta) * phi_3_up->get(n, nu, delta);
        }

        for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
          dt_pi->get(mu, nu) -=
              2. * christoffel_first_kind_3_up->get(mu, alpha, delta) *
              christoffel_first_kind_3_up->get(nu, delta, alpha);
        }
      }

      for (size_t m = 0; m < Dim; ++m) {
        dt_pi->get(mu, nu) -=
            pi_one_normal->get(m + 1) * phi_1_up->get(m, mu, nu);

        for (size_t n = 0; n < Dim; ++n) {
          dt_pi->get(mu, nu) -=
              inverse_spatial_metric->get(m, n) * d_phi.get(m, n, mu, nu);
          // std::cout << "d_phi = d_m d_n g_mu_nu " << m << " " << n << " " <<
          // mu
          //           << " " << nu << " " << d_phi.get(m, n, mu, nu) << "\n";
        }
      }

      dt_pi->get(mu, nu) *= get(*lapse);

      dt_pi->get(mu, nu) +=
          gamma12 * shift_dot_three_index_constraint->get(mu, nu);
      if (mesh_velocity.has_value()) {
        dt_pi->get(mu, nu) +=
            gamma12 * mesh_velocity_dot_three_index_constraint->get(mu, nu);
      }

      for (size_t m = 0; m < Dim; ++m) {
        // DualFrame term
        dt_pi->get(mu, nu) += shift->get(m) * d_pi.get(m, mu, nu);
      }
    }
  }

  // Equation for dt_phi
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        dt_phi->get(i, mu, nu) =
            pi.get(mu, nu) * half_phi_two_normals->get(i) -
            d_pi.get(i, mu, nu) +
            get(gamma2) * three_index_constraint->get(i, mu, nu);
        for (size_t n = 0; n < Dim; ++n) {
          dt_phi->get(i, mu, nu) +=
              phi_one_normal->get(i, n + 1) * phi_1_up->get(n, mu, nu);
        }

        dt_phi->get(i, mu, nu) *= get(*lapse);
        for (size_t m = 0; m < Dim; ++m) {
          dt_phi->get(i, mu, nu) += shift->get(m) * d_phi.get(m, i, mu, nu);
        }
      }
    }
  }
}
}  // namespace gh

// Explicit instantiations of structs defined in `Equations.cpp` as well as of
// `partial_derivatives` function for use in the computation of spatial
// derivatives of `gradients_tags`, and of the initial gauge source function
// (needed in `Initialize.hpp`).
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"

using derivative_frame = Frame::Inertial;

template <size_t Dim>
using derivative_tags_initial_gauge =
    tmpl::list<gh::Tags::InitialGaugeH<DataVector, Dim, derivative_frame>>;

template <size_t Dim>
using variables_tags_initial_gauge =
    tmpl::list<gh::Tags::InitialGaugeH<DataVector, Dim, derivative_frame>>;

template <size_t Dim>
using derivative_tags = typename gh::System<Dim>::gradients_tags;

template <size_t Dim>
using variables_tags = typename gh::System<Dim>::variables_tag::tags_list;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define INSTANTIATE(_, data)                                                 \
  template struct gh::TimeDerivative<DIM(data)>;                             \
  template Variables<                                                        \
      db::wrap_tags_in<::Tags::deriv, derivative_tags<DIM(data)>,            \
                       tmpl::size_t<DIM(data)>, derivative_frame>>           \
  partial_derivatives<derivative_tags<DIM(data)>, variables_tags<DIM(data)>, \
                      DIM(data), derivative_frame>(                          \
      const Variables<variables_tags<DIM(data)>>& u,                         \
      const Mesh<DIM(data)>& mesh,                                           \
      const InverseJacobian<DataVector, DIM(data), Frame::ElementLogical,    \
                            derivative_frame>& inverse_jacobian);            \
  template Variables<db::wrap_tags_in<                                       \
      ::Tags::deriv, derivative_tags_initial_gauge<DIM(data)>,               \
      tmpl::size_t<DIM(data)>, derivative_frame>>                            \
  partial_derivatives<derivative_tags_initial_gauge<DIM(data)>,              \
                      variables_tags_initial_gauge<DIM(data)>, DIM(data),    \
                      derivative_frame>(                                     \
      const Variables<variables_tags_initial_gauge<DIM(data)>>& u,           \
      const Mesh<DIM(data)>& mesh,                                           \
      const InverseJacobian<DataVector, DIM(data), Frame::ElementLogical,    \
                            derivative_frame>& inverse_jacobian);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
