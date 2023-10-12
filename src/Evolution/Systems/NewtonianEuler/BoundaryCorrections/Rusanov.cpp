// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/NewtonianEuler/BoundaryCorrections/Rusanov.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/NewtonianEuler/SoundSpeedSquared.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/NormalDotFlux.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace NewtonianEuler::BoundaryCorrections {
template <size_t Dim>
Rusanov<Dim>::Rusanov(CkMigrateMessage* msg) : BoundaryCorrection<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<BoundaryCorrection<Dim>> Rusanov<Dim>::get_clone() const {
  return std::make_unique<Rusanov>(*this);
}

template <size_t Dim>
void Rusanov<Dim>::pup(PUP::er& p) {
  BoundaryCorrection<Dim>::pup(p);
}

template <size_t Dim>
template <size_t ThermodynamicDim>
double Rusanov<Dim>::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> packaged_mass_density,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        packaged_momentum_density,
    const gsl::not_null<Scalar<DataVector>*> packaged_energy_density,
    const gsl::not_null<Scalar<DataVector>*>
        packaged_normal_dot_flux_mass_density,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        packaged_normal_dot_flux_momentum_density,
    const gsl::not_null<Scalar<DataVector>*>
        packaged_normal_dot_flux_energy_density,
    const gsl::not_null<Scalar<DataVector>*> packaged_abs_char_speed,

    const Scalar<DataVector>& mass_density,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& momentum_density,
    const Scalar<DataVector>& energy_density,

    const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_mass_density,
    const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_momentum_density,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_energy_density,

    const tnsr::I<DataVector, Dim, Frame::Inertial>& velocity,
    const Scalar<DataVector>& specific_internal_energy,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity,
    const EquationsOfState::EquationOfState<false, ThermodynamicDim>&
        equation_of_state) const {
  {
    // Compute sound speed
    Scalar<DataVector>& sound_speed = *packaged_mass_density;
    Scalar<DataVector>& normal_dot_velocity = *packaged_energy_density;
    sound_speed_squared(make_not_null(&sound_speed), mass_density,
                        specific_internal_energy, equation_of_state);
    get(sound_speed) = sqrt(get(sound_speed));
    dot_product(make_not_null(&normal_dot_velocity), velocity, normal_covector);
    if (normal_dot_mesh_velocity.has_value()) {
      get(*packaged_abs_char_speed) =
          max(abs(get(normal_dot_velocity) + get(sound_speed) -
                  get(*normal_dot_mesh_velocity)),
              abs(get(normal_dot_velocity) - get(sound_speed) -
                  get(*normal_dot_mesh_velocity)));
    } else {
      get(*packaged_abs_char_speed) =
          max(abs(get(normal_dot_velocity) + get(sound_speed)),
              abs(get(normal_dot_velocity) - get(sound_speed)));
    }
  }

  *packaged_mass_density = mass_density;
  *packaged_momentum_density = momentum_density;
  *packaged_energy_density = energy_density;

  normal_dot_flux(packaged_normal_dot_flux_mass_density, normal_covector,
                  flux_mass_density);
  normal_dot_flux(packaged_normal_dot_flux_momentum_density, normal_covector,
                  flux_momentum_density);
  normal_dot_flux(packaged_normal_dot_flux_energy_density, normal_covector,
                  flux_energy_density);

  return max(get(*packaged_abs_char_speed));
}

template <size_t Dim>
void Rusanov<Dim>::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_mass_density,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        boundary_correction_momentum_density,
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_energy_density,
    const Scalar<DataVector>& mass_density_int,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& momentum_density_int,
    const Scalar<DataVector>& energy_density_int,
    const Scalar<DataVector>& normal_dot_flux_mass_density_int,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        normal_dot_flux_momentum_density_int,
    const Scalar<DataVector>& normal_dot_flux_energy_density_int,
    const Scalar<DataVector>& abs_char_speed_int,
    const Scalar<DataVector>& mass_density_ext,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& momentum_density_ext,
    const Scalar<DataVector>& energy_density_ext,
    const Scalar<DataVector>& normal_dot_flux_mass_density_ext,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        normal_dot_flux_momentum_density_ext,
    const Scalar<DataVector>& normal_dot_flux_energy_density_ext,
    const Scalar<DataVector>& abs_char_speed_ext,
    const dg::Formulation dg_formulation) const {
  if (dg_formulation == dg::Formulation::WeakInertial) {
    get(*boundary_correction_mass_density) =
        0.5 * (get(normal_dot_flux_mass_density_int) -
               get(normal_dot_flux_mass_density_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(mass_density_ext) - get(mass_density_int));
    for (size_t i = 0; i < Dim; ++i) {
      boundary_correction_momentum_density->get(i) =
          0.5 * (normal_dot_flux_momentum_density_int.get(i) -
                 normal_dot_flux_momentum_density_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (momentum_density_ext.get(i) - momentum_density_int.get(i));
    }
    get(*boundary_correction_energy_density) =
        0.5 * (get(normal_dot_flux_energy_density_int) -
               get(normal_dot_flux_energy_density_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(energy_density_ext) - get(energy_density_int));
  } else {
    get(*boundary_correction_mass_density) =
        -0.5 * (get(normal_dot_flux_mass_density_int) +
                get(normal_dot_flux_mass_density_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(mass_density_ext) - get(mass_density_int));
    for (size_t i = 0; i < Dim; ++i) {
      boundary_correction_momentum_density->get(i) =
          -0.5 * (normal_dot_flux_momentum_density_int.get(i) +
                  normal_dot_flux_momentum_density_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (momentum_density_ext.get(i) - momentum_density_int.get(i));
    }
    get(*boundary_correction_energy_density) =
        -0.5 * (get(normal_dot_flux_energy_density_int) +
                get(normal_dot_flux_energy_density_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(energy_density_ext) - get(energy_density_int));
  }
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID Rusanov<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(_, data) template class Rusanov<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION

#define THERMODIM(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATION(_, data)                                                 \
  template double Rusanov<DIM(data)>::dg_package_data<THERMODIM(data)>(        \
      gsl::not_null<Scalar<DataVector>*> packaged_mass_density,                \
      gsl::not_null<tnsr::I<DataVector, DIM(data), Frame::Inertial>*>          \
          packaged_momentum_density,                                           \
      gsl::not_null<Scalar<DataVector>*> packaged_energy_density,              \
      gsl::not_null<Scalar<DataVector>*>                                       \
          packaged_normal_dot_flux_mass_density,                               \
      gsl::not_null<tnsr::I<DataVector, DIM(data), Frame::Inertial>*>          \
          packaged_normal_dot_flux_momentum_density,                           \
      gsl::not_null<Scalar<DataVector>*>                                       \
          packaged_normal_dot_flux_energy_density,                             \
      gsl::not_null<Scalar<DataVector>*> packaged_abs_char_speed,              \
                                                                               \
      const Scalar<DataVector>& mass_density,                                  \
      const tnsr::I<DataVector, DIM(data), Frame::Inertial>& momentum_density, \
      const Scalar<DataVector>& energy_density,                                \
                                                                               \
      const tnsr::I<DataVector, DIM(data), Frame::Inertial>&                   \
          flux_mass_density,                                                   \
      const tnsr::IJ<DataVector, DIM(data), Frame::Inertial>&                  \
          flux_momentum_density,                                               \
      const tnsr::I<DataVector, DIM(data), Frame::Inertial>&                   \
          flux_energy_density,                                                 \
                                                                               \
      const tnsr::I<DataVector, DIM(data), Frame::Inertial>& velocity,         \
      const Scalar<DataVector>& specific_internal_energy,                      \
                                                                               \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& normal_covector,  \
      const std::optional<tnsr::I<DataVector, DIM(data), Frame::Inertial>>&    \
          mesh_velocity,                                                       \
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity,       \
      const EquationsOfState::EquationOfState<false, THERMODIM(data)>&         \
          equation_of_state) const;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3), (1, 2))

#undef INSTANTIATION
#undef THERMODIM
#undef DIM
}  // namespace NewtonianEuler::BoundaryCorrections
