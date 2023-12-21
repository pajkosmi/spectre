// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/DgSubcell/GhostZoneLogicalCoordinates.hpp"
#include "Evolution/DgSubcell/Tags/Coordinates.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/FiniteDifference/Factory.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/FiniteDifference/Reconstructor.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/FiniteDifference/Tag.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/Tags.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/BoundaryConditions/HydroFreeOutflow.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/ConservativeFromPrimitive.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Fluxes.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
/// \endcond

namespace grmhd::GhValenciaDivClean::BoundaryConditions {
/*!
 * \brief Sets Dirichlet boundary conditions using the analytic solution or
 * analytic data on the spacetime variables and hydro free outflow on the GRMHD
 * variables.
 */
class DirichletFreeOutflow final : public BoundaryCondition {
 private:
  template <typename T>
  using Flux = ::Tags::Flux<T, tmpl::size_t<3>, Frame::Inertial>;
  std::unique_ptr<evolution::initial_data::InitialData> analytic_prescription_;

 public:
  /// \brief What analytic solution/data to prescribe.
  struct AnalyticPrescription {
    static constexpr Options::String help =
        "What analytic solution/data to prescribe.";
    using type = std::unique_ptr<evolution::initial_data::InitialData>;
  };
  using options = tmpl::list<AnalyticPrescription>;
  static constexpr Options::String help{
      "DirichletFreeOutflow boundary conditions using either analytic solution "
      "or analytic data for GH variables and hydro free outflow for GRMHD"};

  DirichletFreeOutflow() = default;
  DirichletFreeOutflow(DirichletFreeOutflow&&) = default;
  DirichletFreeOutflow& operator=(DirichletFreeOutflow&&) = default;
  DirichletFreeOutflow(const DirichletFreeOutflow&);
  DirichletFreeOutflow& operator=(const DirichletFreeOutflow&);
  ~DirichletFreeOutflow() override = default;

  explicit DirichletFreeOutflow(CkMigrateMessage* msg);

  explicit DirichletFreeOutflow(
      std::unique_ptr<evolution::initial_data::InitialData>
          analytic_prescription);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, DirichletFreeOutflow);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags =
      tmpl::list<domain::Tags::Coordinates<3, Frame::Inertial>,
                 ::gh::ConstraintDamping::Tags::ConstraintGamma1,
                 ::gh::ConstraintDamping::Tags::ConstraintGamma2>;
  using dg_interior_primitive_variables_tags =
      tmpl::list<hydro::Tags::RestMassDensity<DataVector>,
                 hydro::Tags::ElectronFraction<DataVector>,
                 hydro::Tags::SpecificInternalEnergy<DataVector>,
                 hydro::Tags::SpatialVelocity<DataVector, 3>,
                 hydro::Tags::MagneticField<DataVector, 3>,
                 hydro::Tags::LorentzFactor<DataVector>,
                 hydro::Tags::Pressure<DataVector>>;
  using dg_gridless_tags = tmpl::list<::Tags::Time>;
  std::optional<std::string> dg_ghost(
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
      gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> phi,
      gsl::not_null<Scalar<DataVector>*> tilde_d,
      gsl::not_null<Scalar<DataVector>*> tilde_ye,
      gsl::not_null<Scalar<DataVector>*> tilde_tau,
      gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> tilde_s,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> tilde_b,
      gsl::not_null<Scalar<DataVector>*> tilde_phi,

      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> tilde_d_flux,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> tilde_ye_flux,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> tilde_tau_flux,
      gsl::not_null<tnsr::Ij<DataVector, 3, Frame::Inertial>*> tilde_s_flux,
      gsl::not_null<tnsr::IJ<DataVector, 3, Frame::Inertial>*> tilde_b_flux,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> tilde_phi_flux,

      gsl::not_null<Scalar<DataVector>*> gamma1,
      gsl::not_null<Scalar<DataVector>*> gamma2,
      gsl::not_null<Scalar<DataVector>*> lapse,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
      gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
          inv_spatial_metric,

      const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,

      const Scalar<DataVector>& interior_rest_mass_density,
      const Scalar<DataVector>& interior_electron_fraction,
      const Scalar<DataVector>& interior_specific_internal_energy,
      const tnsr::I<DataVector, 3, Frame::Inertial>& interior_spatial_velocity,
      const tnsr::I<DataVector, 3, Frame::Inertial>& interior_magnetic_field,
      const Scalar<DataVector>& interior_lorentz_factor,
      const Scalar<DataVector>& interior_pressure,

      const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
      const Scalar<DataVector>& interior_gamma1,
      const Scalar<DataVector>& interior_gamma2, double time) const;

  using fd_interior_evolved_variables_tags = tmpl::list<>;
  using fd_interior_temporary_tags =
      tmpl::list<evolution::dg::subcell::Tags::Mesh<3>>;
  using fd_interior_primitive_variables_tags =
      tmpl::list<hydro::Tags::RestMassDensity<DataVector>,
                 hydro::Tags::ElectronFraction<DataVector>,
                 hydro::Tags::Temperature<DataVector>,
                 hydro::Tags::Pressure<DataVector>,
                 hydro::Tags::SpecificInternalEnergy<DataVector>,
                 hydro::Tags::LorentzFactor<DataVector>,
                 hydro::Tags::SpatialVelocity<DataVector, 3>,
                 hydro::Tags::MagneticField<DataVector, 3>,
                 gr::Tags::SpacetimeMetric<DataVector, 3>,
                 ::gh::Tags::Pi<DataVector, 3>, ::gh::Tags::Phi<DataVector, 3>>;
  using fd_gridless_tags =
      tmpl::list<::Tags::Time, ::domain::Tags::FunctionsOfTime,
                 domain::Tags::ElementMap<3, Frame::Grid>,
                 domain::CoordinateMaps::Tags::CoordinateMap<3, Frame::Grid,
                                                             Frame::Inertial>,
                 fd::Tags::Reconstructor>;
  void fd_ghost(

      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
      gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> phi,
      gsl::not_null<Scalar<DataVector>*> rest_mass_density,
      gsl::not_null<Scalar<DataVector>*> electron_fraction,
      gsl::not_null<Scalar<DataVector>*> temperature,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
          lorentz_factor_times_spatial_velocity,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> magnetic_field,
      gsl::not_null<Scalar<DataVector>*> divergence_cleaning_field,
      const Direction<3>& direction,

      // fd_interior_temporary_tags
      const Mesh<3>& subcell_mesh,

      // interior prim vars tags
      const Scalar<DataVector>& interior_rest_mass_density,
      const Scalar<DataVector>& interior_electron_fraction,
      const Scalar<DataVector>& interior_temperature,
      const Scalar<DataVector>& interior_pressure,
      const Scalar<DataVector>& interior_specific_internal_energy,
      const Scalar<DataVector>& interior_lorentz_factor,
      const tnsr::I<DataVector, 3, Frame::Inertial>& interior_spatial_velocity,
      const tnsr::I<DataVector, 3, Frame::Inertial>& interior_magnetic_field,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& interior_spacetime_metric,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& interior_pi,
      const tnsr::iaa<DataVector, 3, Frame::Inertial>& interior_phi,

      // fd_gridless_tags
      double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<::domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time,
      const ElementMap<3, Frame::Grid>& logical_to_grid_map,
      const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>&
          grid_to_inertial_map,
      const fd::Reconstructor& reconstructor,
      const AnalyticSolutionOrData& analytic_solution_or_data) const {
    const size_t ghost_zone_size{reconstructor.ghost_zone_size()};

    const auto ghost_logical_coords =
        evolution::dg::subcell::fd::ghost_zone_logical_coordinates(
            subcell_mesh, ghost_zone_size, direction);

    const auto ghost_inertial_coords = grid_to_inertial_map(
        logical_to_grid_map(ghost_logical_coords), time, functions_of_time);

    using SpacetimeMetric = gr::Tags::SpacetimeMetric<DataVector, 3>;
    using Pi = gh::Tags::Pi<DataVector, 3>;
    using Phi = gh::Tags::Phi<DataVector, 3>;
    using spacetime_tags = tmpl::list<SpacetimeMetric, Pi, Phi>;

    const size_t buffer_size_per_grid_pts =
        Variables<spacetime_tags>::number_of_independent_components;

    const size_t dim_direction{direction.dimension()};

    const auto subcell_extents{subcell_mesh.extents()};

    const size_t num_face_pts{
        subcell_extents.slice_away(dim_direction).product()};

    DataVector buffer_for_vars{
        num_face_pts * ((1 + ghost_zone_size) * (buffer_size_per_grid_pts)),
        0.0};

    Variables<spacetime_tags> outermost_prim_vars{
        buffer_for_vars.data(), num_face_pts * buffer_size_per_grid_pts};
    Variables<spacetime_tags> ghost_prim_vars{
        outermost_prim_vars.data() + outermost_prim_vars.size(),
        num_face_pts * buffer_size_per_grid_pts * ghost_zone_size};

    auto get_boundary_val = [&direction, &subcell_extents](auto volume_tensor) {
      return evolution::dg::subcell::slice_tensor_for_subcell(
          volume_tensor, subcell_extents, 1, direction, {});
    };

    // ensuring derivative of spacetime metric at origin is 0
    get<SpacetimeMetric>(outermost_prim_vars) =
        get_boundary_val(interior_spacetime_metric);
    // ensuring derivative of Pi at origin is 0.  Mike: should this happen?
    get<Pi>(outermost_prim_vars) = get_boundary_val(interior_pi);

    // Phi is antisymmetric.
    for (size_t i = 0; i < 3; i++) {
      for (size_t a = 0; a < 4; a++) {
        for (size_t b = 0; b < 4; b++) {
          get<Phi>(outermost_prim_vars).get(i, a, b) =
              -1.0 * get_boundary_val(interior_phi).get(i, a, b);
        }
      }
    }

    // Now copy `outermost_prim_vars` into each slices of `ghost_prim_vars`.
    Index<3> ghost_data_extents = subcell_extents;
    ghost_data_extents[dim_direction] = ghost_zone_size;

    for (size_t i_ghost = 0; i_ghost < ghost_zone_size; ++i_ghost) {
      add_slice_to_data(make_not_null(&ghost_prim_vars), outermost_prim_vars,
                        ghost_data_extents, dim_direction, i_ghost);
    }

    // move data from buffer to guard cells
    *spacetime_metric = get<SpacetimeMetric>(ghost_prim_vars);
    *pi = get<Pi>(ghost_prim_vars);
    *phi = get<Phi>(ghost_prim_vars);

    // Mike: apply hydrofree outflow logic to spacetime variables instead of
    // picking analytic solutions Compute FD ghost data with the analytic data
    // or solution
    // auto boundary_values = [&analytic_solution_or_data,
    // &ghost_inertial_coords,
    //                         &time]() {
    //   if constexpr (std::is_base_of_v<MarkAsAnalyticData,
    //                                   AnalyticSolutionOrData>) {
    //     (void)time;
    //     return analytic_solution_or_data.variables(
    //         ghost_inertial_coords,
    //         tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
    //                    ::gh::Tags::Pi<DataVector, 3>,
    //                    ::gh::Tags::Phi<DataVector, 3>>{});
    //   } else {
    //     return analytic_solution_or_data.variables(
    //         ghost_inertial_coords, time,
    //         tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
    //                    ::gh::Tags::Pi<DataVector, 3>,
    //                    ::gh::Tags::Phi<DataVector, 3>>{});
    //   }
    // }();

    // MIKE: BCs
    // *spacetime_metric =
    //     get<gr::Tags::SpacetimeMetric<DataVector, 3>>(boundary_values);
    // *pi = get<::gh::Tags::Pi<DataVector, 3>>(boundary_values);
    // *phi = get<::gh::Tags::Phi<DataVector, 3>>(boundary_values);

    // Note: Once we support high-order fluxes with GHMHD we will need to
    // handle this correctly.
    std::optional<Variables<db::wrap_tags_in<
        Flux, typename grmhd::ValenciaDivClean::System::flux_variables>>>
        cell_centered_ghost_fluxes{std::nullopt};
    // Set to zero since it shouldn't be used
    Scalar<DataVector> pressure{};
    Scalar<DataVector> specific_internal_energy{};
    tnsr::I<DataVector, 3> spatial_velocity{};
    Scalar<DataVector> lorentz_factor{};
    const tnsr::I<DataVector, 3> interior_shift{};
    const Scalar<DataVector> interior_lapse{};
    const tnsr::ii<DataVector, 3> interior_spatial_metric{};
    tnsr::ii<DataVector, 3> spatial_metric{};
    tnsr::II<DataVector, 3> inv_spatial_metric{};
    Scalar<DataVector> sqrt_det_spatial_metric{};
    Scalar<DataVector> lapse{};
    tnsr::I<DataVector, 3> shift{};

    grmhd::ValenciaDivClean::BoundaryConditions::HydroFreeOutflow::
        fd_ghost_impl(
            rest_mass_density, electron_fraction, temperature,
            make_not_null(&pressure), make_not_null(&specific_internal_energy),
            lorentz_factor_times_spatial_velocity,
            make_not_null(&spatial_velocity), make_not_null(&lorentz_factor),
            magnetic_field, divergence_cleaning_field,

            make_not_null(&spatial_metric), make_not_null(&inv_spatial_metric),
            make_not_null(&sqrt_det_spatial_metric), make_not_null(&lapse),
            make_not_null(&shift),

            direction,

            // fd_interior_temporary_tags
            subcell_mesh,

            // fd_interior_primitive_variables_tags
            interior_rest_mass_density, interior_electron_fraction,
            interior_temperature, interior_pressure,
            interior_specific_internal_energy, interior_lorentz_factor,
            interior_spatial_velocity, interior_magnetic_field,
            // Note: metric vars are empty because they shouldn't be used
            interior_spatial_metric, interior_lapse, interior_shift,

            // fd_gridless_tags
            reconstructor.ghost_zone_size(),
            cell_centered_ghost_fluxes.has_value());
  }
};
}  // namespace grmhd::GhValenciaDivClean::BoundaryConditions
