// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/IdealFluid.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma:  no_include <pup.h>

/// \cond
namespace PUP {
class er;  // IWYU pragma: keep
}  // namespace PUP
/// \endcond

namespace grmhd::AnalyticData {

/*!
 * \brief Analytic initial data for a cylindrical or spherical blast wave.
 *
 * This class implements analytic initial data for a cylindrical blast wave,
 * as described, e.g., in \cite Kidder2016hev Sec. 6.2.3.
 * A uniform magnetic field threads an ideal fluid. The solution begins with
 * material at fixed (typically high) density and pressure at rest inside a
 * cylinder of radius \f$r < r_{\rm in}\f$ and material at fixed (typically low)
 * density and pressure at rest in a cylindrical shell with radius
 * \f$r > r_{\rm out}\f$. In the region \f$ r_{\rm in} < r < r_{\rm out}\f$,
 * the solution transitions such that the logarithms of the density and
 * pressure vary linearly. E.g., if \f$\rho(r < r_{\rm in}) = \rho_{\rm in}\f$
 * and \f$\rho(r > r_{\rm out}) = \rho_{\rm out}\f$, then
 * \f[
 * \log \rho = [(r_{\rm in} - r) \log(\rho_{\rm out})
 *              + (r - r_{\rm out}) \log(\rho_{\rm in})]
 *              / (r_{\rm in} - r_{\rm out}).
 * \f]
 * Note that the cylinder's axis is the \f$z\f$ axis. To evolve this analytic
 * initial data, use a cubic or cylindrical domain with periodic boundary
 * conditions applied to the outer boundaries whose normals are parallel or
 * antiparallel to the z axis. In the transverse (e.g., x and y) dimensions, the
 * domain should be large enough that the blast wave doesn't reach the boundary
 * at the final time. E.g., if `InnerRadius = 0.8`, `OuterRadius = 1.0`, and
 * the final time is 4.0, a good domain extends from `(x,y)=(-6.0, -6.0)` to
 * `(x,y)=(6.0, 6.0)`.
 *
 * An analogous problem with spherical geometry has also been used
 * \cite CerdaDuran2007 \cite CerdaDuran2008 \cite Cipolletta2019. The magnetic
 * field is chosen to be in the z-direction instead of the x-direction.
 */
class BlastWave : public evolution::initial_data::InitialData,
                  public MarkAsAnalyticData,
                  public AnalyticDataBase {
 public:
  using equation_of_state_type = EquationsOfState::IdealFluid<true>;

  enum class Geometry { Cylindrical, Spherical };

  /// Inside InnerRadius, density is InnerDensity.
  struct InnerRadius {
    using type = double;
    static constexpr Options::String help = {
        "Inside InnerRadius, density is InnerDensity."};
    static type lower_bound() { return 0.0; }
  };
  /// Outside OuterRadius, density is OuterDensity.
  struct OuterRadius {
    using type = double;
    static constexpr Options::String help = {
        "Outside OuterRadius, density is OuterDensity."};
    static type lower_bound() { return 0.0; }
  };
  /// Density at radii less than InnerRadius.
  struct InnerDensity {
    using type = double;
    static constexpr Options::String help = {
        "Density at radii less than InnerRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// Density at radii greater than OuterRadius.
  struct OuterDensity {
    using type = double;
    static constexpr Options::String help = {
        "Density at radii greater than OuterRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// Pressure at radii less than InnerRadius.
  struct InnerPressure {
    using type = double;
    static constexpr Options::String help = {
        "Pressure at radii less than InnerRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// Pressure at radii greater than OuterRadius.
  struct OuterPressure {
    using type = double;
    static constexpr Options::String help = {
        "Pressure at radii greater than OuterRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// The x,y,z components of the uniform magnetic field threading the matter.
  struct MagneticField {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "The x,y,z components of the uniform magnetic field."};
  };
  /// The adiabatic index of the ideal fluid.
  struct AdiabaticIndex {
    using type = double;
    static constexpr Options::String help = {
        "The adiabatic index of the ideal fluid."};
    static type lower_bound() { return 1.0; }
  };
  /// The geometry of the blast wave, i.e. Cylindrical or Spherical.
  struct GeometryOption {
    static std::string name() { return "Geometry"; }
    using type = Geometry;
    static constexpr Options::String help = {
        "The geometry of the blast wave, i.e. Cylindrical or Spherical."};
  };

  /// Radial velocity of the fluid
  struct RadialVelocity {
    using type = double;
    static constexpr Options::String help = {
        "The radial velocity of the ideal fluid."};
    //static type lower_bound() { return -2.989e10; }
  };

  /// Radial velocity fraction of the fluid
  struct RadialVelocityFraction {
    using type = double;
    static constexpr Options::String help = {
      "The inner radius, as a fraction of the outer radius"};
    static type lower_bound() { return 0; }
  };

  using options = tmpl::list<InnerRadius, OuterRadius, InnerDensity,
                             OuterDensity, InnerPressure, OuterPressure,
                             MagneticField, AdiabaticIndex, RadialVelocity,
                             RadialVelocityFraction, GeometryOption>;

  static constexpr Options::String help = {
      "Cylindrical or spherical blast wave analytic initial data."};

  BlastWave() = default;
  BlastWave(const BlastWave& /*rhs*/) = default;
  BlastWave& operator=(const BlastWave& /*rhs*/) = default;
  BlastWave(BlastWave&& /*rhs*/) = default;
  BlastWave& operator=(BlastWave&& /*rhs*/) = default;
  ~BlastWave() override = default;

  BlastWave(double inner_radius, double outer_radius, double inner_density,
            double outer_density, double inner_pressure, double outer_pressure,
            const std::array<double, 3>& magnetic_field, double adiabatic_index,
            double radial_velocity, double radial_velocity_fraction,
            Geometry geometry, const Options::Context& context = {});

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit BlastWave(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(BlastWave);
  /// \endcond

  /// @{
  /// Retrieve the GRMHD variables at a given position.
  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::RestMassDensity<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::RestMassDensity<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::ElectronFraction<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::ElectronFraction<DataType>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x,
      tmpl::list<hydro::Tags::SpecificInternalEnergy<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::SpecificInternalEnergy<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::Pressure<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::Pressure<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::SpatialVelocity<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataType, 3>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::MagneticField<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::MagneticField<DataType, 3>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x,
      tmpl::list<hydro::Tags::DivergenceCleaningField<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::DivergenceCleaningField<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::LorentzFactor<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::SpecificEnthalpy<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpecificEnthalpy<DataType>>;
  /// @}

  /// Retrieve a collection of hydrodynamic variables at position x
  template <typename DataType, typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataType, 3>& x,
                                         tmpl::list<Tags...> /*meta*/) const {
    static_assert(sizeof...(Tags) > 1,
                  "The generic template will recurse infinitely if only one "
                  "tag is being retrieved.");
    return {tuples::get<Tags>(variables(x, tmpl::list<Tags>{}))...};
  }

  /// Retrieve the metric variables
  template <typename DataType, typename Tag>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataType, 3>& x,
                                     tmpl::list<Tag> /*meta*/) const {
    constexpr double dummy_time = 0.0;
    return background_spacetime_.variables(x, dummy_time, tmpl::list<Tag>{});
  }

  const EquationsOfState::IdealFluid<true>& equation_of_state() const {
    return equation_of_state_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 private:
  double inner_radius_ = std::numeric_limits<double>::signaling_NaN();
  double outer_radius_ = std::numeric_limits<double>::signaling_NaN();
  double inner_density_ = std::numeric_limits<double>::signaling_NaN();
  double outer_density_ = std::numeric_limits<double>::signaling_NaN();
  double inner_pressure_ = std::numeric_limits<double>::signaling_NaN();
  double outer_pressure_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 3> magnetic_field_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  double adiabatic_index_ = std::numeric_limits<double>::signaling_NaN();
  double radial_velocity_ = std::numeric_limits<double>::signaling_NaN();
  double radial_velocity_fraction_ =
  std::numeric_limits<double>::signaling_NaN();
  Geometry geometry_ = Geometry::Cylindrical;
  EquationsOfState::IdealFluid<true> equation_of_state_{};
  gr::Solutions::Minkowski<3> background_spacetime_{};

  friend bool operator==(const BlastWave& lhs, const BlastWave& rhs);

  friend bool operator!=(const BlastWave& lhs, const BlastWave& rhs);
};

}  // namespace grmhd::AnalyticData

/// \cond
template <>
struct Options::create_from_yaml<grmhd::AnalyticData::BlastWave::Geometry> {
  template <typename Metavariables>
  static grmhd::AnalyticData::BlastWave::Geometry create(
      const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
grmhd::AnalyticData::BlastWave::Geometry
Options::create_from_yaml<grmhd::AnalyticData::BlastWave::Geometry>::create<
    void>(const Options::Option& options);
/// \endcond
