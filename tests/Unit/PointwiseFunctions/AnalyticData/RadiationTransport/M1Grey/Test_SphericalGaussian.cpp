// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <memory>

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Informer/InfoFromBuild.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/RotationallySymmetric.tpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {
namespace {

// Tests, equality code from homogeneous sphere
template <size_t Dim>
void test_equality(double gaussian_width,
                   double emissivity_and_opacity_amplitude, double outer_radius,
                   double outer_opacity) {
  register_classes_with_charm<
      RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<Dim>>();

  // Inner radius should be smaller than outer radius
  CHECK_THROWS_WITH(
      SphericalGaussian<Dim>(gaussian_width, emissivity_and_opacity_amplitude,
                             gaussian_width / 2.0, outer_opacity),
      Catch::Matchers::ContainsSubstring("is greater than the"));

  const SphericalGaussian<Dim> spherical_gaussian_original{
      gaussian_width, emissivity_and_opacity_amplitude, outer_radius,
      outer_opacity};

  const auto spherical_gaussian =
      serialize_and_deserialize(spherical_gaussian_original);

  // serialize/deserialize should still be the same
  CHECK(spherical_gaussian ==
        SphericalGaussian<Dim>(gaussian_width, emissivity_and_opacity_amplitude,
                               outer_radius, outer_opacity));

  CHECK(spherical_gaussian !=
        SphericalGaussian<Dim>(gaussian_width + 0.2,
                               emissivity_and_opacity_amplitude, outer_radius,
                               outer_opacity));

  CHECK(spherical_gaussian !=
        SphericalGaussian<Dim>(gaussian_width,
                               emissivity_and_opacity_amplitude + 0.2,
                               outer_radius, outer_opacity));

  CHECK(spherical_gaussian !=
        SphericalGaussian<Dim>(gaussian_width, emissivity_and_opacity_amplitude,
                               outer_radius + 0.2, outer_opacity));

  CHECK(spherical_gaussian !=
        SphericalGaussian<Dim>(gaussian_width, emissivity_and_opacity_amplitude,
                               outer_radius, outer_opacity + 0.2));
}

template <size_t Dim>
void test_consistency(double gaussian_width,
                      double emissivity_and_opacity_amplitude,
                      double outer_radius, double outer_opacity) {
  register_classes_with_charm<
      RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<Dim>>();

  const std::unique_ptr<evolution::initial_data::InitialData> option_solution =
      TestHelpers::test_option_tag_factory_creation<
          evolution::initial_data::OptionTags::InitialData,
          RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<Dim>>(
          RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<
              Dim>::name() +
          ":\n"
          "  Radius: " +
          std::to_string(gaussian_width) +
          "\n"
          "  EmissivityAndOpacityAmplitude: " +
          std::to_string(emissivity_and_opacity_amplitude) +
          "\n"
          "  OuterRadius: " +
          std::to_string(outer_radius) +
          "\n"
          "  OuterOpacity: " +
          std::to_string(outer_opacity) + "\n")
          ->get_clone();

  const auto deserialized_option_solution =
      serialize_and_deserialize(option_solution);

  const auto& spherical_gaussian = dynamic_cast<
      const RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<Dim>&>(
      *deserialized_option_solution);

  // establish sample grid
  const Mesh<3> mesh{
      {{7, 7, 7}},
      {{Spectral::Basis::Legendre, Spectral::Basis::Legendre,
        Spectral::Basis::Legendre}},
      {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::GaussLobatto,
        Spectral::Quadrature::GaussLobatto}}};
  const auto log_coords = logical_coordinates(mesh);

  // Coordinates where we check the data. Includes inside and outside sphere
  tnsr::I<DataVector, 3> in_coords{mesh.number_of_grid_points(), 0.0};
  for (size_t i = 0; i < 3; ++i) {
    in_coords.get(i) = (log_coords.get(i) + 1.);
  }

  using NeutrinoSpecies = neutrinos::ElectronNeutrinos<1>;

  const auto vars = spherical_gaussian.variables(
      in_coords,
      tmpl::list<
          RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial,
                                                   NeutrinoSpecies>,
          RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial,
                                                   NeutrinoSpecies>,
          RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>,
          RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<
              NeutrinoSpecies>,
          RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<
              NeutrinoSpecies>,
          hydro::Tags::LorentzFactor<DataVector>,
          hydro::Tags::SpatialVelocity<DataVector, 3>>{});

  // Access the variables
  const auto tilde_e =
      get<RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial,
                                                   NeutrinoSpecies>>(vars);

  const auto tilde_s =
      get<RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial,
                                                   NeutrinoSpecies>>(vars);

  const auto grey_emissivity =
      get<RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>>(
          vars);

  const auto grey_absorption_opacity = get<
      RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<NeutrinoSpecies>>(
      vars);

  const auto grey_scattering_opacity = get<
      RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<NeutrinoSpecies>>(
      vars);

  const auto lorentz_factor = get<hydro::Tags::LorentzFactor<DataVector>>(vars);

  const auto velocity = get<hydro::Tags::SpatialVelocity<DataVector, 3>>(vars);

  const auto radius =
      RadiationTransport::M1Grey::AnalyticData::radius<Dim == 3>(in_coords);

  // Dimensionally dependent checks
  if (Dim == 2) {
    CHECK(RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<
              Dim>::name() == "CylindricalGaussian");
    CHECK(radius ==
          sqrt(square(get<0>(in_coords)) + square(get<1>(in_coords))));
  } else {
    CHECK(RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<
              Dim>::name() == "SphericalGaussian");
    CHECK(radius == sqrt(square(get<0>(in_coords)) + square(get<1>(in_coords)) +
                         square(get<2>(in_coords))));
  }

  // Spherical Gaussian functions TODO: these are private
  //   CHECK(RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<
  //               Dim>::energy_profile(in_coords));

  // Check points are not uniform for tildeE
  CHECK(min(get(tilde_e)) < max(get(tilde_e)));

  // All components of momentum density should be zero for this problem
  CHECK(min(get<0>(tilde_s)) == 0.0);
  CHECK(min(get<0>(tilde_s)) == max(get<0>(tilde_s)));
  CHECK(min(get<1>(tilde_s)) == 0.0);
  CHECK(min(get<1>(tilde_s)) == max(get<1>(tilde_s)));
  CHECK(min(get<2>(tilde_s)) == 0.0);
  CHECK(min(get<2>(tilde_s)) == max(get<2>(tilde_s)));

  // Grey emissivity in sphere should be greater than outer region
  CHECK(min(get(grey_emissivity)) < max(get(grey_emissivity)));

  // Absorption opacity in sphere should be greater than outer region
  CHECK(min(get(grey_absorption_opacity)) < max(get(grey_absorption_opacity)));

  // Scattering opacity should be zero
  CHECK(min(get(grey_scattering_opacity)) == 0.0);
  CHECK(min(get(grey_scattering_opacity)) == max(get(grey_scattering_opacity)));

  // Lorentz factor should be 1
  CHECK(min(get(lorentz_factor)) == 1.0);
  CHECK(min(get(lorentz_factor)) == max(get(lorentz_factor)));

  // Velocity components should be zero
  CHECK(min(get<0>(velocity)) == 0.0);
  CHECK(min(get<0>(velocity)) == max(get<0>(velocity)));
  CHECK(min(get<1>(velocity)) == 0.0);
  CHECK(min(get<1>(velocity)) == max(get<1>(velocity)));
  CHECK(min(get<2>(velocity)) == 0.0);
  CHECK(min(get<2>(velocity)) == max(get<2>(velocity)));
}

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.RadiationTransport.M1Grey."
    "SphericalGaussian",
    "[Unit][PointwiseFunctions]") {
  const double gaussian_width = 1.0;
  const double emissivity_and_opacity_amplitude = 1.0;
  const double outer_radius = 2.0;
  const double outer_opacity = 0.5;

  test_equality<2>(gaussian_width, emissivity_and_opacity_amplitude,
                   outer_radius, outer_opacity);
  // test_equality<3>(gaussian_width, emissivity_and_opacity_amplitude,
  //  outer_radius, outer_opacity);
  test_consistency<2>(gaussian_width, emissivity_and_opacity_amplitude,
                      outer_radius, outer_opacity);
}

}  // namespace
}  // namespace RadiationTransport::M1Grey::AnalyticData
