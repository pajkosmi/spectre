// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <memory>

#include <iostream>
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Informer/InfoFromBuild.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/HomogeneousSphere.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {
namespace {

void test_equality(double sphere_radius, double emissivity_and_opacity,
                   double outer_radius, double outer_opacity) {
  register_classes_with_charm<
      RadiationTransport::M1Grey::AnalyticData::HomogeneousSphereImpl>();

  // Inner radius should be smaller than outer radius
  CHECK_THROWS_WITH(HomogeneousSphereImpl(sphere_radius, emissivity_and_opacity,
                                          sphere_radius / 2.0, outer_opacity),
                    Catch::Matchers::ContainsSubstring("is greater than the"));

  const HomogeneousSphereImpl homogeneous_original{
      sphere_radius, emissivity_and_opacity, outer_radius, outer_opacity};

  const auto homogeneous_sphere =
      serialize_and_deserialize(homogeneous_original);

  // serialize/deserialize should still be the same
  CHECK(homogeneous_sphere ==
        HomogeneousSphereImpl(sphere_radius, emissivity_and_opacity,
                              outer_radius, outer_opacity));

  CHECK(homogeneous_sphere !=
        HomogeneousSphereImpl(sphere_radius + 0.1, emissivity_and_opacity,
                              outer_radius, outer_opacity));
  CHECK(homogeneous_sphere !=
        HomogeneousSphereImpl(sphere_radius, emissivity_and_opacity + 0.1,
                              outer_radius, outer_opacity));
  CHECK(homogeneous_sphere !=
        HomogeneousSphereImpl(sphere_radius, emissivity_and_opacity,
                              outer_radius + 0.1, outer_opacity));

  CHECK(homogeneous_sphere !=
        HomogeneousSphereImpl(sphere_radius, emissivity_and_opacity,
                              outer_radius, outer_opacity + 0.1));
}

void test_consistency(double sphere_radius, double emissivity_and_opacity,
                      double outer_radius, double outer_opacity) {
  register_classes_with_charm<
      RadiationTransport::M1Grey::AnalyticData::HomogeneousSphereImpl>();

  const std::unique_ptr<evolution::initial_data::InitialData> option_solution =
      TestHelpers::test_option_tag_factory_creation<
          evolution::initial_data::OptionTags::InitialData,
          RadiationTransport::M1Grey::AnalyticData::HomogeneousSphereImpl>(
          "HomogeneousSphereImpl:\n"
          "  Radius: " +
          std::to_string(sphere_radius) +
          "\n"
          "  EmissivityAndOpacity: " +
          std::to_string(emissivity_and_opacity) +
          "\n"
          "  OuterRadius: " +
          std::to_string(outer_radius) +
          "\n"
          "  OuterOpacity: " +
          std::to_string(outer_opacity) + "\n")
          ->get_clone();

  const auto deserialized_option_solution =
      serialize_and_deserialize(option_solution);

  const auto& homogeneous_sphere = dynamic_cast<
      const RadiationTransport::M1Grey::AnalyticData::HomogeneousSphereImpl&>(
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

  const auto vars = homogeneous_sphere.variables(
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
    "HomogeneousSphere",
    "[Unit][PointwiseFunctions]") {
  const double sphere_radius = 1.0;
  const double emissivity_and_opacity = 1.0;
  const double outer_radius = 1.5;
  const double outer_opacity = 0.5;

  test_equality(sphere_radius, emissivity_and_opacity, outer_radius,
                outer_opacity);

  test_consistency(sphere_radius, emissivity_and_opacity, outer_radius,
                   outer_opacity);
}

}  // namespace
}  // namespace RadiationTransport::M1Grey::AnalyticData
