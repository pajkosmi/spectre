// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <memory>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Informer/InfoFromBuild.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/CcsnCollapse.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Utilities/GetOutput.hpp"

namespace grmhd::AnalyticData {
namespace {

static_assert(
    not is_analytic_solution_v<CcsnCollapse>,
    "CcsnCollapse should be analytic_data and not an analytic_solution");
static_assert(
    is_analytic_data_v<CcsnCollapse>,
    "CcsnCollapse should be analytic_data and not an analytic_solution");

// Make sure different parameters give different answers
// void test_equality(std::string progenitor_filename, double
// polytropic_constant,
//                    double adiabatic_index, double central_angular_velocity,
//                    double diff_rot_parameter, double max_dens_ratio)

void test_equality(std::string progenitor_filename,
                   double polytropic_transition_density,
                   double polytropic_constant_low, double adiabatic_index_low,
                   double adiabatic_index_high, double thermal_adiabatic_index,
                   double central_angular_velocity, double diff_rot_parameter,
                   double max_dens_ratio_interp) {
  Parallel::register_classes_with_charm<grmhd::AnalyticData::CcsnCollapse>();
  // Base case for comparison
  //   const CcsnCollapse ccsn_progenitor_original{
  //       progenitor_filename,      polytropic_constant, adiabatic_index,
  //       central_angular_velocity, diff_rot_parameter,  max_dens_ratio};

  const CcsnCollapse ccsn_progenitor_original{
      progenitor_filename,      polytropic_transition_density,
      polytropic_constant_low,  adiabatic_index_low,
      adiabatic_index_high,     thermal_adiabatic_index,
      central_angular_velocity, diff_rot_parameter,
      max_dens_ratio_interp};

  //   CcsnCollapse(std::string progenitor_filename,
  //                double polytropic_transition_density,
  //                double polytropic_constant_low, double adiabatic_index_low,
  //                double adiabatic_index_high, double thermal_adiabatic_index,
  //                double central_angular_velocity, double diff_rot_parameter,
  //                double max_dens_ratio_interp)

  CHECK_THROWS_WITH(
      CcsnCollapse(progenitor_filename + "badname",
                   polytropic_transition_density, polytropic_constant_low,
                   adiabatic_index_low, adiabatic_index_high,
                   thermal_adiabatic_index, central_angular_velocity,
                   diff_rot_parameter, max_dens_ratio_interp),
      Catch::Matchers::Contains("Data file not found"));

  const auto ccsn_progenitor =
      serialize_and_deserialize(ccsn_progenitor_original);
  CHECK(ccsn_progenitor ==
        CcsnCollapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low, adiabatic_index_low,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity, diff_rot_parameter,
                     max_dens_ratio_interp));
  CHECK(ccsn_progenitor !=
        CcsnCollapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low + 0.1, adiabatic_index_low,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity, diff_rot_parameter,
                     max_dens_ratio_interp));
  CHECK(ccsn_progenitor !=
        CcsnCollapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low + 0.1, adiabatic_index_low + 0.1,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity, diff_rot_parameter,
                     max_dens_ratio_interp));
  CHECK(ccsn_progenitor !=
        CcsnCollapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low + 0.1, adiabatic_index_low,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity + 0.1, diff_rot_parameter,
                     max_dens_ratio_interp));
  CHECK(ccsn_progenitor !=
        CcsnCollapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low + 0.1, adiabatic_index_low,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity, diff_rot_parameter + 0.1,
                     max_dens_ratio_interp));
  // Add more differences?
}

void test_ccsn_collapse(std::string progenitor_filename,
                        double polytropic_transition_density,
                        double polytropic_constant_low,
                        double adiabatic_index_low, double adiabatic_index_high,
                        double thermal_adiabatic_index,
                        double central_angular_velocity,
                        double diff_rot_parameter,
                        double max_dens_ratio_interp) {
  Parallel::register_classes_with_charm<grmhd::AnalyticData::CcsnCollapse>();

  const std::unique_ptr<evolution::initial_data::InitialData> option_solution =
      TestHelpers::test_option_tag_factory_creation<
          evolution::initial_data::OptionTags::InitialData,
          grmhd::AnalyticData::CcsnCollapse>(
          "CcsnCollapse:\n"
          "  ProgenitorFilename: " +
          std::to_string(progenitor_filename) +
          "\n"
          "  PiecewisePolytropicTransitionDensity: " +
          std::to_string(polytropic_transition_densit) +
          "\n"
          "  ThermalAdiabaticIndex: " +
          std::to_string(thermal_adiabatic_index) +
          "\n"
          "  PolytropicAdiabaticIndexLow: " +
          std::to_string(adiabatic_index_low) +
          "\n"
          "  PolytropicAdiabaticIndexHigh: " +
          std::to_string(adiabatic_index_high) +
          "\n"
          "  CentralAngularVelocity: " +
          std::to_string(central_angular_velocity) +
          "\n"
          "  DifferentialRotationParameter: " +
          std::to_string(diff_rot_parameter) +
          "\n"
          "  PolytropicConstantLow: " +
          std::to_string(polytropic_constant_low) +
          "\n"
          "  MaxDensityRatioForLinearInterpolation: " +
          std::to_string(max_dens_ratio_interp) + "\n")
          ->get_clone();

  const auto deserialized_option_solution =
      serialize_and_deserialize(option_solution);
  const auto& ccsn_progenitor =
      dynamic_cast<const grmhd::AnalyticData::CcsnCollapse&>(
          *deserialized_option_solution);

  const Mesh<3> mesh{
      {{5, 5, 5}},
      {{Spectral::Basis::Legendre, Spectral::Basis::Legendre,
        Spectral::Basis::Legendre}},
      {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::GaussLobatto,
        Spectral::Quadrature::GaussLobatto}}};
  const auto log_coords = logical_coordinates(mesh);

  // Coordinates where we check the data. Includes the origin.
  tnsr::I<DataVector, 3, Frame::Inertial> in_coords{
      mesh.number_of_grid_points(), 0.0};
  for (size_t i = 0; i < 3; ++i) {
    in_coords.get(i) = 1.0e-2 * (log_coords.get(i) + 1.);
  }

  const size_t num_radial_points = 1;
  tnsr::I<DataVector, 3, Frame::Inertial> in_coords_large_radius{
      num_radial_points, 1.0e30};

  INFO("Check Physicality of Results");

  const auto vars = ccsn_progenitor.variables(
      in_coords,
      tmpl::list<
          hydro::Tags::RestMassDensity<DataVector>,
          hydro::Tags::ElectronFraction<DataVector>,
          hydro::Tags::SpatialVelocity<DataVector, 3>,
          hydro::Tags::SpecificEnthalpy<DataVector>,
          hydro::Tags::Pressure<DataVector>,
          hydro::Tags::SpecificInternalEnergy<DataVector>,
          hydro::Tags::LorentzFactor<DataVector>,
          hydro::Tags::MagneticField<DataVector, 3>,
          hydro::Tags::DivergenceCleaningField<DataVector>,
          gr::Tags::Lapse<DataVector>,
          gr::Tags::Shift<3, Frame::Inertial, DataVector>,
          gr::Tags::SpatialMetric<3, Frame::Inertial, DataVector>,
          gr::Tags::SqrtDetSpatialMetric<DataVector>,
          gr::Tags::InverseSpatialMetric<3, Frame::Inertial, DataVector>,
          gr::Tags::ExtrinsicCurvature<3, Frame::Inertial, DataVector>,
          ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                        Frame::Inertial>,
          ::Tags::deriv<gr::Tags::Shift<3, Frame::Inertial, DataVector>,
                        tmpl::size_t<3>, Frame::Inertial>,
          ::Tags::deriv<gr::Tags::SpatialMetric<3, Frame::Inertial, DataVector>,
                        tmpl::size_t<3>, Frame::Inertial>>{});

  // Check velocity > c check and provide small interpolation ratio test
  // coverage
#ifdef SPECTRE_DEBUG
  const std::unique_ptr<evolution::initial_data::InitialData>
      v_grtr_than_c_solution =
          TestHelpers::test_option_tag_factory_creation<
              evolution::initial_data::OptionTags::InitialData,
              grmhd::AnalyticData::CcsnCollapse>(
              "CcsnCollapse:\n"
              "  ProgenitorFilename: " +
              std::to_string(progenitor_filename) +
              "\n"
              "  PiecewisePolytropicTransitionDensity: " +
              std::to_string(polytropic_transition_density) +
              "\n"
              "  ThermalAdiabaticIndex: " +
              std::to_string(thermal_adiabatic_index) +
              "\n"
              "  PolytropicAdiabaticIndexLow: " +
              std::to_string(adiabatic_index_low) +
              "\n"
              "  PolytropicAdiabaticIndexHigh: " +
              std::to_string(adiabatic_index_high) +
              "\n"
              "  CentralAngularVelocity: 147669\n"
              "  DifferentialRotationParameter: 1.0e20\n"
              "  PolytropicConstantLow: " +
              std::to_string(polytropic_constant_low) +
              "\n"
              "  MaxDensityRatioForLinearInterpolation: 0.0\n")
              ->get_clone();
  //   "CcsnCollapse:\n"
  //   "  ProgenitorFilename: " +
  //   progenitor_filename +
  //   "\n"
  //   "  AdiabaticIndex: " +
  //   std::to_string(adiabatic_index) +
  //   "\n"
  //   "  CentralAngularVelocity: 147669\n"
  //   "  DifferentialRotationParameter: 1.0e20\n"
  //   "  PolytropicConstant: " +
  //   std::to_string(polytropic_constant) +
  //   "\n"
  //   "  MaxDensityRatioForLinearInterpolation: 0.0\n")
  //   ->get_clone();

  const auto deserialized_v_grtr_than_c_solution =
      serialize_and_deserialize(v_grtr_than_c_solution);
  const auto& ccsn_progenitor_v_grtr_than_c =
      dynamic_cast<const grmhd::AnalyticData::CcsnCollapse&>(
          *deserialized_v_grtr_than_c_solution);

  CHECK_THROWS_WITH(
      ccsn_progenitor_v_grtr_than_c.variables(
          in_coords, tmpl::list<hydro::Tags::SpatialVelocity<DataVector, 3>>{}),
      Catch::Matchers::Contains("Spatial velocity"));
#endif

  // Check radius too large check
  CHECK_THROWS_WITH(ccsn_progenitor.variables(
                        in_coords_large_radius,
                        tmpl::list<hydro::Tags::RestMassDensity<DataVector>>{}),
                    Catch::Matchers::Contains("Requested radius "));

  // Ensure density is positive
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataVector>>(vars);

  CHECK(min(get(rest_mass_density)) > 0.0);

  // Relevant metric variables
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<3, Frame::Inertial, DataVector>>(vars);
  const auto& sqrt_det_spatial_metric_analytic =
      get<gr::Tags::SqrtDetSpatialMetric<DataVector>>(vars);
  const auto& det_spatial_metric_numeric =
      determinant_and_inverse(spatial_metric).first;

  // Parabolic interpolation is used for metric variables,
  // possibly causing slight deviations below 1.0 (flat metric).
  Approx custom_approx = Approx::custom().epsilon(1.e-11).scale(1.0);

  // Metric must always be greater or equal to 1.0 (i.e. flat space)
  CHECK(min(get(sqrt_det_spatial_metric_analytic)) == custom_approx(1.0));

  // Same for numerically calculated determinant case
  CHECK(min(sqrt(get(det_spatial_metric_numeric))) == custom_approx(1.0));

  // Analytic and numeric determinant of the spatial metric
  // should be equal (to machine precision)
  CHECK(max(abs(get(sqrt_det_spatial_metric_analytic) -
                sqrt(get(det_spatial_metric_numeric)))) == approx(0.0));
}

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticData.GrMhd.CcsnCollapse",
                  "[Unit][PointwiseFunctions]") {
  // Load proper data at below file path
  const std::string progenitor_filename =
      unit_test_src_path() +
      "PointwiseFunctions/AnalyticData/GrMhd/"
      "CcsnCollapseID.dat";
  //   const double polytropic_constant = 0.11943;
  //   const double adiabatic_index = 1.3;
  //   const double central_angular_velocity = 1.0e-5;
  //   const double diff_rot_parameter = 339.0;
  //   const double max_dens_ratio = 100.0;

  const double polytropic_transition_density = 1.6e-4;
  const double polytropic_constant_low = 0.11943;
  const double adiabatic_index_low = 1.3;
  const double adiabatic_index_high = 1.666666667;
  const double thermal_adiabatic_index = 1.7;
  const double central_angular_velocity = 1.0e-5;
  const double diff_rot_parameter = 339.0;
  const double max_dens_ratio_interp = 100.0;

  // Test if (de)serialized data are equal
  test_equality(progenitor_filename, polytropic_transition_density,
                polytropic_constant_low, adiabatic_index_low,
                adiabatic_index_high, thermal_adiabatic_index,
                central_angular_velocity, diff_rot_parameter,
                max_dens_ratio_interp);

  // Check physicality of interpolated data
  //   test_ccsn_collapse(progenitor_filename, polytropic_constant,
  //   adiabatic_index,
  //                      central_angular_velocity, diff_rot_parameter,
  //                      max_dens_ratio);
  test_ccsn_collapse(progenitor_filename, polytropic_transition_density,
                     polytropic_constant_low, adiabatic_index_low,
                     adiabatic_index_high, thermal_adiabatic_index,
                     central_angular_velocity, diff_rot_parameter,
                     max_dens_ratio_interp);
}

}  // namespace
}  // namespace grmhd::AnalyticData
