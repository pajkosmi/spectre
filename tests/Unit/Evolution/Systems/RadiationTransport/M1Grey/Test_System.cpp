// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <utility>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/System.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "Evolution/Systems/RadiationTransport/Tags.hpp"
#include "Helpers/Evolution/Imex/TestSector.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

SPECTRE_TEST_CASE("Unit.RadiationTransport.M1Grey.System.Imex",
                  "[Unit][M1Grey]") {
  // struct DummySpecies;
  using DummySpecies = neutrinos::ElectronNeutrinos<1>;

  using system = RadiationTransport::M1Grey::System<tmpl::list<DummySpecies>>;
  using sector = tmpl::front<system::implicit_sectors>;

  // FIXME better values?
  // Metric
  Scalar<DataVector> lapse{};
  get(lapse) = DataVector{1.0};
  tnsr::ii<DataVector, 3> spatial_metric{};
  get<0, 0>(spatial_metric) = DataVector{0.1};
  get<0, 1>(spatial_metric) = DataVector{0.0};
  get<0, 2>(spatial_metric) = DataVector{0.0};
  get<1, 1>(spatial_metric) = DataVector{0.1};
  get<1, 2>(spatial_metric) = DataVector{0.0};
  get<2, 2>(spatial_metric) = DataVector{0.1};
  tnsr::II<DataVector, 3> inverse_spatial_metric{};
  Scalar<DataVector> sqrt_det_spatial_metric{};
  determinant_and_inverse(make_not_null(&sqrt_det_spatial_metric),
                          make_not_null(&inverse_spatial_metric),
                          spatial_metric);
  get(sqrt_det_spatial_metric) = sqrt(get(sqrt_det_spatial_metric));

  // neutrino energy density (fluid frame)
  Scalar<DataVector> tilde_e{};
  get(tilde_e) = DataVector{10.2};
  // neutrino momentum density (fluid frame)
  tnsr::i<DataVector, 3> tilde_s{};
  get<0>(tilde_s) = DataVector{0.1};
  get<1>(tilde_s) = DataVector{0.2};
  get<2>(tilde_s) = DataVector{0.3};
  // neutrino properties
  Scalar<DataVector> emissivity{};
  get(emissivity) = DataVector{2.0};
  Scalar<DataVector> absorption_opacity{};
  get(absorption_opacity) = DataVector{4.0};
  Scalar<DataVector> scattering_opacity{};
  get(scattering_opacity) = DataVector{3.0};
  // fluid velocity
  tnsr::I<DataVector, 3> fluid_velocity{};
  get<0>(fluid_velocity) = DataVector{0.1};
  get<1>(fluid_velocity) = DataVector{0.1};
  get<2>(fluid_velocity) = DataVector{0.1};
  Scalar<DataVector> lorentz_factor{1_st};
  tenex::evaluate(
      make_not_null(&lorentz_factor),
      1.0 / sqrt(1.0 - fluid_velocity(ti::I) * fluid_velocity(ti::J) *
                           spatial_metric(ti::i, ti::j)));
  // recalculate tilde e (energy density)
  tenex::evaluate(make_not_null(&tilde_e),
                  sqrt(tilde_s(ti::i) * tilde_s(ti::j) *
                       inverse_spatial_metric(ti::I, ti::J)) +
                      0.01);  // force nearly thin
  // get(tilde_e) = 4./3. * square(get(lorentz_factor)) - 1./3.; //thick
  // tenex::evaluate<ti::i>(make_not_null(&tilde_s), 4./3. *
  // square(lorentz_factor()) * fluid_velocity(ti::J) * spatial_metric(ti::i,
  // ti::j)); // thick
//   Scalar<DataVector> closure_factor{};
//   get(closure_factor) = DataVector{0.0};
//   Scalar<DataVector> tilde_j{};
//   get(tilde_j) = DataVector{0.0};
//   Scalar<DataVector> tilde_h_normal{};
//   get(tilde_h_normal) = DataVector{0.0};
//   tnsr::i<DataVector, 3> tilde_h_spatial{};
//   get<0>(tilde_h_spatial) = DataVector{1.0};
//   get<1>(tilde_h_spatial) = DataVector{1.0};
//   get<2>(tilde_h_spatial) = DataVector{1.0};
//   tnsr::II<DataVector, 3> tilde_p{};
//   get<0, 0>(tilde_p) = DataVector{0.1};
//   get<0, 1>(tilde_p) = DataVector{0.0};
//   get<0, 2>(tilde_p) = DataVector{0.0};
//   get<1, 1>(tilde_p) = DataVector{0.1};
//   get<1, 2>(tilde_p) = DataVector{0.0};
//   get<2, 2>(tilde_p) = DataVector{0.1};

  const double stencil_size = 1.0e-3;
  const double tolerance = 1.0e-7;

  // has to match tags_from_evolution in System.hpp :: M1Solve

  tuples::TaggedTuple<
      //   RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial,
      //   DummySpecies>,
      //   RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial,
      //   DummySpecies>,
      gr::Tags::Lapse<DataVector>,
      gr::Tags::SpatialMetric<DataVector, 3, Frame::Inertial>,
      gr::Tags::SqrtDetSpatialMetric<DataVector>,
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame::Inertial>,
      RadiationTransport::M1Grey::Tags::GreyEmissivity<DummySpecies>,
      RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<DummySpecies>,
      RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<DummySpecies>,
      hydro::Tags::LorentzFactor<DataVector>,
      hydro::Tags::SpatialVelocity<DataVector, 3>>
      list_of_values = {std::move(lapse),
                        std::move(spatial_metric),
                        std::move(sqrt_det_spatial_metric),
                        std::move(inverse_spatial_metric),
                        std::move(emissivity),
                        std::move(absorption_opacity),
                        std::move(scattering_opacity),
                        std::move(lorentz_factor),
                        std::move(fluid_velocity)};

  using sector_variables_tag = Tags::Variables<sector::tensors>;
  using SectorVariables = sector_variables_tag::type;

  // Mike: I think these are arbitrary?
  // note, here, explicit means not implicit, in the time stepper sense
  SectorVariables explicit_values(1);
  Scalar<DataVector>& tilde_e_vals = get<
      RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial, DummySpecies>>(
      explicit_values);

  get(tilde_e_vals) = DataVector{0.9}; //1.1932159566199232};
  tnsr::i<DataVector, 3>& tilde_s_vals = get<
      RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial, DummySpecies>>(
      explicit_values);
  get<0>(tilde_s_vals) = DataVector{0.01};
  get<1>(tilde_s_vals) = DataVector{0.02};
  get<2>(tilde_s_vals) = DataVector{0.03};

  TestHelpers::imex::test_sector<sector>(stencil_size, tolerance,
                                         explicit_values, list_of_values);
}
