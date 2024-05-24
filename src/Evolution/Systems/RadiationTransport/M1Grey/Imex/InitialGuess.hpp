// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <vector>

#include "Evolution/Imex/GuessResult.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"

/// \cond
class DataVector;
template <typename>
class Variables;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

namespace M1Grey {
namespace Imex {

template <typename NeutrinoSpeciesList>
struct InitialGuess;

template <typename... NeutrinoSpecies>
struct InitialGuess<tmpl::list<NeutrinoSpecies...>> {
  using return_tags =
      tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial,
                                                          NeutrinoSpecies>...,
                 RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial,
                                                          NeutrinoSpecies>...>;
  using argument_tags = tmpl::list<
      RadiationTransport::M1Grey::Tags::TildeJ<NeutrinoSpecies>...,
      RadiationTransport::M1Grey::Tags::TildeHSpatial<Frame::Inertial,
                                                      NeutrinoSpecies>...,
      gr::Tags::Lapse<DataVector>, gr::Tags::SpatialMetric<DataVector, 3>>;

  template <typename return_tags>
  static std::vector<imex::GuessResult> apply(
      gsl::not_null<Scalar<DataVector>*> tilde_e,
      gsl::not_null<tnsr::i<DataVector, 3>*> tilde_s,
      const Scalar<DataVector>& tilde_j,
      const tnsr::i<DataVector, 3>& tilde_h_spatial,
      const Scalar<DataVector>& lapse,
      const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
      const Variables<return_tags>& inhomogeneous_terms,
      double implicit_weight) {
    const size_t num_grid_pts = get(lapse).size();

    Variables<tmpl::list<::Tags::TempScalar<0>, ::Tags::TempScalar<1>,
                         ::Tags::TempScalar<2>, ::Tags::TempScalar<3>>>
        buffer{num_grid_pts};

    // MIKE fix this with proper numbers
    (*tilde_e) = tilde_j;

    for (size_t i = 0; i < 3; ++i) {
      (*tilde_s).get(i) = 1.0 * tilde_h_spatial.get(i);
    }

    std::vector<imex::GuessResult> result{num_grid_pts,
                                          imex::GuessResult::ExactSolution};

    for (size_t i = 0; i < num_grid_pts; ++i) {
      result.at(i) = imex::GuessResult::InitialGuess;
    }

    return result;
  }
};

}  // namespace Imex
}  // namespace M1Grey
