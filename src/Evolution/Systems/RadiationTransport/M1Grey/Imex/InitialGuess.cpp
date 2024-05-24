// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/RadiationTransport/M1Grey/Imex/InitialGuess.hpp"

#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Imex/GuessResult.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

namespace M1Grey::Imex {

static std::vector<imex::GuessResult> apply(
    gsl::not_null<Scalar<DataVector>*> tilde_e,
    gsl::not_null<tnsr::i<DataVector, 3>*> tilde_s,
    const Scalar<DataVector>& tilde_j,
    const tnsr::i<DataVector, 3>& tilde_h_spatial,
    const Scalar<DataVector>& lapse,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Variables<tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<
                                   Frame::Inertial, NeutrinoSpecies>...,
                               RadiationTransport::M1Grey::Tags::TildeS<
                                   Frame::Inertial, NeutrinoSpecies>...>>&
        inhomogeneous_terms,
    const double implicit_weight) {
  const size_t num_grid_pts = get(lapse).size();

  Variables<tmpl::list<::Tags::TempScalar<0>, ::Tags::TempScalar<1>,
                       ::Tags::TempScalar<2>, ::Tags::TempScalar<3>>>
      buffer{num_grid_pts};

  // MIKE fix this with proper numbers
  (*tilde_e) = 1.0 * tilde_j;

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

}  // namespace M1Grey::Imex
