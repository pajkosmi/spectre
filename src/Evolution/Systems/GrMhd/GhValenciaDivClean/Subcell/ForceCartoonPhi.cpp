// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <cstddef>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>

#include "Evolution/Systems/GrMhd/GhValenciaDivClean/Subcell/ForceCartoonPhi.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Initialization/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/FiniteDifference/PartialDerivatives.tpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/TMPL.hpp"

#include <iostream>

namespace grmhd::GhValenciaDivClean::subcell {

template <size_t Dim>
void ForceCartoonPhi<Dim>::apply(
    gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> phi,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coords,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric) {
  ::fd::general_cartoon_deriv(*phi, spacetime_metric, inertial_coords);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) template class ForceCartoonPhi<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
}  // namespace grmhd::GhValenciaDivClean::subcell
