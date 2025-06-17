// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/DgSubcell/Tags/Coordinates.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/PrimitiveFromConservativeOptions.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
/// \endcond

namespace grmhd::GhValenciaDivClean::subcell {
/*!
 * \brief Mutator that forces cartoon Phi.
 *
 * In the DG-subcell `step_actions` list this will normally be called using the
 * `::Actions::MutateApply` action in the following way in the action list:
 * - `Actions::Label<Labels::BeginSubcellAfterDgRollback>`
 * - `Actions::MutateApply<ForceCartoonPhi>`
 */
template <size_t Dim>
struct ForceCartoonPhi {
  using return_tags = tmpl::list<gh::Tags::Phi<DataVector, Dim>>;
  using argument_tags = tmpl::list<
      evolution::dg::subcell::Tags::Coordinates<Dim, Frame::Inertial>,
      gr::Tags::SpacetimeMetric<DataVector, Dim>>;

  static void apply(
      gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> phi,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coordinates,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric);
};
}  // namespace grmhd::GhValenciaDivClean::subcell
