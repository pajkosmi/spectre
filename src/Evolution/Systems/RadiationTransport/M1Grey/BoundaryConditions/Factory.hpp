// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/BoundaryConditions/DirichletAnalytic.hpp"
#include "Utilities/TMPL.hpp"

namespace RadiationTransport::M1Grey::BoundaryConditions {
/// Typelist of standard BoundaryConditions
template <size_t Dim, typename NeutrinoSpeciesList>
using standard_boundary_conditions =
    tmpl::list<DirichletAnalytic<Dim, NeutrinoSpeciesList>,
               domain::BoundaryConditions::Periodic<
                   BoundaryCondition<Dim, NeutrinoSpeciesList>>>;
}  // namespace RadiationTransport::M1Grey::BoundaryConditions
