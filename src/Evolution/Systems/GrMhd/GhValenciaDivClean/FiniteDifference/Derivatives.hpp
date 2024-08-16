// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/GrMhd/GhValenciaDivClean/System.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Mesh;
template <typename TagsList>
class Variables;
namespace evolution::dg::subcell {
class GhostData;
}  // namespace evolution::dg::subcell
/// \endcond

namespace grmhd::GhValenciaDivClean::fd {
/*!
 * \brief Compute partial derivatives of the spacetime variables \f$g_{ab}\f$,
 * \f$\Phi_{iab}\f$, and \f$\Pi_{ab}\f$.
 *
 * The derivatives are computed using FD of order deriv_order
 */
void spacetime_derivatives(
    gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::deriv,
        typename grmhd::GhValenciaDivClean::System::gradients_tags,
        tmpl::size_t<3>, Frame::Inertial>>*>
        result,
    const Variables<
        typename grmhd::GhValenciaDivClean::System::variables_tag::tags_list>&
        volume_evolved_variables,
    const DirectionalIdMap<3, evolution::dg::subcell::GhostData>&
        all_ghost_data,
    const size_t& deriv_order, const Mesh<3>& volume_mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>&
        cell_centered_logical_to_inertial_inv_jacobian,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords);
}  // namespace grmhd::GhValenciaDivClean::fd
