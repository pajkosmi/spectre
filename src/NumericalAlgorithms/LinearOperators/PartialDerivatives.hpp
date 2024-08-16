// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines functions computing partial derivatives.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Mesh;

namespace domain {
namespace Tags {
template <size_t Dim>
struct Mesh;
}  // namespace Tags
}  // namespace domain
namespace Tags {
template <class TagList>
struct Variables;
}  // namespace Tags
/// \endcond

namespace Tags {
/*!
 * \ingroup DataBoxTagsGroup
 * \brief Prefix indicating spatial derivatives
 *
 * Prefix indicating the spatial derivatives of a Tensor.
 *
 * \tparam Tag The tag to wrap
 * \tparam Dim The volume dim as a type (e.g. `tmpl::size_t<Dim>`)
 * \tparam Frame The frame of the derivative index
 *
 * \see Tags::DerivCompute
 */
template <typename Tag, typename Dim, typename Frame, typename = std::nullptr_t>
struct deriv;

template <typename Tag, typename Dim, typename Frame>
struct deriv<Tag, Dim, Frame, Requires<tt::is_a_v<Tensor, typename Tag::type>>>
    : db::PrefixTag, db::SimpleTag {
  using type =
      TensorMetafunctions::prepend_spatial_index<typename Tag::type, Dim::value,
                                                 UpLo::Lo, Frame>;
  using tag = Tag;
};

/*!
 * \ingroup DataBoxTagsGroup
 * \brief Prefix indicating spacetime derivatives
 *
 * Prefix indicating the spacetime derivatives of a Tensor or that a Variables
 * contains spatial derivatives of Tensors.
 *
 * \tparam Tag The tag to wrap
 * \tparam Dim The volume dim as a type (e.g. `tmpl::size_t<Dim>`)
 * \tparam Frame The frame of the derivative index
 */
template <typename Tag, typename Dim, typename Frame, typename = std::nullptr_t>
struct spacetime_deriv;

template <typename Tag, typename Dim, typename Frame>
struct spacetime_deriv<Tag, Dim, Frame,
                       Requires<tt::is_a_v<Tensor, typename Tag::type>>>
    : db::PrefixTag, db::SimpleTag {
  using type =
      TensorMetafunctions::prepend_spacetime_index<typename Tag::type,
                                                   Dim::value, UpLo::Lo, Frame>;
  using tag = Tag;
};

}  // namespace Tags

/// @{
/// \ingroup NumericalAlgorithmsGroup
/// \brief Compute the partial derivatives of each variable with respect to
/// the element logical coordinate.
///
/// \requires `DerivativeTags` to be the head of `VariableTags`
///
/// Returns a `Variables` with a spatial tensor index appended to the front
/// of each tensor within `u` and each `Tag` wrapped with a `Tags::deriv`.
///
/// \tparam DerivativeTags the subset of `VariableTags` for which derivatives
/// are computed.
template <typename DerivativeTags, typename VariableTags, size_t Dim>
void logical_partial_derivatives(
    gsl::not_null<std::array<Variables<DerivativeTags>, Dim>*>
        logical_partial_derivatives_of_u,
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh);

template <typename DerivativeTags, typename VariableTags, size_t Dim>
auto logical_partial_derivatives(const Variables<VariableTags>& u,
                                 const Mesh<Dim>& mesh)
    -> std::array<Variables<DerivativeTags>, Dim>;
/// @}

/// @{
/*!
 * \ingroup NumericalAlgorithmsGroup
 * \brief Computes the logical partial derivative of a tensor, prepending the
 * spatial derivative index, e.g. for \f$\partial_i T_{a}{}^{b}\f$ the C++ call
 * is `get(i, a, b)`.
 *
 * There is an overload that accepts a `buffer` of size
 * `mesh.number_of_grid_points()` or larger. When passed this function performs
 * no memory allocations, which helps improve performance.
 *
 * If you have a `Variables` with several tensors you need to differentiate you
 * should use the `logical_partial_derivatives` function that operates on
 * `Variables` since that'll be more efficient.
 */
template <typename SymmList, typename IndexList, size_t Dim>
void logical_partial_derivative(
    gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>*>
        logical_derivative_of_u,
    gsl::not_null<gsl::span<double>*> buffer,
    const Tensor<DataVector, SymmList, IndexList>& u, const Mesh<Dim>& mesh);

template <typename SymmList, typename IndexList, size_t Dim>
void logical_partial_derivative(
    gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>*>
        logical_derivative_of_u,
    const Tensor<DataVector, SymmList, IndexList>& u, const Mesh<Dim>& mesh);

template <typename SymmList, typename IndexList, size_t Dim>
auto logical_partial_derivative(
    const Tensor<DataVector, SymmList, IndexList>& u, const Mesh<Dim>& mesh)
    -> TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>;
/// @}

/// @{
/// \ingroup NumericalAlgorithmsGroup
/// \brief Compute the partial derivatives of each variable with respect to
/// the coordinates of `DerivativeFrame`.
///
/// Either compute partial derivatives of _all_ variables in `VariableTags`, or
/// of a subset of the `VariablesTags`. The subset of tags (`DerivativeTags`)
/// must be the head of `VariablesTags`.
///
/// The return-by-reference overload infers all template parameters from the
/// arguments. The tensor types in the output buffer must have a spatial index
/// appended to the front.
///
/// The return-by-value overload requires that the `DerivativeTags` are
/// specified explicitly as the first template parameter. It returns a
/// `Variables` with the `DerivativeTags` wrapped in `Tags::deriv`.
template <typename ResultTags, typename DerivativeTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    gsl::not_null<Variables<ResultTags>*> du,
    const std::array<Variables<DerivativeTags>, Dim>&
        logical_partial_derivatives_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian);

template <typename ResultTags, typename DerivativeTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives_cartoon(
    gsl::not_null<Variables<ResultTags>*> du,
    const std::array<Variables<DerivativeTags>, Dim>&
        logical_partial_derivatives_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const gsl::span<const double>& volume_vars,
    const Index<3>& subcell_extents);

template <typename ResultTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    gsl::not_null<Variables<ResultTags>*> du, const Variables<VariableTags>& u,
    const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian);

template <typename DerivativeTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
auto partial_derivatives(
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian)
    -> Variables<db::wrap_tags_in<Tags::deriv, DerivativeTags,
                                  tmpl::size_t<Dim>, DerivativeFrame>>;
/// @}

/// @{
/// \ingroup NumericalAlgorithmsGroup
/// \brief Compute the partial derivative of a `Tensor` with respect to
/// the coordinates of `DerivativeFrame`.
///
/// Returns a `Tensor` with a spatial tensor index appended to the front
/// of the input `Tensor`.
///
/// If you have a `Variables` with several tensors you need to differentiate,
/// you should use the `partial_derivatives` function that operates on
/// `Variables` since that'll be more efficient.
template <typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
void partial_derivative(
    const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        DerivativeFrame>*>
        du,
    const TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>& logical_partial_derivative_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian);

template <typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
void partial_derivative(
    const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        DerivativeFrame>*>
        du,
    const Tensor<DataVector, SymmList, IndexList>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian);

template <typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
auto partial_derivative(
    const Tensor<DataVector, SymmList, IndexList>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian)
    -> TensorMetafunctions::prepend_spatial_index<
        Tensor<DataVector, SymmList, IndexList>, Dim, UpLo::Lo,
        DerivativeFrame>;
/// @}

namespace Tags {

/*!
 * \ingroup DataBoxTagsGroup
 * \brief Compute the spatial derivatives of tags in a Variables
 *
 * Computes the spatial derivatives of the Tensors in the Variables represented
 * by `VariablesTag` in the frame mapped to by `InverseJacobianTag`. To only
 * take the derivatives of a subset of these Tensors you can set the
 * `DerivTags` template parameter. It takes a `tmpl::list` of the desired
 * tags and defaults to the full `tags_list` of the Variables.
 *
 * This tag may be retrieved via `::Tags::Variables<db::wrap_tags_in<deriv,
 * DerivTags, Dim, deriv_frame>`.
 */
template <typename VariablesTag, typename MeshTag, typename InverseJacobianTag,
          typename DerivTags = typename VariablesTag::type::tags_list>
struct DerivCompute
    : db::add_tag_prefix<
          deriv, ::Tags::Variables<DerivTags>,
          tmpl::size_t<
              tmpl::back<typename InverseJacobianTag::type::index_list>::dim>,
          typename tmpl::back<
              typename InverseJacobianTag::type::index_list>::Frame>,
      db::ComputeTag {
 private:
  using inv_jac_indices = typename InverseJacobianTag::type::index_list;
  static constexpr auto Dim = tmpl::back<inv_jac_indices>::dim;
  using deriv_frame = typename tmpl::back<inv_jac_indices>::Frame;

 public:
  using base = db::add_tag_prefix<
      deriv, ::Tags::Variables<DerivTags>,
      tmpl::size_t<
          tmpl::back<typename InverseJacobianTag::type::index_list>::dim>,
      typename tmpl::back<
          typename InverseJacobianTag::type::index_list>::Frame>;
  using return_type = typename base::type;
  static constexpr void (*function)(
      gsl::not_null<return_type*>, const typename VariablesTag::type&,
      const Mesh<Dim>&,
      const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                            deriv_frame>&) =
      partial_derivatives<typename return_type::tags_list,
                          typename VariablesTag::type::tags_list, Dim,
                          deriv_frame>;
  using argument_tags = tmpl::list<VariablesTag, MeshTag, InverseJacobianTag>;
};

/*!
 * \ingroup DataBoxTagsGroup
 * \brief Computes the spatial derivative of a single tensor tag not in a
 * Variables.
 *
 * Computes the spatial derivative for a single tensor represented by
 * 'TensorTag' in the frame mapped to by 'InverseJacobianTag'. It takes a
 * single Tensor designated by 'TensorTag', the inverse Jacobian, and a mesh.
 */
template <typename TensorTag, typename InverseJacobianTag, typename MeshTag>
struct DerivTensorCompute
    : ::Tags::deriv<TensorTag,
                    tmpl::size_t<tmpl::back<
                        typename InverseJacobianTag::type::index_list>::dim>,
                    typename tmpl::back<
                        typename InverseJacobianTag::type::index_list>::Frame>,
      db::ComputeTag {
 private:
  using inv_jac_indices = typename InverseJacobianTag::type::index_list;
  static constexpr auto Dim = tmpl::back<inv_jac_indices>::dim;
  using deriv_frame = typename tmpl::back<inv_jac_indices>::Frame;

 public:
  using base = ::Tags::deriv<TensorTag, tmpl::size_t<Dim>, deriv_frame>;
  using return_type = typename base::type;
  static constexpr void (*function)(
      gsl::not_null<return_type*>, const typename TensorTag::type&,
      const Mesh<Dim>&,
      const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                            deriv_frame>&) =
      partial_derivative<typename TensorTag::type::symmetry,
                         typename TensorTag::type::index_list, Dim,
                         deriv_frame>;
  using argument_tags =
      tmpl::list<TensorTag, domain::Tags::Mesh<Dim>, InverseJacobianTag>;
};
}  // namespace Tags
