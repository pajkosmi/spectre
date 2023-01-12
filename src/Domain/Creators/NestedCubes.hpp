// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <memory>
#include <vector>

#include "DataStructures/Index.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/BoundaryConditions/GetBoundaryConditionsBase.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Options/Options.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace CoordinateMaps {
class Affine;
template <typename Map1, typename Map2>
class ProductOf2Maps;
template <typename Map1, typename Map2, typename Map3>
class ProductOf3Maps;
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain
/// \endcond

namespace domain {
namespace creators {

template <size_t VolumeDim>
class NestedCubes : public DomainCreator<VolumeDim> {
 public:
  using maps_list = tmpl::list<
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                            CoordinateMaps::Affine>,
      domain::CoordinateMap<
          Frame::BlockLogical, Frame::Inertial,
          CoordinateMaps::ProductOf2Maps<CoordinateMaps::Affine,
                                         CoordinateMaps::Affine>>,
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                            CoordinateMaps::ProductOf3Maps<
                                CoordinateMaps::Affine, CoordinateMaps::Affine,
                                CoordinateMaps::Affine>>>;

  struct DomainBounds {
    using type = std::array<double, 2>;
    static constexpr Options::String help = {
        "Min and amx coordinates of domain outer boundary. The min and max "
        "will be applied to each direction."};
  };

  struct BlocksPerDim {
    using type = size_t;
    static constexpr Options::String help = {
        "Number of nested cubes. A single cube (block) with has zero nested "
        "levels. A Rubix cube (3 blocks in each dimension where one of those "
        "blocks is the central cube) has 1 nested level."};
  };

  struct IsPeriodicIn {
    using type = std::array<bool, VolumeDim>;
    static constexpr Options::String help = {
        "Whether the domain is periodic in each dimension."};
  };

  struct MaxRefinementLevel {
    using type = size_t;
    static constexpr Options::String help = {
        "Max refinement level in the inner-most block."};
  };

  struct NumberOfGridPoints {
    using type = size_t;
    static constexpr Options::String help = {
        "Initial number of grid points at each nested level."};
  };

  template <typename BoundaryConditionsBase>
  struct BoundaryCondition {
    static std::string name() { return "BoundaryCondition"; }
    static constexpr Options::String help =
        "The boundary condition to impose on all sides.";
    using type = std::unique_ptr<BoundaryConditionsBase>;
  };

  using common_options = tmpl::list<DomainBounds, MaxRefinementLevel,
                                    BlocksPerDim, NumberOfGridPoints>;
  using options_periodic = tmpl::list<IsPeriodicIn>;

  template <typename Metavariables>
  using options = tmpl::append<
      common_options,
      tmpl::conditional_t<
          domain::BoundaryConditions::has_boundary_conditions_base_v<
              typename Metavariables::system>,
          tmpl::list<BoundaryCondition<
              domain::BoundaryConditions::get_boundary_conditions_base<
                  typename Metavariables::system>>>,
          options_periodic>>;

  static constexpr Options::String help = {"Nested cubes."};

  NestedCubes(typename DomainBounds::type domain_bounds,
              typename MaxRefinementLevel::type max_refinement_level,
              typename BlocksPerDim::type blocks_per_dim,
              typename NumberOfGridPoints::type grid_points,
              typename IsPeriodicIn::type is_periodic_in,
              const Options::Context& context = {});

  NestedCubes(typename DomainBounds::type domain_bounds,
              typename MaxRefinementLevel::type max_refinement_level,
              typename BlocksPerDim::type blocks_per_dim,
              typename NumberOfGridPoints::type grid_points,
              std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
                  boundary_condition = nullptr,
              const Options::Context& context = {});

  NestedCubes() = default;
  NestedCubes(const NestedCubes&) = delete;
  NestedCubes(NestedCubes&&) = default;
  NestedCubes& operator=(const NestedCubes&) = delete;
  NestedCubes& operator=(NestedCubes&&) = default;
  ~NestedCubes() override = default;

  Domain<VolumeDim> create_domain() const override;

  std::vector<DirectionMap<
      VolumeDim,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
  external_boundary_conditions() const override;

  std::vector<std::array<size_t, VolumeDim>> initial_extents() const override;

  std::vector<std::array<size_t, VolumeDim>> initial_refinement_levels()
      const override;

 private:
  void common_construction(const typename DomainBounds::type& domain_bounds,
                           const typename Options::Context& context);

  typename std::array<std::vector<double>, VolumeDim> block_bounds_{};
  typename IsPeriodicIn::type is_periodic_in_{make_array<VolumeDim>(false)};
  typename MaxRefinementLevel::type max_refinement_{};
  typename BlocksPerDim::type blocks_per_dim_{};
  typename NumberOfGridPoints::type grid_points_{};
  domain::creators::AlignedLattice<VolumeDim> aligned_lattice_;
  Index<VolumeDim> number_of_blocks_by_dim_{};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      boundary_condition_{};
};
}  // namespace creators
}  // namespace domain
