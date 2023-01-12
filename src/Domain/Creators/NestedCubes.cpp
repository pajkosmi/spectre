// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/NestedCubes.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "Domain/BoundaryConditions/None.hpp"
#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Domain/Domain.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdArrayHelpers.hpp"

#include "Parallel/Printf.hpp"

namespace domain::creators {
template <size_t VolumeDim>
void NestedCubes<VolumeDim>::common_construction(
    const typename DomainBounds::type& domain_bounds,
    const Options::Context& context) {
  if (domain_bounds[1] <= domain_bounds[0]) {
    PARSE_ERROR(context, "Max domain bound "
                             << domain_bounds[1]
                             << " must be greater than min domain bound "
                             << domain_bounds[0]);
  }

  Parallel::printf("Max refinement %d\n", max_refinement_);
  Parallel::printf("Blocks per dim %d\n", blocks_per_dim_);
  Parallel::printf("Grid points %d\n", grid_points_);
  Parallel::printf("Number of blocks by dim %s\n", number_of_blocks_by_dim_);
  // Parallel::printf("Domain bounds %s\n", block_bounds_);

  if (blocks_per_dim_ % 2 == 0) {
    PARSE_ERROR(context, "Blocks per dim must be odd.");
  }

  const double domain_width = domain_bounds[1] - domain_bounds[0];
  for (size_t i = 0; i < VolumeDim; i++) {
    std::vector<double> bounds(blocks_per_dim_ + 1);
    for (size_t j = 0; j < blocks_per_dim_ + 1; j++) {
      bounds[j] =
          domain_bounds[0] + double(j) / double(blocks_per_dim_) * domain_width;
    }
    // This should only be off by roundoff, but we set it exactly anyways
    bounds[blocks_per_dim_] = domain_bounds[1];

    gsl::at(block_bounds_, i) = std::move(bounds);
  }
}

template <size_t VolumeDim>
NestedCubes<VolumeDim>::NestedCubes(
    typename DomainBounds::type domain_bounds,
    typename MaxRefinementLevel::type max_refinement,
    typename BlocksPerDim::type blocks_per_dim,
    typename NumberOfGridPoints::type grid_points,
    typename IsPeriodicIn::type is_periodic_in, const Options::Context& context)
    // clang-tidy: trivially copyable
    : is_periodic_in_(std::move(is_periodic_in)),  // NOLINT
      max_refinement_(                             // NOLINT
          std::move(max_refinement)),              // NOLINT
      blocks_per_dim_(                             // NOLINT
          std::move(blocks_per_dim)),              // NOLINT
      grid_points_(std::move(grid_points)),
      number_of_blocks_by_dim_(blocks_per_dim_),
      boundary_condition_(nullptr) {
  common_construction(domain_bounds, context);
}

template <size_t VolumeDim>
NestedCubes<VolumeDim>::NestedCubes(
    typename DomainBounds::type domain_bounds,
    typename MaxRefinementLevel::type max_refinement,
    typename BlocksPerDim::type blocks_per_dim,
    typename NumberOfGridPoints::type grid_points,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        boundary_condition,
    const Options::Context& context)
    : is_periodic_in_(make_array<VolumeDim>(false)),  // NOLINT
      max_refinement_(                                // NOLINT
          std::move(max_refinement)),                 // NOLINT
      blocks_per_dim_(                                // NOLINT
          std::move(blocks_per_dim)),                 // NOLINT
      grid_points_(std::move(grid_points)),
      number_of_blocks_by_dim_(blocks_per_dim_),
      boundary_condition_(std::move(boundary_condition)) {
  using domain::BoundaryConditions::is_none;
  if (is_none(boundary_condition_)) {
    PARSE_ERROR(
        context,
        "None boundary condition is not supported. If you would like an "
        "outflow-type boundary condition, you must use that.");
  }

  using domain::BoundaryConditions::is_periodic;
  if (is_periodic(boundary_condition_)) {
    is_periodic_in_[0] = true;
    is_periodic_in_[1] = true;
    is_periodic_in_[2] = true;
    boundary_condition_ = nullptr;
  }

  common_construction(domain_bounds, context);
}

template <size_t VolumeDim>
Domain<VolumeDim> NestedCubes<VolumeDim>::create_domain() const {
  return rectilinear_domain<VolumeDim>(number_of_blocks_by_dim_, block_bounds_,
                                       {}, {}, is_periodic_in_, {}, false);
}

template <size_t VolumeDim>
std::vector<DirectionMap<
    VolumeDim, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
NestedCubes<VolumeDim>::external_boundary_conditions() const {
  if (boundary_condition_ == nullptr) {
    return {};
  }
  // Set boundary conditions by using the computed domain's external
  // boundaries
  const auto domain = create_domain();
  const auto& blocks = domain.blocks();
  std::vector<DirectionMap<
      VolumeDim,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      boundary_conditions{blocks.size()};
  for (size_t i = 0; i < blocks.size(); ++i) {
    for (const Direction<VolumeDim>& external_direction :
         blocks[i].external_boundaries()) {
      boundary_conditions[i][external_direction] =
          boundary_condition_->get_clone();
    }
  }
  return boundary_conditions;
}

template <size_t VolumeDim>
std::vector<std::array<size_t, VolumeDim>>
NestedCubes<VolumeDim>::initial_extents() const {
  const auto& indices =
      indices_for_rectilinear_domains(number_of_blocks_by_dim_);
  std::vector<std::array<size_t, VolumeDim>> result(
      indices.size(), make_array<VolumeDim, size_t>(grid_points_));
  return result;
}

template <size_t VolumeDim>
std::vector<std::array<size_t, VolumeDim>>
NestedCubes<VolumeDim>::initial_refinement_levels() const {
  const auto& indices =
      indices_for_rectilinear_domains(number_of_blocks_by_dim_);
  std::vector<std::array<size_t, VolumeDim>> result(indices.size());

  // Floor divide
  const int center_block_index = static_cast<int>(blocks_per_dim_) / 2;

  const auto get_refinement = [](const std::vector<int>& in) {
    return std::max(*std::min_element(std::begin(in), std::end(in)), 0);
  };

  const auto get_distance = [this,
                             &center_block_index](const int coord) -> int {
    return static_cast<int>(max_refinement_) -
           std::abs(coord - center_block_index);
  };

  size_t index = 0;
  std::vector<int> distances(VolumeDim);

  if constexpr (VolumeDim == 3) {
    for (int z = 0; z < static_cast<int>(blocks_per_dim_); z++) {
      for (int y = 0; y < static_cast<int>(blocks_per_dim_); y++) {
        for (int x = 0; x < static_cast<int>(blocks_per_dim_); x++) {
          distances[0] = get_distance(x);
          distances[1] = get_distance(y);
          distances[2] = get_distance(z);

          result[index] = make_array<VolumeDim, size_t>(
              static_cast<size_t>(get_refinement(distances)));
          ++index;
        }
      }
    }
  } else if constexpr (VolumeDim == 2) {
    for (int y = 0; y < static_cast<int>(blocks_per_dim_); y++) {
      for (int x = 0; x < static_cast<int>(blocks_per_dim_); x++) {
        distances[0] = get_distance(x);
        distances[1] = get_distance(y);

        result[index] = make_array<VolumeDim, size_t>(
            static_cast<size_t>(get_refinement(distances)));
        ++index;
      }
    }
  } else if constexpr (VolumeDim == 1) {
    for (int x = 0; x < static_cast<int>(blocks_per_dim_); x++) {
      distances[0] = get_distance(x);

      result[index] = make_array<VolumeDim, size_t>(
          static_cast<size_t>(get_refinement(distances)));
      ++index;
    }
  }

  return result;
}  // namespace domain::creators

template class NestedCubes<1>;
template class NestedCubes<2>;
template class NestedCubes<3>;
}  // namespace domain::creators
