// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/RotationallySymmetric.tpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {

template <size_t Dim>
SphericalGaussian<Dim>::SphericalGaussian(CkMigrateMessage* msg)
    : evolution::initial_data::InitialData(msg) {}

template <size_t Dim>
SphericalGaussian<Dim>::SphericalGaussian(
    const double radius, const double emissivity_and_opacity_amplitude,
    const double outer_radius, const double outer_opacity)
    : radius_(radius),
      emissivity_and_opacity_amplitude_(emissivity_and_opacity_amplitude),
      outer_radius_(outer_radius),
      outer_opacity_(outer_opacity) {
  if (UNLIKELY(radius_ > outer_radius_)) {
    ERROR("Radius " << radius_
                    << " is greater than the  "
                       "outer radius: "
                    << outer_radius_);
  }
}

template <size_t Dim>
std::unique_ptr<evolution::initial_data::InitialData>
SphericalGaussian<Dim>::get_clone() const {
  return std::make_unique<SphericalGaussian<Dim>>(*this);
}

template <size_t Dim>
void SphericalGaussian<Dim>::pup(PUP::er& p) {
  evolution::initial_data::InitialData::pup(p);
  p | radius_;
  p | emissivity_and_opacity_amplitude_;
  p | outer_radius_;
  p | outer_opacity_;
}

template <size_t Dim>
PUP::able::PUP_ID SphericalGaussian<Dim>::my_PUP_ID = 0;

template <size_t LocalDim>
bool operator==(const SphericalGaussian<LocalDim>& lhs,
                const SphericalGaussian<LocalDim>& rhs) {
  return lhs.radius_ == rhs.radius_ and
         lhs.emissivity_and_opacity_amplitude_ ==
             rhs.emissivity_and_opacity_amplitude_ and
         lhs.outer_radius_ == rhs.outer_radius_ and
         lhs.outer_opacity_ == rhs.outer_opacity_;
}

template <size_t Dim>
bool operator!=(const SphericalGaussian<Dim>& lhs,
                const SphericalGaussian<Dim>& rhs) {
  return not(lhs == rhs);
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::energy_profile(const DataVector& r) const {
  // Just a smooth function.  Not particularly close to the actual
  // solution.
  return exp(-square(1.0 / radius_) * square(r));
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::emission_profile(const DataVector& r) const {
  return emissivity_and_opacity_amplitude_ *
         exp(-square(1.0 / radius_) * square(r));
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::absorption_profile(
    const DataVector& r) const {
  return emissivity_and_opacity_amplitude_ *
             exp(-square(1.0 / radius_) * square(r)) +
         outer_opacity_ * step_function(r - outer_radius_);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                \
  template class SphericalGaussian<DIM(data)>;                              \
  template bool operator==<DIM(data)>(const SphericalGaussian<DIM(data)>&,  \
                                      const SphericalGaussian<DIM(data)>&); \
  template bool operator!=<DIM(data)>(const SphericalGaussian<DIM(data)>&,  \
                                      const SphericalGaussian<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (2, 3))
#undef INSTANTIATE
#undef DIM

INSTANTIATE_ROTATIONALLY_SYMMETRIC(SphericalGaussian<2>)
INSTANTIATE_ROTATIONALLY_SYMMETRIC(SphericalGaussian<3>)
}  // namespace RadiationTransport::M1Grey::AnalyticData
