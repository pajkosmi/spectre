// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <memory>

#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Informer/InfoFromBuild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace radiationtransport::AnalyticData {
namespace {

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.RadiationTransport.M1Grey."
    "RotationallySymmetric",
    "[Unit][PointwiseFunctions]") {
    1 == 1;
    }

}  // namespace
}  // namespace radiationtransport::AnalyticData
