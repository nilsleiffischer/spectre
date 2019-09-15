// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/Xcts/NeutronStarBinary.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "ErrorHandling/Error.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/TovStar.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace Xcts {
namespace AnalyticData {

NeutronStarBinary::NeutronStarBinary(
    const double separation, const double eccentricity,
    const std::array<double, 2> central_rest_mass_densities,
    const double polytropic_constant, const double polytropic_exponent) noexcept
    : separation_(separation),
      eccentricity_(eccentricity),
      central_rest_mass_densities_(std::move(central_rest_mass_densities)),
      polytropic_constant_(polytropic_constant),
      polytropic_exponent_(polytropic_exponent),
      equation_of_state_{polytropic_constant_, polytropic_exponent_},
      stars_{{{central_rest_mass_densities_[0], polytropic_constant,
               polytropic_exponent},
              {central_rest_mass_densities_[1], polytropic_constant,
               polytropic_exponent}}} {}

void NeutronStarBinary::pup(PUP::er& p) noexcept {
  p | separation_;
  p | eccentricity_;
  p | central_rest_mass_densities_;
  p | polytropic_constant_;
  p | polytropic_exponent_;
  p | equation_of_state_;
  p | stars_;
}

bool operator==(const NeutronStarBinary& lhs,
                const NeutronStarBinary& rhs) noexcept {
  // there is no comparison operator for the EoS, but should be okay as
  // the `polytropic_exponent`s and `polytropic_constant`s are compared
  return lhs.central_rest_mass_densities_ ==
             rhs.central_rest_mass_densities_ and
         lhs.separation_ == rhs.separation_ and
         lhs.eccentricity_ == rhs.eccentricity_ and
         lhs.polytropic_constant_ == rhs.polytropic_constant_ and
         lhs.polytropic_exponent_ == rhs.polytropic_exponent_;
}

bool operator!=(const NeutronStarBinary& lhs,
                const NeutronStarBinary& rhs) noexcept {
  return not(lhs == rhs);
}

}  // namespace AnalyticData
}  // namespace Xcts
/// \endcond
