// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Elliptic/Systems/Xcts/Equations.hpp"

namespace Poisson {
template <size_t Dim>
struct FirstOrderCorrectionSystem;
namespace Solutions {
template <size_t Dim>
struct ProductOfSinusoids;
}  // namespace Solutions
}  // namespace Poisson

namespace Xcts {
template <size_t Dim, Equations EnabledEquations>
struct FirstOrderSystem;
namespace Solutions {
struct ConstantDensityStar;
}
}  // namespace Xcts

template <typename System, typename InitialGuess>
struct Metavariables;
