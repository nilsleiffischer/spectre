// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <pup_stl.h>
#include <utility>
#include <variant>

#include "Utilities/TMPL.hpp"

namespace PUP {
/// \ingroup ParallelGroup
/// Serialization of std::optional for Charm++
template <typename T>
void pup(PUP::er& p, std::optional<T>& t) noexcept {  // NOLINT
  bool valid = false;
  if (p.isUnpacking()) {
    p | valid;
    if (valid) {
      T value{};
      p | value;
      t = std::move(value);
    } else {
      t.reset();
    }
  } else {
    valid = t.has_value();
    p | valid;
    if (valid) {
      p | *t;
    }
  }
}

/// \ingroup ParallelGroup
/// Serialization of std::optional for Charm++
template <typename T>
void operator|(er& p, std::optional<T>& t) noexcept {  // NOLINT
  pup(p, t);
}

/// \ingroup ParallelGroup
/// Serialization of std::variant for Charm++
template <typename... Ts>
void pup(PUP::er& p, std::variant<Ts...>& t) noexcept {  // NOLINT
  // clang tidy: complains about pointer decay, I don't understand the error
  // since nothing here is doing anything with pointers.
  assert(not t.valueless_by_exception());  // NOLINT

  size_t current_index = 0;
  size_t send_index = t.index();
  p | send_index;

  const auto pup_helper = [&current_index, &p, send_index,
                           &t](auto type_v) noexcept {
    using type = tmpl::type_from<decltype(type_v)>;
    if (current_index == send_index) {
      if (p.isUnpacking()) {
        type temp{};
        p | temp;
        t = std::move(temp);
      } else {
        p | std::get<type>(t);
      }
    }
    current_index++;
  };
  EXPAND_PACK_LEFT_TO_RIGHT(pup_helper(tmpl::type_<Ts>{}));
}

/// \ingroup ParallelGroup
/// Serialization of std::variant for Charm++
template <typename... Ts>
void operator|(PUP::er& p, std::variant<Ts...>& t) noexcept {  // NOLINT
  pup(p, t);
}
}  // namespace PUP
