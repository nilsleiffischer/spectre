// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "@CMAKE_SOURCE_DIR@/src/Informer/InfoFromBuild.hpp"

#include <boost/preprocessor.hpp>
#include <sstream>

namespace {
std::string link_date() { return std::string(__TIMESTAMP__); }

std::string git_description() {
  return std::string("@GIT_DESCRIPTION@");
}

std::string git_branch() { return std::string("@GIT_BRANCH@"); }
}  // namespace

std::string info_from_build() {
  static const std::string info = [] {
    std::ostringstream os;
    os << "SpECTRE Build Information:\n";
    os << "Version:                      " << spectre_version() << "\n";
    os << "Compiled on host:             @HOSTNAME@\n";
    os << "Compiled in directory:        @CMAKE_BINARY_DIR@\n";
    os << "Source directory is:          @CMAKE_SOURCE_DIR@\n";
    os << "Compiled on git branch:       " << git_branch() << "\n";
    os << "Compiled on git revision:     " << git_description() << "\n";
    os << "Linked on:                    " << link_date() << "\n";
    return os.str();
  }();
  return info;
}
