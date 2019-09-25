// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <boost/python.hpp>

namespace gr {
namespace Solutions {
namespace py_bindings {
void bind_tov();
void bind_tov_isotropic();
}  // namespace py_bindings
}  // namespace Solutions
}  // namespace gr

BOOST_PYTHON_MODULE(_GeneralRelativitySolutions) {
  Py_Initialize();
  gr::Solutions::py_bindings::bind_tov();
  gr::Solutions::py_bindings::bind_tov_isotropic();
}
