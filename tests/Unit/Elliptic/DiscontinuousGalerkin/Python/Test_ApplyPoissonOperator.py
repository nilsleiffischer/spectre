# Distributed under the MIT License.
# See LICENSE.txt for details.

import spectre.elliptic.dg as elliptic_dg
import unittest
import numpy as np
import numpy.testing as npt


class TestApplyPoissonOperator(unittest.TestCase):
    def test_apply_poisson_operator(self):
        elliptic_dg.apply_poisson_operator_flat_cartesian_1d()


if __name__ == '__main__':
    unittest.main(verbosity=2)
