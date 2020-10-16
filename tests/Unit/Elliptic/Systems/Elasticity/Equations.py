# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def constitutive_relation_2d(strain, bulk_modulus, shear_modulus):
    lame_constant = bulk_modulus - 2. / 3. * shear_modulus
    return -2. * shear_modulus * lame_constant / (
        lame_constant + 2. * shear_modulus
    ) * np.trace(strain) * np.eye(2) - 2. * shear_modulus * strain


def constitutive_relation_3d(strain, bulk_modulus, shear_modulus):
    lame_constant = bulk_modulus - 2. / 3. * shear_modulus
    return -2. * shear_modulus * strain - lame_constant * np.trace(
        strain) * np.eye(3)


def primal_fluxes_2d(strain, coordinates, bulk_modulus, shear_modulus):
    return -constitutive_relation_2d(strain, bulk_modulus, shear_modulus)


def primal_fluxes_3d(strain, coordinates, bulk_modulus, shear_modulus):
    return -constitutive_relation_3d(strain, bulk_modulus, shear_modulus)


def add_non_euclidean_sources(christoffel_second_kind, christoffel_contracted,
                              stress):
    return (-np.einsum('i,ij', christoffel_contracted, stress) -
            np.einsum('ijk,jk', christoffel_second_kind, stress))


def auxiliary_fluxes(displacement):
    dim = len(displacement)
    # Compute the tensor product with a Kronecker delta and symmetrize the last
    # two indices.
    tensor_product = np.tensordot(np.eye(dim), displacement, axes=0)
    return 0.5 * (tensor_product + np.transpose(tensor_product, (0, 2, 1)))


def non_euclidean_auxiliary_fluxes(metric, displacement):
    co_displacement = np.einsum('ij,j', metric, displacement)
    return auxiliary_fluxes(co_displacement)
