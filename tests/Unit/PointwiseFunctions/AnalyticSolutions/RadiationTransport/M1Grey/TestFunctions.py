# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np
import scipy.optimize as opt


# Functions for testing ConstantM1.cpp
def constant_m1_spatial_velocity(x, t, mean_velocity, comoving_energy_density):
    return np.asarray(mean_velocity)


def constant_m1_lorentz_factor(x, t, mean_velocity, comoving_energy_density):
    return 1.0 / np.sqrt(1.0 - np.linalg.norm(np.asarray(mean_velocity)) ** 2)


def constant_m1_tildeE(x, t, mean_velocity, comoving_energy_density):
    w_sqr = 1.0 / (1.0 - np.linalg.norm(np.asarray(mean_velocity)) ** 2)
    return comoving_energy_density / 3.0 * (4.0 * w_sqr - 1.0)


def constant_m1_tildeS(x, t, mean_velocity, comoving_energy_density):
    w_sqr = 1.0 / (1.0 - np.linalg.norm(np.asarray(mean_velocity)) ** 2)
    prefactor = 4.0 / 3.0 * comoving_energy_density * w_sqr
    return np.asarray(mean_velocity) * prefactor


# End Functions for testing ConstantM1.cpp


def homogen_sphere_m1_tildeE(
    x, t, radius, emissivity_and_opacity, outer_radius, outer_opacity
):
    # how sharp/rounded the edges of the sphere are
    # the closer to 0 this becomes, the sharper the discontinuity
    sharpness = -0.03

    radii = np.linalg.norm(np.asarray(x))
    normalized_radii = (radii - radius) / sharpness

    energy_difference = 1.0 - 1.0e-12
    energy_sum = 1.0 + 1.0e-12

    e_tilde = (
        energy_difference / np.pi * np.arctan(normalized_radii)
        + 0.5 * energy_sum
    )

    return e_tilde


def homogen_sphere_m1_tildeS(x, t, mean_velocity, comoving_energy_density):
    return np.asarray(mean_velocity) * 0.0


# End Functions for testing HomogeneousSphere
