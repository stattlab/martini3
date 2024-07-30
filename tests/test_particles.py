import pytest
import sys
import os

sys.path.append("../")

from martini3 import particles


def test_particle_constructor_name():
    assert particles.Particle("Test", 2, (0, 0), charge=1).name == "Test"


def test_particle_constructor_mass():
    assert particles.Particle("Test", 2, (0, 0), charge=1).mass == 2


def test_particle_constructor_charge():
    assert particles.Particle("Test", 2, (0, 0), charge=1).charge == 1


def test_particle_constructor_lj():
    assert particles.Particle("Test", 2, (0, 0), charge=1).lj_params == (0, 0)


def test_mass_lookup_Q4():
    assert particles.mass_lookup("Q4") == 72.0


def test_lj_lookup_P6_TC1():
    particle_name = "P6"
    particle_names = ["TC1"]
    assert particles.lj_lookup(particle_name, particle_names) == {
        "TC1": (4.650000e-01, 7.030000e-01)
    }


def test_lj_lookup_P6_list():
    particle_name = "P6"
    particle_names = ["TC1", "TP1d", "P6"]
    assert particles.lj_lookup(particle_name, particle_names) == {
        "TC1": (4.650000e-01, 7.030000e-01),
        "TP1d": (0.395, 3.19),
        "P6": (0.47, 4.99),
    }


def test_init_particles():
    particle_names = ["P6"]
    particle_list = particles.init_particles(particle_names)
    assert (
        particle_list[0].name == "P6"
        and particle_list[0].mass == 72
        and particle_list[0].lj_params == {"P6": (0.47, 4.99)}
    )
