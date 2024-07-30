import pytest
import sys
import os
import hoomd.md

sys.path.append("../")

from martini3 import force_fields as f
from martini3 import particles as p


@pytest.fixture
def p6_particle_list():
    return p.init_particles(["P6"])


@pytest.fixture
def q5_particle_list():
    return p.init_particles(["Q5"])


@pytest.fixture
def cell():
    return hoomd.md.nlist.Cell(buffer=0.4)


@pytest.fixture
def q5_p6_particle_list():
    return p.init_particles(["Q5", "P6", "SP6", "W", "SN4a"])


def test_LJ_P6(p6_particle_list, cell):
    assert (
        f.init_lj_potentials(p6_particle_list, cell).params[("P6", "P6")]["epsilon"]
        == 4.99
    )


def test_LJ_Q5(q5_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_particle_list, cell).params[("Q5", "Q5")]["epsilon"]
        == 6.45
    )


def test_LJ_Q5_P6_1(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "Q5")]["epsilon"]
        == 6.45
    )


def test_LJ_Q5_P6_2(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("P6", "P6")]["epsilon"]
        == 4.99
    )


def test_LJ_Q5_P6_3_lower(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "P6")]["epsilon"]
        > 4.99
    )


def test_LJ_Q5_P6_3_upper(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "P6")]["epsilon"]
        < 10
    )


def test_LJ_Q5_P6_4_lower(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "SP6")]["epsilon"]
        > 4.99
    )


def test_LJ_Q5_P6_4_upper(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "SP6")]["epsilon"]
        < 9
    )


def test_LJ_Q5_P6_5_lower(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "SN4a")][
            "epsilon"
        ]
        > 2
    )


def test_LJ_Q5_P6_5_upper(q5_p6_particle_list, cell):
    assert (
        f.init_lj_potentials(q5_p6_particle_list, cell).params[("Q5", "SN4a")][
            "epsilon"
        ]
        < 9
    )
