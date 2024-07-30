# Contains all particles and their interrelation matrix.
import csv
import os

script_path = os.path.abspath(__file__)
root = script_path.split("/martini3")[0]


class Particle:
    def __init__(self, name, mass, lj, charge=None):
        self.name = name
        self.mass = mass
        self.charge = charge
        self.lj_params = lj


def nth_particle(row_text, item_number):
    """
    Splits the text in a row of text into an array, finding the nth piece of information from that row

    Args:
      row_text: One row of text document (in this case the itp)
      item_number (int): if mass, item number is one. Otherwise it is whatever.

    Returns:
      item (str): nth piece of text present in row
    """
    iter = 0
    for item in row_text:
        if item != "" and iter == item_number:
            return item
        elif item != "":
            iter = iter + 1


def mass_lookup(particle_name):
    """
    Returns the mass of a particle given the particle name. Uses Martini general mass (R = 72, S = 54? T = 36?).

    Args:
      particle_name (string): name of particle

    Returns:
      mass (float): mass of particle
    """
    masses = root + "/martini3/" + "martini_v3.0.0.itp"
    if particle_name == "cholesterol":
        mass = 432
    else:
        with open(masses, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                if row != []:
                    row_text = row[0].split(" ")
                if row_text[0] == particle_name:
                    mass = float(nth_particle(row_text, 1))
                    break
    return mass

def cholesterol_lj(particle_names):
    """
    Creates non-interacting leonard jones potentials between cholesterol and all particles in simulation

    Args:
      particle names (list of strings)

    Returns:
      lj_potentials (dict):dict of all particle names and the tuple of LJ potentials
    """
    #there is no interaction with the central rigid body cholesterol bead
    lj_potentials = {}
    for particle in particle_names:
        lj_potentials.update({particle: (0, 0)})
    return lj_potentials

def lj_lookup(particle_name, particle_names):
    """
    For a given particle, a dict of all the other particle names (key) and their lj params (sigma,epsilon) with the given particle

    Args:
      particle name (string)
      particle names (list of strings)

    Returns:
      lj_potentials (dict): dict of all particle names and the tuple of LJ potentials
    """
    forces = root + "/martini3/" + "martini_v3.0.0.itp"
    # lj  potentials is a dict with dict{key = name of other particle: value = (sigma,eps)}
    with open(forces, newline="") as f:
        reader = csv.reader(f)
        lj_potentials = {}
        nonbond_params_point = False
        for row in reader:
            if row != []:
                if row[0] == "[ nonbond_params ]":
                    nonbond_params_point = True
                if nonbond_params_point:
                    row_text = row[0].split(" ")
                    p1_name = nth_particle(row_text, 0)
                    p2_name = nth_particle(row_text, 1)
                    if p1_name == particle_name:
                        if p2_name in particle_names:
                            sigma = float(nth_particle(row_text, 3))
                            eps = float(nth_particle(row_text, 4))
                            lj_potentials.update({p2_name: (sigma, eps)})
                    if p2_name == particle_name:
                        if p1_name in particle_names:
                            sigma = float(nth_particle(row_text, 3))
                            eps = float(nth_particle(row_text, 4))
                            lj_potentials.update({p1_name: (sigma, eps)})
    if particle_name =="cholesterol":
        lj_potentials =  cholesterol_lj(particle_names)
    elif "cholesterol" in particle_names:
        lj_potentials.update({"cholesterol": (0, 0)})
    return lj_potentials


def init_particles(particle_list_string):
    """
    Creates list of particle (type particle) for all particle types in simulation

    Args:
      particle_list_string (list of strings): list of all unique particle types in a simulation

    Returns:
      particle_list (list of Particles): Particle is a class containing a name, mass, and dict of all leonard jones interaction potentials
    """
    particle_names = particle_list_string
    particle_list = []
    for i, particle in enumerate(particle_names):
        mass = mass_lookup(particle)
        lj = lj_lookup(particle, particle_names)
        particle_list.append(Particle(particle, mass, lj))

    return particle_list
