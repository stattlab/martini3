import hoomd

# from particles import Particle
from martini3 import q_adj
from martini3 import particles
import numpy as np
import csv
import gsd.hoomd

def init_dpd_potentials(types, cell):
    """
    initialize dpd potentials of all types of beads in your simulation
        according to the sigma of Martini3 specifications. All pairs have an interaction 
        strength of 20, gamma 4.5, and assumes thermostat temperature of 1.0
    
    Useful when relaxing a high energy/random state

    Args:
      types (list of Particles): List of all particles types in simulation.
      cell (hoomd neighborlist): HOOMD neighborlist

    Returns:
      hoomd.md.pair.dpd: Pair potential between all particle types.
    """
    dpd = hoomd.md.pair.DPD(nlist=cell,kT=1.0)

    for i in range(len(types)):
        for j in range(len(types)):
            particle_a = types[i]
            particle_b = types[j]
            dpd.params[(particle_a.name, particle_b.name)] = dict(
                A=20, gamma=4.5,
            )
            dpd.r_cut[(particle_a.name, particle_b.name)] = \
                                    particle_a.lj_params[particle_b.name][0]  # nm
    return dpd

def init_lj_potentials(types, cell):
    """
    initialize the lj potentials of all types of beads in your simulation
    according to the Martini3 specifications

    Args:
      types (list of Particles): List of all particles types in simulation.
      cell (hoomd neighborlist): HOOMD neighborlist

    Returns:
      hoomd.md.pair.LJ: Pair potential between all particle types.
    """
    lj = hoomd.md.pair.LJ(nlist=cell)

    for i in range(len(types)):
        for j in range(len(types)):
            particle_a = types[i]
            particle_b = types[j]

            # Adjust for Q bead interactions if one of the beads is Q
            if (particle_a.lj_params[particle_b.name][0] !=0) and "Q" in (particle_a.name) and "Q" not in (particle_b.name):
                if q_adj.eps_b.get(particle_b.name) != None:
                    eps_b = q_adj.eps_b.get(particle_b.name)
                    b_size = "R"
                else:
                    if q_adj.eps_b.get(particle_b.name[0:2]) != None:
                        b_size = "R"
                        eps_b = q_adj.eps_b.get(particle_b.name[0:2])
                    elif q_adj.eps_b.get(particle_b.name[1:3]) != None:
                        b_size = particle_b.name[0]
                        eps_b = q_adj.eps_b.get(particle_b.name[1:3])
                    else:
                        print(
                            "error lipid.py line 78 idk what would cause this case but writing this "
                        )

                # This line assumes that particle a is a Q bead and not a D bead and that it has no letter modifiers
                if len(particle_a.name) > 2:
                    a_size = particle_a.name[0]
                    a_name = particle_a.name[1:3]
                else:
                    a_size = "R"
                    a_name = particle_a.name

                eps_qb = particle_a.lj_params[particle_b.name][1]
                eps_w = q_adj.eps_b.get("W")
                eps_c1 = q_adj.eps_b.get("C1")
                gamma = q_adj.gamma.get(a_name).get(b_size)
                eps_inside_w = q_adj.eps_w.get(a_name).get(b_size)
                eps_inside_c1 = q_adj.eps_c1.get(a_name).get(b_size)
                p_qb = q_adj.p_qb.get(a_size).get(b_size)

                sigma_qw = particle_a.lj_params["W"][0]
                sigma_qb = particle_a.lj_params[particle_b.name][0]

                eps_b_inside = eps_inside_w + (eps_b - eps_w) / (eps_c1 - eps_w) * (
                    eps_inside_c1 - eps_inside_w
                )

                eps_final = (
                    eps_qb
                    + p_qb
                    * gamma
                    * sigma_qw
                    / sigma_qb
                    * eps_w
                    * eps_inside_w
                    / eps_b
                    / eps_b_inside
                    * (eps_b - eps_b_inside)
                    / (eps_w - eps_inside_w)
                )

                lj.params[(particle_a.name, particle_b.name)] = dict(
                    epsilon=eps_final, sigma=particle_a.lj_params[particle_b.name][0]
                )
                lj.r_cut[(particle_a.name, particle_b.name)] = 1.1  # nm
            elif (particle_a.lj_params[particle_b.name][0] !=0)  and "Q" not in (particle_a.name) and "Q" in (particle_b.name):
                if q_adj.eps_b.get(particle_a.name) != None:
                    eps_b = q_adj.eps_b.get(particle_a.name)
                    a_size = "R"
                else:
                    if q_adj.eps_b.get(particle_a.name[0:2]) != None:
                        a_size = "R"
                        eps_b = q_adj.eps_b.get(particle_a.name[0:2])
                    elif q_adj.eps_b.get(particle_a.name[1:3]) != None:
                        a_size = particle_a.name[0]
                        eps_b = q_adj.eps_b.get(particle_a.name[1:3])
                    else:
                        print(
                            "error lipid.py line 78 idk what would cause this case but writing this "
                        )

                # This line assumes that particle b is a Q bead and not a D bead
                # this also assumes Q beads have no name modifiers
                if len(particle_b.name) > 2:
                    b_size = particle_b.name[0]
                    b_name = particle_b.name[1:3]
                else:
                    b_size = "R"
                    b_name = particle_b.name

                eps_qb = particle_b.lj_params[particle_a.name][1]
                eps_w = q_adj.eps_b.get("W")
                eps_c1 = q_adj.eps_b.get("C1")
                gamma = q_adj.gamma.get(b_name).get(a_size)
                eps_inside_w = q_adj.eps_w.get(b_name).get(a_size)
                eps_inside_c1 = q_adj.eps_c1.get(b_name).get(a_size)
                p_qb = q_adj.p_qb.get(b_size).get(a_size)

                sigma_qw = particle_b.lj_params["W"][0]
                sigma_qb = particle_b.lj_params[particle_a.name][0]

                eps_b_inside = eps_inside_w + (eps_b - eps_w) / (eps_c1 - eps_w) * (
                    eps_inside_c1 - eps_inside_w
                )

                eps_final = (
                    eps_qb
                    + p_qb
                    * gamma
                    * sigma_qw
                    / sigma_qb
                    * eps_w
                    * eps_inside_w
                    / eps_b
                    / eps_b_inside
                    * (eps_b - eps_b_inside)
                    / (eps_w - eps_inside_w)
                )

                lj.params[(particle_a.name, particle_b.name)] = dict(
                    epsilon=eps_final, sigma=particle_a.lj_params[particle_b.name][0]
                )
                lj.r_cut[(particle_a.name, particle_b.name)] = 1.1  # nm
            else:
                lj.params[(particle_a.name, particle_b.name)] = dict(
                    epsilon=particle_a.lj_params[particle_b.name][1],
                    sigma=particle_a.lj_params[particle_b.name][0],
                )
                lj.r_cut[(particle_a.name, particle_b.name)] = 1.1  # nm
    return lj


def init_coulomb_potentials(types, cell):
    """
    initialize the coulomb potentials of all charged species of beads in your simulation
    according to the Martini3 specifications

    Args:
      types (list of Particles or list of strings): List of all particles types in simulation.
      cell (hoomd neighborlist): hoomd neighborlist

    Returns:
      hoomd.md.pair.ReactionField: Pair potential between all particle types using ReactionField Method implementation.
    """

    # eps = eps_0*eps_r
    # U = -charge/4pi*eps_0*r
    # F = charge/4pi*eps_0*r**2
    # for our simulation model, charge = sqrt(4pi*eps_0)
    eps_r = 15
    coulomb = hoomd.md.pair.ReactionField(cell, default_r_cut=1.1)
    for type_i in types:
        for type_j in types:
            name_i = type_i.name
            name_j = type_j.name

            if "Q" in (name_i) and "Q" in (name_j):
                coulomb.params[(name_i, name_j)] = dict(
                    epsilon=1.0, eps_rf=eps_r, use_charge=True
                )
            else:
                coulomb.params[(name_i, name_j)] = dict(epsilon=0, eps_rf=eps_r)

    return coulomb


def init_harmonic_bonds(contents, name):
    """
    initialize the bonded potentials of all bonds present in simulation
    according to the Martini3 specifications. Also writes identities of bonds 
    to bonds.csv

    Args:
      contents (list of Molecules): List of all molecules in simulation.
      name (string): path to folder where gsd is saved.

    Returns:
      hoomd.md.bond.Harmonic: Contains all bonds present in simulation.
    """
    bond_harmonic = hoomd.md.bond.Harmonic()

    bond_set = set()
    for molecule in contents.contents:
        for bond in molecule.bonds:
            if type(bond_harmonic.params[str(bond.idx)].get("k")) != type(10):
                bond_harmonic.params[str(bond.idx)] = dict(
                    k=bond.force, r0=bond.spatial
                )
                bond_set.add((bond.idx, bond.force, bond.spatial))

    with open(name + "bonds.csv", "w") as file:
        writer = csv.writer(file)
        for bond in bond_set:
            writer.writerow([bond[0], bond[1], bond[2]])
    return bond_harmonic


def init_angles(contents, name):
    """
    initialize the angled potentials of all angles present in simulation
    according to the Martini3 specifications. Also saves identities of angles
    to angles.csv

    Args:
      contents (list of molecules): List of all molecules in simulation.
      name (string): path to folder where gsd is saved.

    Returns:
      hoomd.md.angle.Harmonic: Contains all angles present in simulation.
    """
    angle_bonding = hoomd.md.angle.CosineSquared()
    angle_set = set()
    for molecule in contents.contents:
        for angle in molecule.angles:
            if type(angle_bonding.params[str(angle.idx)].get("k")) != type(10):
                angle_bonding.params[str(angle.idx)] = dict(
                    k=angle.force, t0=float(angle.spatial) / 180 * np.pi
                )
                angle_set.add(
                    (angle.idx, angle.force, float(angle.spatial) / 180 * np.pi)
                )
    with open(name + "angles.csv", "w") as file:
        writer = csv.writer(file)
        for angle in angle_set:
            writer.writerow([angle[0], angle[1], angle[2]])
    return angle_bonding


def init_dihedrals(contents, name):
    """
    initialize the OPLS dihedral potentials of all angles present in simulation
    according to the Martini3 specifications. also saves iddentities of dihedrals to dihedrals.csv

    Args:
      contents (list of molecules): List of all molecules in simulation.
      name (string): path to folder where gsd is saved.

    Returns:
      hoomd.md.diheadral.Periodic: Contains all dihedrals present in simulation.
    """
    dihedral_set = set()
    dihedral_bonding = hoomd.md.dihedral.OPLS()
    for molecule in contents.contents:
        for dihedral in molecule.dihedrals:
            if type(dihedral_bonding.params[str(dihedral.idx)].get("k1")) != type(10):
                dihedral_bonding.params[str(dihedral.idx)] = dict(
                    k1=dihedral.k1, k2=dihedral.k2, k3=dihedral.k3, k4=dihedral.k4
                )
                dihedral_set.add(
                    (dihedral.idx, dihedral.k1, dihedral.k2, dihedral.k3, dihedral.k4)
                )
    with open(name + "dihedrals.csv", "w") as file:
        writer = csv.writer(file)
        for dihedral in dihedral_set:
            writer.writerow(
                [dihedral[0], dihedral[1], dihedral[2], dihedral[3], dihedral[4]]
            )
    return dihedral_bonding

def init_improper_dihedrals(contents, name):
    """
    initialize the improper dihedral potentials of all dihdeals present in simulation
    according to the Martini3 specifications. Also writes identiteis of improper dihedrals
    to improper_dihedrals.csv

    Args:
      contents (list of molecules): List of all molecules in simulation.
      name (string): path to folder where gsd is saved.

    Returns:
      hoomd.md.improper.Harmonic: Contains all impropers present in simulation.
    """
    improper_dihedral_set = set()
    improper_dihedral_bonding = hoomd.md.improper.Harmonic()
    for molecule in contents.contents:
        for improper_dihedral in molecule.improper_dihedrals:
            # psome weirdness made me drop the if statement (which I believe was protecting us from writing the same thing many times, but I have decided that is fine for these)
            improper_dihedral_bonding.params[str(improper_dihedral.idx)] = dict(
                k=improper_dihedral.k, chi0=str(float(improper_dihedral.chi0)*np.pi/180)
            )
            improper_dihedral_set.add(
                (improper_dihedral.idx, improper_dihedral.k, float(improper_dihedral.chi0)*np.pi/180)
            )
    with open(name + "improper_dihedrals.csv", "w") as file:
        writer = csv.writer(file)
        for improper_dihedral in improper_dihedral_set:
            writer.writerow(
                [improper_dihedral[0], improper_dihedral[1], improper_dihedral[2]]
            )
    return improper_dihedral_bonding


def read_harmonic_bonds(bonds_path):
    """
    read the bonded potentials of all bonds present in simulation from bonds.csv
    according to the Martini3 specifications

    Args:
      bonds_path (string): path to bonds.csv file.

    Returns:
      hoomd.md.bond.Harmonic: Contains all bonds present in simulation.
    """
    bond_harmonic = hoomd.md.bond.Harmonic()

    with open(bonds_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            bond_harmonic.params[str(row[0])] = dict(k=row[1], r0=row[2])

    return bond_harmonic


def read_angles(angle_path):
    """
    read the angled potentials of all angles present in simulation
    according to the Martini3 specifications

    Args:
      angle_path (string): path to angles.csv file.

    Returns:
      hoomd.md.angle.Harmonic: Contains all angles present in simulation.
    """
    angle_bonding = hoomd.md.angle.CosineSquared()

    with open(angle_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            angle_bonding.params[str(row[0])] = dict(k=row[1], t0=float(row[2]))

    return angle_bonding


def read_dihedrals(dihedral_path):
    """
    read the angled potentials of all angles present in simulation
    according to the Martini3 specifications

    Args:
      dihedral_path (string): path to dihedral.csv file.

    Returns:
      hoomd.md.diheadral.Periodic: Contains all dihedrals present in simulation.
    """
    dihedral_bonding = hoomd.md.dihedral.OPLS()
    with open(dihedral_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            dihedral_bonding.params[str(row[0])] = dict(
                k1=row[1], k2=row[2], k3=row[3], k4=row[4]
            )
    return dihedral_bonding

def read_improper_dihedrals(improper_dihedral_path):
    """
    read the improper potentials of all angles present in simulation
    according to the Martini3 specifications

    Args:
      improper_dihedral_path (string): path to improper_dihedral.csv

    Returns:
      hoomd.md.diheadral.Periodic: Contains all improper_dihedrals present in simulation.
    """
    improper_dihedral_bonding = hoomd.md.improper.Harmonic()
    with open(improper_dihedral_path, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            improper_dihedral_bonding.params[str(row[0])] = dict(
                k=row[1], chi0=row[2]
            )
    return improper_dihedral_bonding

def make_rigid():
    #Right now this is just cholesterol
    """
    Make rigid body entity for cholesterol with correct constituent beads and positions. 
    Assigns subbeads to the rigid body entity.
    Args:
    Returns:
      hoomd.md.constrain.Rigid: Contains rigid body entity
    """
    rigid = hoomd.md.constrain.Rigid()
    rigid.body['cholesterol'] = {
        'constituent_types': ["P1","SC4","SC3","SC3","SC3","TC2","TC2","C2"],
        'positions': [[ 0.02610333, -0.01810241, -0.53121775],
        [-0.22690158, -0.01555445, -0.18325152],
        [ 0.11430364, -0.08128266, -0.20192171],
        [-0.20322607, -0.01623255,  0.24130113],
        [ 0.13670035, -0.02708035,  0.1005572 ],
        [ 0.08236502,  0.16761672, -0.20705957],
        [0.00900036, 0.15953455, 0.25212949],
        [ 0.06255672, -0.04036072,  0.54116896]
                    ],
        'orientations': [(1.0, 0.0, 0.0, 0.0)]*8,
    }
    return rigid


def init_all_potentials(types, contents, name, pair_on, return_dpd=False):
    """
    initialize all potentials using functions written above. also makes bond.csv, angle.csv, etc. so the simulation can be started from a gsd.

    Args:
      types (list of types): list of all types present in simulation
      contents (list of molecules): List of all molecules in simulation.
      name (string): Path to file where csv's are stored
      pair_on (bool): if pair_on = false, do not compute the pair potentials (saves a few seconds if not needed)
      return_dpd (bool): if True, init_all_potentials will also return a hoomd.md.pair.DPD
                    potential. Defaults to False. Useful when relaxing a high energy and/or random state

    Returns:
      hoomd.md.pair.LJ: Contains all LJ pair potentials present in simulation.
      hoomd.md.pair.ReactionField: Contains all Coulomb pair potentials present in simulation.
      hoomd.md.bond.Harmonic: Contains all bond potentials present in simulation.
      hoomd.md.angle.Harmonic: Contains angle potentials present in simulation.
      hoomd.md.dihedral.OPLS: Contains all dihedral potentials present in simulation.
      hoomd.md.improper.Harmonic: Contains all improper potentials present in simulation.
      hoomd.md.constrain.Rigid: Contains all improper potentials present in simulation.
      hoomd.md.pair.DPD: Contains all DPD pair potentials present in simulation, based on LJ parameters.

    """
    cell = hoomd.md.nlist.Cell(buffer=0.4,exclusions = ('bond','body'))
    if pair_on:
        lj = init_lj_potentials(types, cell)
        coulomb = init_coulomb_potentials(types, cell)
    else:
        coulomb = hoomd.md.pair.ReactionField(cell, default_r_cut=1.1)
        lj = hoomd.md.pair.LJ(nlist=cell)
    bond_harmonic = init_harmonic_bonds(contents, name)
    angle_bonding = init_angles(contents, name)
    dihedral_bonding = init_dihedrals(contents, name)
    improper_dihedral_bonding = init_improper_dihedrals(contents, name)
    rigid = make_rigid()
    if (return_dpd):
        dpd = init_dpd_potentials(types,cell)
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid,dpd
    else:
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid


def forces_from_gsd(path, gsd_name, return_dpd=False):
    """
    initialize all potentials using functions written above.
    also makes bond.csv, angle.csv, etc. so the simulation can be started from a gsd.

    Args:
      path (string): path to file where the gsd is
      gsd_name (string): gsd file name
      return_dpd (bool): if True, init_all_potentials will also return a hoomd.md.pair.DPD
                    potential. Defaults to False. Useful when relaxing a high energy and/or random state

    Returns:
      hoomd.md.pair.LJ: Contains all LJ pair potentials present in simulation.
      hoomd.md.pair.ReactionField: Contains all Coulomb pair potentials present in simulation.
      hoomd.md.bond.Harmonic: Contains all bond potentials present in simulation.
      hoomd.md.angle.Harmonic: Contains angle potentials present in simulation.
      hoomd.md.dihedral.OPLS: Contains all dihedral potentials present in simulation.
      hoomd.md.improper.Harmonic: Contains all improper potentials present in simulation.
      hoomd.md.constrain.Rigid: Contains all improper potentials present in simulation.
      hoomd.md.pair.DPD: Contains all DPD pair potentials present in simulation, based on LJ parameters.
    """
    gsd_path = path + gsd_name
    bonds_path = path + "bonds.csv"
    dihedrals_path = path + "dihedrals.csv"
    improper_dihedrals_path = path + "improper_dihedrals.csv"
    angles_path = path + "angles.csv"
    traj = gsd.hoomd.open(gsd_path, mode="r")
    cell = hoomd.md.nlist.Cell(buffer=0.4,exclusions = ('bond','body'))
    frame = traj[-1]
    string_types = frame.particles.types
    particle_types = particles.init_particles(string_types)
    lj = init_lj_potentials(particle_types, cell)
    coulomb = init_coulomb_potentials(particle_types, cell)
    bond_harmonic = read_harmonic_bonds(bonds_path)
    angle_bonding = read_angles(angles_path)
    dihedral_bonding = read_dihedrals(dihedrals_path)
    improper_dihedral_bonding = read_improper_dihedrals(improper_dihedrals_path)
    rigid = make_rigid()
    if (return_dpd):
        dpd = init_dpd_potentials(particle_types,cell)
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid,dpd
    else:
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid