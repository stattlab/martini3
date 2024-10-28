import numpy as np
from martini3 import particles
from martini3 import force_fields
import gsd.hoomd

# Cleanest way to do  I think...
# 1. write table when making force which is bond.idx and forces
# 2. Save table next to init.gsd (or any gsd's)
# 3. Run sim.create_state_from_gsd(filename=file_name + "_init.gsd")
# 4. for all bond.idx, angle.idx, dihedral.idx, create the Potentials
# 5. For LJ, coulomb, iter over particle.types

def wrap_pbc(box,x):
        np_box = np.asarray(box[0:3])
        delta = np.where(x > 0.5 * np_box, x - np_box, x)
        delta = np.where(delta < - 0.5 * np_box, np_box + delta, delta)
        return delta

def init_cell(contents, name, box_size=20, pair_on=True,file_name = "init", return_dpd=False):
    """
    initialize a simulation given the list of beads to place into it with correct bonding and nonbonding potentials.

    Args:
      contents (list of molecules): List of all molecules in simulation.
      name (string): Path to file where csv's/gsd are stored
      box_size (int/[a,b,c]): if int, square box. Otherwise, the box lengths. Only supports square boxes.
      pair_on (bool): if pair_on = false, do not compute the pair potentials (saves a few seconds if not needed)
      file_name: name of saved gsd
      return_dpd (bool): if True, will also return a hoomd.md.pair.DPD
          potential. Defaults to False. Useful when relaxing a high energy and/or random state

    Returns:
      hoomd.md.pair.LJ: Contains all LJ pair potentials present in simulation.
      hoomd.md.pair.ReactionField: Contains all Coulomb pair potentials present in simulation.
      hoomd.md.bond.Harmonic: Contains all bond  potentials present in simulation.
      hoomd.md.angle.Harmonic: Contains angle dihedral  potentials present in simulation.
      hoomd.md.dihedral.OPLS: Contains all ddihedral  potentials present in simulation.
      hoomd.md.improper.Harmonic: Contains all improper  potentials present in simulation.
      hoomd.md.constrain.Rigid: Contains all improper  potentials present in simulation.

    """
    # initialize frame
    frame = gsd.hoomd.Frame()

    # box set-up
    try:
        box = []
        box.append(float(box_size[0]))
        box.append(float(box_size[1]))
        box.append(float(box_size[2]))
        box.append(float(0))
        box.append(float(0))
        box.append(float(0))
        frame.configuration.box = box
    except:
        L = box_size
        frame.configuration.box = [L, L, L, 0, 0, 0]

    # types set-up (used for instantiating coulomb and lj potentials where interactons are based on the types of beads present)
    types = particles.init_particles(contents.bead_types)
    name_list = [particle.name for particle in types]
    mass_list = [particle.mass for particle in types]

    # Each bead type ('w','Q1' used in the simulation gets assigned an index. Here,
    # I go through all the beads present and make sure they have the correct index assigned)
    index_list = range(len(name_list))

    # flatten out contents to get particle specific parameters
    type_id = []
    positions = []
    charges = []
    mass_pairing = []
    bodies = []
    orientations = []
    moment_inertias = []
    for molecule in contents.contents:
        for bead in molecule.bead_types:
            for i in index_list:
                if bead == name_list[i]:
                    type_id.append(i)
                    mass_pairing.append(mass_list[i])
        for position in molecule.positions:
            positions.append(position)
        for charge in molecule.charges:
            charges.append(charge)
        for body in molecule.body:
            bodies.append(body)
        for orientation in molecule.orientation:
            orientations.append(orientation)
            
        for moment_inertia in molecule.moment_inertia:
            moment_inertias.append(moment_inertia)
    frame.particles.N = len(positions)
    frame.particles.position = wrap_pbc(frame.configuration.box,positions)
    frame.particles.typeid = type_id
    frame.particles.types = name_list
    frame.particles.mass = mass_pairing
    frame.particles.diameter = [0.3] * len(positions)
    frame.particles.charge = charges
    frame.particles.body = bodies
    frame.particles.orientation = orientations
    frame.particles.moment_inertia = moment_inertias

    # Get bond specific parameters for hoomd simulation
    bond_types = []
    bond_type_id = []
    bond_group = []
    angle_types = []
    angle_type_id = []
    angle_group = []
    dihedral_types = []
    dihedral_type_id = []
    dihedral_group = []
    improper_dihedral_types = []
    improper_dihedral_type_id = []
    improper_dihedral_group = []
    for molecule in contents.contents:
        for bond in molecule.bonds:
            if str(bond.idx) not in bond_types:
                bond_types.append(str(bond.idx))
            bond_type_id.append(bond.idx)
            bond_group.append(bond.bead_indices)
        for angle in molecule.angles:
            if str(angle.idx) not in angle_types:
                angle_types.append(str(angle.idx))
            angle_type_id.append(angle.idx)
            angle_group.append(angle.bead_indices)
        for dihedral in molecule.dihedrals:
            if str(dihedral.idx) not in dihedral_types:
                dihedral_types.append(str(dihedral.idx))
            dihedral_type_id.append(dihedral.idx)
            dihedral_group.append(dihedral.bead_indices)
        for improper_dihedral in molecule.improper_dihedrals:
            if str(improper_dihedral.idx) not in improper_dihedral_types:
                improper_dihedral_types.append(str(improper_dihedral.idx))
            improper_dihedral_type_id.append(improper_dihedral.idx)
            improper_dihedral_group.append(improper_dihedral.bead_indices)

    frame.bonds.N = len(bond_type_id)
    frame.bonds.types = bond_types
    frame.bonds.typeid = bond_type_id
    frame.bonds.group = bond_group

    frame.angles.N = len(angle_type_id)
    frame.angles.types = angle_types
    frame.angles.typeid = angle_type_id
    frame.angles.group = angle_group

    frame.dihedrals.N = len(dihedral_type_id)
    frame.dihedrals.types = dihedral_types
    frame.dihedrals.typeid = dihedral_type_id
    frame.dihedrals.group = dihedral_group

    frame.impropers.N = len(improper_dihedral_type_id)
    frame.impropers.types = improper_dihedral_types
    frame.impropers.typeid = improper_dihedral_type_id
    frame.impropers.group = improper_dihedral_group

    with gsd.hoomd.open(name=name + file_name + ".gsd", mode="w") as f:
        f.append(frame)

    if return_dpd:
        (
            lj,
            coulomb,
            bond_harmonic,
            angle_bonding,
            dihedral_bonding, improper_dihedral_bonding,rigid,
            dpd,
        ) = force_fields.init_all_potentials(types, contents, name, pair_on, return_dpd)
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid,dpd
    else:
        (
            lj,
            coulomb,
            bond_harmonic,
            angle_bonding,
            dihedral_bonding, improper_dihedral_bonding,rigid,
        ) = force_fields.init_all_potentials(types, contents, name, pair_on, return_dpd)
        return lj, coulomb, bond_harmonic, angle_bonding, dihedral_bonding,improper_dihedral_bonding,rigid
