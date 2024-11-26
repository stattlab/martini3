import csv
import os
import sys
from martini3 import polymers
import numpy as np

script_path = os.path.abspath(__file__)
root = script_path.split("/martini3")[0]


class Bond:
    def __init__(self, bond_idx, bead_indices, spatial, force):
        self.idx = bond_idx
        self.bead_indices = bead_indices
        self.spatial = spatial
        self.force = force

    def update_idx(self, new_idx):
        self.idx = new_idx

    def __eq__(self, other):
        if other != None:
            return (self.spatial, self.force) == (other.spatial, other.force)
        else:
            return False

    def __hash__(self):
        return hash((self.spatial, self.force))

    def find_existing_bond(self, bond_list):
        for bond in bond_list:
            if bond == self:
                return bond
        return None


def count_unique_bonds(bond_list):
    unique_bonds = set(bond_list)
    return len(unique_bonds)


class Angle:
    def __init__(self, bond_idx, bead_indices, spatial, force):
        self.idx = bond_idx
        self.bead_indices = bead_indices
        self.spatial = spatial
        self.force = force

    def update_idx(self, new_idx):
        self.idx = new_idx

    def __eq__(self, other):
        if other != None:
            return (self.spatial, self.force) == (other.spatial, other.force)
        else:
            return False

    def __hash__(self):
        return hash((self.spatial, self.force))

    def find_existing_angle(self, angle_list):
        for angle in angle_list:
            if angle == self:
                return angle
        return None


def count_unique_angles(angle_list):
    unique_angles = set(angle_list)
    return len(unique_angles)

def quaternion_multiply(quaternion1, quaternion0):
    #copid from https://stackoverflow.com/questions/39000758/how-to-multiply-two-quaternions-by-python-or-numpy
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)

class Dihedral:
    def __init__(self, bond_idx, bead_indices, k1, k2, k3, k4):
        self.idx = bond_idx
        self.bead_indices = bead_indices  # dihedrals have 4 beads
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.k4 = k4

    def update_idx(self, new_idx):
        self.idx = new_idx

    def __eq__(self, other):
        if other != None:
            return (self.k1, self.k2, self.k3, self.k4) == (
                other.k1,
                other.k2,
                other.k3,
                other.k4,
            )
        else:
            return False

    def __hash__(self):
        return hash((self.k1, self.k2, self.k3, self.k4))

    def find_existing_dihedral(self, dihedral_list):
        for dihedral in dihedral_list:
            if dihedral == self:
                return dihedral
        return None

class ImproperDihedral:
    def __init__(self, bond_idx, bead_indices, chi0,k ):
        self.idx = bond_idx
        self.bead_indices = bead_indices  # dihedrals have 4 beads
        self.k = k
        self.chi0 = chi0


    def update_idx(self, new_idx):
        self.idx = new_idx

    def __eq__(self, other):
        if other != None:
            return (self.k,self.chi0) == (
                other.k,
                other.chi0)
        else:
            return False

    def __hash__(self):
        return hash((self.k,self.chi0))

    def find_existing_improper_dihedral(self, dihedral_list):
        for dihedral in dihedral_list:
            if dihedral == self:
                return dihedral
        return None

def count_unique_dihedrals(dihedral_list):
    unique_dihedrals = set(dihedral_list)
    return len(unique_dihedrals)

def count_unique_improper_dihedrals(dihedral_list):
    unique_dihedrals = set(dihedral_list)
    return len(unique_dihedrals)


class Contents:
    def __init__(self):
        self.contents = []
        self.bead_types = []
        self.bond_types = (
            {}
        )  # dict of dicts. molecule name, bond index, tuple of information
        self.angle_types = {}
        self.dihedral_types = {}
        self.improper_dihedral_types = {}
        self.current_index = [0, 0, 0, 0,0,0] # 0 = bead index, 1 = bond index, 2 = angle index, 3 = dihedral index, 4 = improper index, 5 = body index

    def add_molecule(self, molecule):
        for bead in molecule.types:
            if bead not in self.bead_types:
                self.bead_types.append(bead)
        molecule_len = len(molecule.id)
        if self.bond_types.get(molecule.name) == None:
            if len(molecule.bonds) > 0:
                self.bond_types.update({molecule.name: molecule.bonds[0].idx})
                self.current_index[1] = self.current_index[1] + count_unique_bonds(
                    molecule.bonds
                )
        if self.angle_types.get(molecule.name) == None:
            if len(molecule.angles) > 0:
                self.angle_types.update({molecule.name: molecule.angles[0].idx})
                self.current_index[2] = self.current_index[2] + count_unique_angles(
                    molecule.angles
                )
        if self.dihedral_types.get(molecule.name) == None:
            if len(molecule.dihedrals) > 0:
                self.dihedral_types.update({molecule.name: molecule.dihedrals[0].idx})
                self.current_index[3] = self.current_index[3] + count_unique_dihedrals(
                    molecule.dihedrals
                )
        if self.improper_dihedral_types.get(molecule.name) == None:
            if len(molecule.improper_dihedrals) > 0:
                self.improper_dihedral_types.update({molecule.name: molecule.improper_dihedrals[0].idx})
                self.current_index[4] = self.current_index[4] + count_unique_improper_dihedrals(
                    molecule.improper_dihedrals
                )
        self.contents.append(molecule)
        self.current_index[0] = self.current_index[0] + molecule_len
        if molecule.name == "cholesterol":
            self.current_index[5]+=1


class Molecule:
    def __init__(self, name, molecule_beads, molecule_bonds, contents,rigid_bodies = None):
        if molecule_beads == "W":
            self.id = [contents.current_index[0]]
            self.types = ["W"]
            self.position = [[0.0, 0.0, 0.0]]
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.improper_dihedrals = []
            self.charges = [None]
            self.name = name
            self.body = [-1]
            self.orientation = [(1,0,0,0)]
            self.moment_inertia = [[0,0,0]]
        elif molecule_beads == "Q5-":
            self.id = [contents.current_index[0]]
            self.types = ["Q5"]
            self.position = [[0.0, 0.0, 0.0]]
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.improper_dihedrals = []
            self.charges = [-1]
            self.name = name
            self.body = [-1]
            self.orientation = [(1,0,0,0)]
            self.moment_inertia = [[0,0,0]]
        elif molecule_beads == "Q5+":
            self.id = [contents.current_index[0]]
            self.types = ["Q5"]
            self.position = [[0.0, 0.0, 0.0]]
            self.bonds = []
            self.angles = []
            self.dihedrals = []            
            self.improper_dihedrals = []
            self.charges = [1]
            self.name = name
            self.body = [-1]
            self.orientation = [(1,0,0,0)]
            self.moment_inertia = [[0,0,0]]
        else:
            
            idx, bead_types, positions, charges = extract_beads(
                molecule_beads, contents.current_index
            )
            bonds, angles, dihedrals,improper_dihedrals = extract_bonds(name, molecule_bonds, contents)
            ##TODO: Think about generalizing this to any body though I don't forsee it being needed
            if name == "cholesterol":
                self.id =  idx
                self.types = bead_types
                self.position = positions
                self.bonds = bonds
                self.angles = angles
                self.dihedrals = dihedrals
                self.improper_dihedrals = improper_dihedrals
                self.charges = [None]*(len(idx))
                self.name = name
                self.body = [int((idx[0]))-1]*((len(idx))-1) + [-1]
                
                
                self.orientation = [[1,0,0,0]]*(len(idx))
                self.moment_inertia = [[63.07779361300712, 67.88855195275258, 10.147523004966846]]+ ([[0,0,0]]*((len(idx)-1)))
            else:
                self.id = idx
                self.types = bead_types
                self.position = positions
                self.bonds = bonds
                self.angles = angles
                self.dihedrals = dihedrals
                self.improper_dihedrals = improper_dihedrals
                self.charges = charges
                self.name = name
                self.body = [-1]*len(charges)
                self.orientation = [[1,0,0,0]]*len(charges)
                self.moment_inertia = [[0,0,0]]*len(charges)

    def shift_positions(self, x_shift, y_shift, z_shift):
        for position in self.position:
            position[0] = position[0] + x_shift
            position[1] = position[1] + y_shift
            position[2] = position[2] + z_shift

    def invert_positions(self, is_inverted):
        if is_inverted == True:
            for position in self.position:
                position[0] = position[0]
                position[1] = -position[1]
                position[2] = -position[2]
            orientations = [] 
            for orientation in self.orientation:
                orientations.append([0,1,0,0])
            self.orientation = orientations    

    def rotate(self, theta):
        for position in self.position:
            x = position[0] * np.cos(theta) - position[1] * np.sin(theta)
            y = position[0] * np.sin(theta) + position[1] * np.cos(theta)
            position[0] = x
            position[1] = y
            position[2] = position[2]

    def expand(self, expansion):
        for position in self.position:
            position[0] = position[0] * expansion
            position[1] = position[1] * expansion
            position[2] = position[2] * expansion


def extract_beads(molecule_beads, current_index):
    """
    Read a csv with information about beads present in a molecule and
    translate it into information python can read

    Args:
      molecule_beads (string): string detailing csv with molecule information.
      current_index (list of Int): current_index[0] = bead index. current_index[1] = bond index. current_index[2] = angle_index


    Returns:
      idx (list of ints): bead index
      bead_types (list of strings): bead types (i.e. 'Q1', 'W')
      positions (tuple of (x,y,z)): tuple containing relative x,y,z positions for the initialization of a molecule

    """
    idx = []
    bead_types = []
    positions = []
    charges = []
    with open(molecule_beads, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] != "Index":
                idx.append(str(int(row[0]) + current_index[0]))
                bead_types.append(row[1])
                positions.append([float(row[2]), float(row[3]), float(row[4])])
                if "Q" in row[1]:
                    charges.append(row[5])
                else:
                    charges.append(None)
    return idx, bead_types, positions, charges


def extract_bonds(name, molecule_bonds, contents):
    """
    Read a csv with information about beads present in a molecule and
    translate it into information python can read

    Args:
      name (string): name of molecule
      molecule_bonds (string): string detailing csv with molecule information.
      current_index (list of Int): current_index[0] = bead index. current_index[1] = bond index. current_index[2] = angle_index

    Returns:
      bonds (list of Bond): bonds present in molecule (indexed from 0 up)
      angles (list of Angle): bonds present in molecule (indexed from 0 up)
    """
    current_index = contents.current_index

    # if the molecule doesnt exist, make it from scratch
    if contents.bond_types.get(name) == None:
        bond_idx = current_index[1]
        angle_idx = current_index[2]
        dihedral_idx = current_index[3]
        improper_dihedral_idx = current_index[4]
        with open(molecule_bonds, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                #Based on which elements are present in the bonds file, instantiate that thing
                if row[0] == "i":
                    bonds = []
                    angles = []
                    dihedrals = []
                    improper_dihedrals = []
                elif row[2] == "":
                    # what I want to do is something along the lines of if the bond already
                    # exists, everything stays the same except bond_idx which is not \
                    # incremented and instead the bond appended is the one that already exists
                    test_bond = Bond(
                        bond_idx,
                        [
                            str(int(row[0]) + current_index[0]),
                            str(int(row[1]) + current_index[0]),
                        ],
                        row[4],
                        row[5],
                    )
                    existing_bond = test_bond.find_existing_bond(bonds)
                    if existing_bond == None:
                        bonds.append(test_bond)
                        bond_idx = bond_idx + 1
                    else:
                        #Even if molecule doesn't exist, the bond might not
                        test_bond.update_idx(existing_bond.idx)
                        bonds.append(test_bond)

                elif row[3] == "":
                    test_angle = Angle(
                        angle_idx,
                        [
                            str(int(row[0]) + current_index[0]),
                            str(int(row[1]) + current_index[0]),
                            str(int(row[2]) + current_index[0]),
                        ],
                        row[4],
                        row[5],
                    )
                    existing_angle = test_angle.find_existing_angle(angles)
                    if existing_angle == None:
                        angles.append(test_angle)
                        angle_idx = angle_idx + 1
                    else:
                        test_angle.update_idx(existing_angle.idx)
                        angles.append(test_angle)
                else:
                    if row[6]=="":
                        test_improper_dihedral = ImproperDihedral(
                            improper_dihedral_idx,
                            [
                                str(int(row[0]) + current_index[0]),
                                str(int(row[1]) + current_index[0]),
                                str(int(row[2]) + current_index[0]),
                                str(int(row[3]) + current_index[0]),
                            ],
                            row[4],
                            row[5],
                        )
                        
                        existing_improper_dihedral = test_improper_dihedral.find_existing_improper_dihedral(improper_dihedrals)
                        if existing_improper_dihedral == None:
                            improper_dihedrals.append(test_improper_dihedral)
                            improper_dihedral_idx = improper_dihedral_idx + 1
                        else:
                            test_improper_dihedral.update_idx(existing_improper_dihedral.idx)
                            improper_dihedrals.append(test_improper_dihedral)
                    else:
                        test_dihedral = Dihedral(
                            dihedral_idx,
                            [
                                str(int(row[0]) + current_index[0]),
                                str(int(row[1]) + current_index[0]),
                                str(int(row[2]) + current_index[0]),
                                str(int(row[3]) + current_index[0]),
                            ],
                            row[4],
                            row[5],
                            row[6],
                            row[7],
                        )
                        existing_dihedral = test_dihedral.find_existing_dihedral(dihedrals)
                        if existing_dihedral == None:
                            dihedrals.append(test_dihedral)
                            dihedral_idx = dihedral_idx + 1
                        else:
                            test_dihedral.update_idx(existing_dihedral.idx)
                            dihedrals.append(test_dihedral)

    # otherwise, figure out what bond_idx and angle indexes this molecule started with and fill the molecule w them
    else:
        bond_idx = contents.bond_types.get(name)
        angle_idx = contents.angle_types.get(name)
        dihedral_idx = contents.dihedral_types.get(name)
        improper_dihedral_idx = contents.improper_dihedral_types.get(name)
        with open(molecule_bonds, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == "i":
                    bonds = []
                    angles = []
                    dihedrals = []
                    improper_dihedrals = []
                elif row[2] == "":
                    test_bond = Bond(
                        bond_idx,
                        [
                            str(int(row[0]) + current_index[0]),
                            str(int(row[1]) + current_index[0]),
                        ],
                        row[4],
                        row[5],
                    )
                    existing_bond = test_bond.find_existing_bond(bonds) #I wonder if this is redundant
                    if existing_bond == None:
                        bonds.append(test_bond)
                        bond_idx = bond_idx + 1
                    else:
                        test_bond.update_idx(existing_bond.idx)
                        bonds.append(test_bond)
                elif row[3] == "":
                    test_angle = Angle(
                        angle_idx,
                        [
                            str(int(row[0]) + current_index[0]),
                            str(int(row[1]) + current_index[0]),
                            str(int(row[2]) + current_index[0]),
                        ],
                        row[4],
                        row[5],
                    )
                    existing_angle = test_angle.find_existing_angle(angles)
                    if existing_angle == None:
                        angles.append(test_angle)
                        angle_idx = angle_idx + 1
                    else:
                        test_angle.update_idx(existing_angle.idx)
                        angles.append(test_angle)
                elif row[6]=="":
                    if improper_dihedral_idx != None:
                        test_improper_dihedral = ImproperDihedral(
                            improper_dihedral_idx,
                            [
                                str(int(row[0]) + current_index[0]),
                                str(int(row[1]) + current_index[0]),
                                str(int(row[2]) + current_index[0]),
                                str(int(row[3]) + current_index[0]),
                            ],
                            row[4],
                            row[5],
                        )
                        existing_improper_dihedral = test_improper_dihedral.find_existing_improper_dihedral(improper_dihedrals)
                        if existing_improper_dihedral == None:
                            improper_dihedrals.append(test_improper_dihedral)
                            improper_dihedral_idx = improper_dihedral_idx + 1
                        else:
                            test_improper_dihedral.update_idx(existing_improper_dihedral.idx)
                            improper_dihedrals.append(test_improper_dihedral)

                elif dihedral_idx != None:
                    test_dihedral = Dihedral(
                        dihedral_idx,
                        [
                            str(int(row[0]) + current_index[0]),
                            str(int(row[1]) + current_index[0]),
                            str(int(row[2]) + current_index[0]),
                            str(int(row[3]) + current_index[0]),
                        ],
                        row[4],
                        row[5],
                        row[6],
                        row[7],
                    )
                    existing_dihedral = test_dihedral.find_existing_dihedral(dihedrals)
                    if existing_dihedral == None:
                        dihedrals.append(test_dihedral)
                        dihedral_idx = dihedral_idx + 1
                    else:
                        test_dihedral.update_idx(existing_dihedral.idx)
                        dihedrals.append(test_dihedral)
    return bonds, angles, dihedrals,improper_dihedrals


def path_to_beads(name):
    beads = root + "/martini3/molecules/" + name + "_bead.csv"
    bonds = root + "/martini3/molecules/" + name + "_bonds.csv"
    return beads, bonds

#Instead of writing pydocs for all the add molecules and make molecules I will write a big comment.

#add_MOLECULE makes a molecule and adds it to contents for you. 
#make_MOLECULE makes a molecule and returns said molecule.
    #THis is useful if you want to see if a molecule will fit in a given place before adding it

def add_DOPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DOPC")
    DOPC = Molecule("DOPC", bead_path, bond_path, contents)
    DOPC.expand(expansion)
    DOPC.rotate(theta)
    DOPC.shift_positions(x_shift, y_shift, z_shift)
    DOPC.invert_positions(is_inverted)
    contents.add_molecule(DOPC)
    return contents


def make_DOPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DOPC")
    DOPC = Molecule("DOPC", bead_path, bond_path, contents)
    DOPC.expand(expansion)
    DOPC.rotate(theta)
    DOPC.shift_positions(x_shift, y_shift, z_shift)
    DOPC.invert_positions(is_inverted)
    return DOPC


def add_DPPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DPPC")
    DPPC = Molecule("DPPC", bead_path, bond_path, contents)
    DPPC.expand(expansion)
    DPPC.rotate(theta)
    DPPC.shift_positions(x_shift, y_shift, z_shift)
    DPPC.invert_positions(is_inverted)
    contents.add_molecule(DPPC)

    return contents


def make_DPPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DPPC")
    DPPC = Molecule("DPPC", bead_path, bond_path, contents)
    DPPC.expand(expansion)
    DPPC.rotate(theta)
    DPPC.shift_positions(x_shift, y_shift, z_shift)
    DPPC.invert_positions(is_inverted)
    return DPPC


def add_DOPE(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DOPE")
    DOPE = Molecule("DOPE", bead_path, bond_path, contents)
    DOPE.expand(expansion)
    DOPE.rotate(theta)
    DOPE.shift_positions(x_shift, y_shift, z_shift)
    DOPE.invert_positions(is_inverted)
    contents.add_molecule(DOPE)
    return contents


def make_DOPE(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DOPE")
    DOPE = Molecule("DOPE", bead_path, bond_path, contents)
    DOPE.expand(expansion)
    DOPE.rotate(theta)
    DOPE.shift_positions(x_shift, y_shift, z_shift)
    DOPE.invert_positions(is_inverted)
    return DOPE

def add_DLPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DLPC")
    DLPC = Molecule("DLPC", bead_path, bond_path, contents)
    DLPC.expand(expansion)
    DLPC.rotate(theta)
    DLPC.shift_positions(x_shift, y_shift, z_shift)
    DLPC.invert_positions(is_inverted)
    contents.add_molecule(DLPC)
    return contents


def make_DLPC(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DLPC")
    DLPC = Molecule("DLPC", bead_path, bond_path, contents)
    DLPC.expand(expansion)
    DLPC.rotate(theta)
    DLPC.shift_positions(x_shift, y_shift, z_shift)
    DLPC.invert_positions(is_inverted)
    return DLPC


def add_DPPE(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DPPE")
    DPPE = Molecule("DPPE", bead_path, bond_path, contents)
    DPPE.expand(expansion)
    DPPE.rotate(theta)
    DPPE.shift_positions(x_shift, y_shift, z_shift)
    DPPE.invert_positions(is_inverted)
    contents.add_molecule(DPPE)

    return contents


def make_DPPE(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("DPPE")
    DPPE = Molecule("DPPE", bead_path, bond_path, contents)
    DPPE.expand(expansion)
    DPPE.rotate(theta)
    DPPE.shift_positions(x_shift, y_shift, z_shift)
    DPPE.invert_positions(is_inverted)
    return DPPE



def add_water(contents, x_shift, y_shift, z_shift):
    water = Molecule("W", "W", "", contents)
    water.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(water)
    return contents


def add_neg_ion(contents, x_shift, y_shift, z_shift):
    ion = Molecule("Q5-", "Q5-", "", contents)
    ion.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(ion)
    return contents


def add_pos_ion(contents, x_shift, y_shift, z_shift):
    ion = Molecule("Q5+", "Q5+", "", contents)
    ion.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(ion)
    return contents


def add_dioxane(contents, x_shift, y_shift, z_shift):
    bead_path, bond_path = path_to_beads("Dioxane")
    dioxane = Molecule("Dioxane", bead_path, bond_path, contents)
    dioxane.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(dioxane)
    return contents


def add_cyclohexane(contents, x_shift, y_shift, z_shift):
    bead_path, bond_path = path_to_beads("Cyclohexane")
    cyclohexane = Molecule("Cyclohexane", bead_path, bond_path, contents)
    cyclohexane.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(cyclohexane)
    return contents


def add_PEO(contents, repeats, x_shift=0, y_shift=0, z_shift=0):
    name = "PEO" + str(repeats)
    bead_path, bond_path = path_to_beads(name)
    if not os.path.isfile(bead_path):
        polymers.make_polym("PEO", repeats)
    PEO = Molecule("PEO", bead_path, bond_path, contents)
    PEO.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(PEO)
    return contents


def add_PBD(contents, repeats, x_shift=0, y_shift=0, z_shift=0):
    name = "PBD" + str(repeats)
    bead_path, bond_path = path_to_beads(name)
    if not os.path.isfile(bead_path):
        polymers.make_polym("PBD", repeats)
    PBD = Molecule("PBD", bead_path, bond_path, contents)
    PBD.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(PBD)
    return contents


def add_PBDbPEO(
    contents,
    PBD_repeats,
    PEO_repeats,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
):
    name_PBDbPEO = "PBD" + str(PBD_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PBDbPEO)
    # add PBD_b_PEO bonds
    # write Molecule function join_polym()
    # shift molecules by user defined amount
    # Add them to contents

    # Make Polymers
    if not os.path.isfile(bead_path):
        polymers.make_block_polym("PBDbPEO", PBD_repeats, PEO_repeats)
    PBDbPEO = Molecule("PBDbPEO", bead_path, bond_path, contents)
    PBDbPEO.shift_positions(x_shift, y_shift, z_shift)
    PBDbPEO.invert_positions(is_inverted)
    contents.add_molecule(PBDbPEO)
    return contents


def make_PBDbPEO(
    contents,
    PBD_repeats,
    PEO_repeats,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
    bond_angle=None,
):
    name_PBDbPEO = "PBD" + str(PBD_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PBDbPEO)
    if not os.path.isfile(bead_path):
        polymers.make_block_polym(
            "PBDbPEO", PBD_repeats, PEO_repeats, bond_angle=bond_angle
        )
    PBDbPEO = Molecule("PBDbPEO", bead_path, bond_path, contents)
    PBDbPEO.shift_positions(x_shift, y_shift, z_shift)
    PBDbPEO.invert_positions(is_inverted)
    return PBDbPEO

def make_seq_def(
    contents,
    name_a,
    name_b,
    seq,
    id,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
    bond_angle=None,
):
    name_polym= "seq_" + name_a + "_"+ name_b + "_" + str(id)
    bead_path, bond_path = path_to_beads(name_polym)
    if not os.path.isfile(bead_path):
        polymers.make_sequence(
            name_a,name_b,seq,id,bond_angle=bond_angle
        )
        
    polym = Molecule(seq, bead_path, bond_path, contents)
    polym.shift_positions(x_shift, y_shift, z_shift)
    polym.invert_positions(is_inverted)
    return polym


def make_PCLbPEO(
    contents,
    PCL_repeats,
    PEO_repeats,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
    bond_angle=None,
):
    name_PBDbPEO = "PCL" + str(PCL_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PBDbPEO)
    if not os.path.isfile(bead_path):
        polymers.make_block_polym(
            "PCLbPEO", PCL_repeats, PEO_repeats, bond_angle=bond_angle
        )
    PCLbPEO = Molecule("PCLbPEO", bead_path, bond_path, contents)
    PCLbPEO.shift_positions(x_shift, y_shift, z_shift)
    PCLbPEO.invert_positions(is_inverted)
    return PCLbPEO

def make_PCL1bPEO(
    contents,
    PCL_repeats,
    PEO_repeats,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
    bond_angle=None,
):
    name_PBDbPEO = "PCL1" + str(PCL_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PBDbPEO)
    if not os.path.isfile(bead_path):
        polymers.make_block_polym(
            "PCL1bPEO", PCL_repeats, PEO_repeats, bond_angle=bond_angle
        )
    PCLbPEO = Molecule("PCL1bPEO", bead_path, bond_path, contents)
    PCLbPEO.shift_positions(x_shift, y_shift, z_shift)
    PCLbPEO.invert_positions(is_inverted)
    return PCLbPEO

def make_PCL2bPEO(
    contents,
    PCL_repeats,
    PEO_repeats,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
    bond_angle=None,
):
    name_PBDbPEO = "PCL2" + str(PCL_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PBDbPEO)
    if not os.path.isfile(bead_path):
        polymers.make_block_polym(
            "PCL2bPEO", PCL_repeats, PEO_repeats, bond_angle=bond_angle
        )
    PCLbPEO = Molecule("PCL2bPEO", bead_path, bond_path, contents)
    PCLbPEO.shift_positions(x_shift, y_shift, z_shift)
    PCLbPEO.invert_positions(is_inverted)
    return PCLbPEO

def add_PDMAEMA(contents, repeats, charge_frac=50, x_shift=0, y_shift=0, z_shift=0):
    name = "PDMAEMA" + str(repeats)
    bead_path, bond_path = path_to_beads(name)
    if os.path.isfile(bond_path):
        os.remove(bond_path)
        os.remove(bead_path)

    polymers.make_polym("PDMAEMA", repeats)
    bead_path, num_charged = polymers.add_charges(
        "PDMAEMA" + str(repeats), charge_frac, "SN1", "SQ1", 1
    )
    PDMAEMA = Molecule("PDMAEMA", bead_path, bond_path, contents)
    PDMAEMA.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(PDMAEMA)
    return contents, num_charged


def add_PDMAEMAbPEO(
    contents,
    PDMAEMA_repeats,
    PEO_repeats,
    charge_frac=50,
    x_shift=0,
    y_shift=0,
    z_shift=0,
    is_inverted=False,
):
    name_PDMAEMAbPEO = "PDMAEMA" + str(PDMAEMA_repeats) + "PEO" + str(PEO_repeats)
    bead_path, bond_path = path_to_beads(name_PDMAEMAbPEO)
    if not os.path.isfile(bead_path):
        polymers.make_block_polym("PDMAEMAbPEO", PDMAEMA_repeats, PEO_repeats)
    bead_path, num_charged = polymers.add_charges(
        name_PDMAEMAbPEO, charge_frac, "SN1", "Q1", 1
    )

    PDMAEMAbPEO = Molecule("PDMAEMAbPEO", bead_path, bond_path, contents)

    PDMAEMAbPEO.shift_positions(x_shift, y_shift, z_shift)
    PDMAEMAbPEO.invert_positions(is_inverted)
    contents.add_molecule(PDMAEMAbPEO)
    return contents, num_charged


def add_slab(
    contents,
    name,
    repeats,
    box_size,
    charge_frac=0,
    amorphous=False,
    dipole=False,
    x_shift=0,
    y_shift=0,
    z_shift=0,
):
    nameID = name + str(repeats)
    bead_path, bond_path = path_to_beads(nameID)
    csv_names, num_charged = polymers.make_slab(
        name,
        repeats,
        box_size,
        charge_frac=charge_frac,
        amorphous=amorphous,
        dipole=dipole,
    )
    slab = Molecule(name + "Slab", bead_path, bond_path, contents)
    slab.shift_positions(x_shift, y_shift, z_shift)
    contents.add_molecule(slab)
    return contents, num_charged


def add_PTX(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("PTX")
    PTX = Molecule("PTX", bead_path, bond_path, contents)
    PTX.expand(expansion)
    PTX.rotate(theta)
    PTX.shift_positions(x_shift, y_shift, z_shift)
    PTX.invert_positions(is_inverted)
    contents.add_molecule(PTX)
    return contents

def make_PTX(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False, expansion=1, theta=0
):
    bead_path, bond_path = path_to_beads("PTX")
    PTX = Molecule("PTX", bead_path, bond_path, contents)
    PTX.expand(expansion)
    PTX.rotate(theta)
    PTX.shift_positions(x_shift, y_shift, z_shift)
    PTX.invert_positions(is_inverted)
    return PTX
def add_cholesterol(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False
):
    bead_path, bond_path = path_to_beads("cholesterol")
    cholesterol = Molecule("cholesterol", bead_path, bond_path, contents)
    cholesterol.shift_positions(x_shift, y_shift, z_shift)
    cholesterol.invert_positions(is_inverted)

    contents.add_molecule(cholesterol)

    return contents

def make_cholesterol(
    contents, x_shift=0, y_shift=0, z_shift=0, is_inverted=False
):
    bead_path, bond_path = path_to_beads("cholesterol")
    cholesterol = Molecule("cholesterol", bead_path, bond_path, contents)
    cholesterol.shift_positions(x_shift, y_shift, z_shift)
    cholesterol.invert_positions(is_inverted)


    return cholesterol