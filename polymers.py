# Goal: define a generic function which takes arg n (number of repeats), ref to polymer file, and returns a Molecule of length n centered at the origin
import sys

sys.path.append("../")
import csv
import os
import numpy as np

script_path = os.path.abspath(__file__)
root = script_path.split("/martini3")[0]


def count_bonds(reader):
    """
    Counts number of bonds in an _bonds.csv

    Args:
      reader (list): open file _bonds.csv

    Returns:
      num_bonds (int):
    """
    num_bonds = 0
    for row in reader:
        if row[2] == "":
            num_bonds = num_bonds + 1
    return num_bonds


def count_angles(reader):
    """
    Counts number of angles in an _bonds.csv

    Args:
      reader (list): open file _bonds.csv

    Returns:
      num_angles (int):
    """
    num_angles = 0
    for row in reader:
        if row[2] != "" and row[3] == "":
            num_angles = num_angles + 1
    return num_angles


def count_dihedrals(reader):
    """
    Counts number of dihedrals in an _bonds.csv

    Args:
      reader (list): open file _bonds.csv

    Returns:
      num_dihedrals (int):
    """
    num_dihedrals = 0
    for row in reader:
        try:
            if row[6] != "" and row[0] != "i":
                num_dihedrals = num_dihedrals + 1
        except:
            1 + 1
    return num_dihedrals


def count_beads(reader):
    """
    Counts number of beads in an _bead.csv

    Args:
      reader (list): open file _bead.csv

    Returns:
      num_dihedrals (int):
    """
    num_beads = 0
    for row in reader:
        if row[0] != "Index":
            num_beads = num_beads + 1
    return num_beads


def make_bonds(name, csv_name_bonds, repeats, block=(False, None), a_block_repeats=0):
    """
    Creates (or appends to) file (named csv_name_bonds) in martini3/molecules/ 
    which contains a polymer with N repeats and assigns angles/bonds according 
    to the name (name) of the molecule which is being added. If adding to a BCP,
    indexing starts after a_block_repeats. Only tested for PCL, PEO, PBD. Previous 
    version worked for PDMAEMA but untested.

    Args:
      name (string): name of polymer added (PEO,PCL,PBD)
      csv_name_bonds (string): name of file to save to
      repeats (int): N
      block (tuple): (bool: is_blocked, string: position ("first","middle", "second"))
      a_block_repeats (int): If blocked, number of repeats in the a-block

    Returns:
      bondlengths: list of all unique bond lengths in str(name)_bonds.csv
      bondangles: list of all unique bond angles str(name)_bonds.csv
    """
    blocked = block[0]
    block_pos = block[1]
    bondlengths = []
    bondangles = []
    if repeats ==1 and "PCL" not in name:
        if (block_pos != "first" and block_pos != "second")  or name !="PEO":
            if block_pos == "first":
                with open(csv_name_bonds, "a", newline="") as write_file:
                    writer = csv.writer(write_file)
                    with open(
                        root + "/martini3/molecules/" + name + "_bonds.csv",
                        newline="",
                    ) as f:
                        reader = csv.reader(f)
                        reader = list(reader)
                        writer.writerow(reader[0])
            return bondlengths, bondangles
    with open(csv_name_bonds, "a", newline="") as write_file:
        writer = csv.writer(write_file)
        with open(
            root + "/martini3/molecules/" + name + "_bonds.csv",
            newline="",
        ) as f:
            reader = csv.reader(f)
            reader = list(reader)
            num_bonds = count_bonds(reader)
            num_angles = count_angles(reader)
            num_dihedrals = count_dihedrals(reader)
            if name == "PEO":
                end_beads = 1  # if more complex polymers with end beads that dont add just 1 bond and 1 angle, will have to change this
            elif name == "PBD" or name == "PDMAEMA" or "PCL" in name:
                end_beads = 0
            bonds_per_repeat = num_bonds - end_beads
            angle_per_repeat = num_angles - end_beads
            dihedrals_per_repeat = num_dihedrals

            # write header
            if block_pos == "first" or blocked == False:
                writer.writerow(reader[0])

            # initiate bondangles and bondlengths bonds
            for a in range(num_bonds):
                bondlengths.append(float(reader[a + 1][4]))
            for a in range(num_angles):
                bondangles.append(float(reader[a + num_bonds + 1][4]))

            # write first end beads bond
            if end_beads > 0 and (blocked == False or block_pos == "first"):
                writer.writerow(reader[1])
            monomer = 0
            bond = 0
            i = monomer * bonds_per_repeat
            j = monomer * bonds_per_repeat
            for monomer in range(repeats - 1):
                for bond in range(bonds_per_repeat):
                    i = monomer * bonds_per_repeat
                    j = monomer * bonds_per_repeat
                    writer.writerow(
                        [
                            str(
                                i
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][0])
                            ),
                            str(
                                j
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][1])
                            ),
                            "",
                            "",
                            reader[end_beads + bond + 1][4],
                            reader[end_beads + bond + 1][5],
                        ]
                    )
            # finish last molecule
            if bonds_per_repeat > 1 and repeats >1:
                
                for bond in range(bonds_per_repeat - 1):
                    writer.writerow(
                        [
                            str(
                                i
                                + bonds_per_repeat
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][0])
                            ),
                            str(
                                j
                                + bonds_per_repeat
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][1])
                            ),
                            "",
                            "",
                            reader[end_beads + bond + 1][4],
                            reader[end_beads + bond + 1][5],
                        ]
                    )
            elif bonds_per_repeat>1:
                writer.writerow(
                        [
                            str(
                                bonds_per_repeat -2
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][0])
                            ),
                            str(
                                bonds_per_repeat -2
                                + a_block_repeats
                                + int(reader[end_beads + bond + 1][1])
                            ),
                            "",
                            "",
                            reader[end_beads + bond + 1][4],
                            reader[end_beads + bond + 1][5],
                        ]
                    )
            if end_beads > 0 and (blocked == False or block_pos == "second"):
                i = monomer * bonds_per_repeat
                j = monomer * bonds_per_repeat
                writer.writerow(
                    [
                        str(
                            i
                            + bonds_per_repeat
                            + a_block_repeats
                            + int(reader[end_beads + bond + 1][0])
                        ),
                        str(
                            j
                            + bonds_per_repeat
                            + a_block_repeats
                            + int(reader[end_beads + bond + 1][1])
                        ),
                        "",
                        "",
                        reader[1][4],
                        reader[1][5],
                    ]
                )

            # write angles
            if end_beads > 0 and (blocked == False or block_pos == "first")and repeats >1:
                writer.writerow(reader[end_beads + 2])
            
            for monomer in range((repeats - 2)):
                for angle in range(angle_per_repeat):
                    i = monomer * angle_per_repeat
                    j = monomer * angle_per_repeat
                    k = monomer * angle_per_repeat
                    writer.writerow(
                        [
                            str(
                                i
                                + a_block_repeats
                                + int(reader[end_beads + 1 + num_bonds + angle][0])
                            ),
                            str(
                                j
                                + a_block_repeats
                                + int(reader[end_beads + 1 + num_bonds + angle][1])
                            ),
                            str(
                                k
                                + a_block_repeats
                                + int(reader[end_beads + 1 + num_bonds + angle][2])
                            ),
                            "",
                            reader[end_beads + num_bonds + angle + 1][4],
                            reader[end_beads + num_bonds + angle + 1][5],
                        ]
                    )
            adj = angle_per_repeat if repeats-2>0 else 0
            if angle_per_repeat > 1 and repeats>1:
                for angle in range(angle_per_repeat):
                    writer.writerow(
                        [
                            str(
                                i
                                + a_block_repeats+angle+adj
                                + int(reader[end_beads + num_bonds + 1][0])
                            ),
                            str(
                                i
                                + a_block_repeats+angle+adj
                                + int(reader[end_beads + num_bonds + 1][1])
                            ),
                            str(
                                i + a_block_repeats+angle+adj
                                + int(reader[end_beads + num_bonds + 1][2])
                            ),
                            "",
                            reader[num_bonds + 1+angle][4],
                            reader[num_bonds + 1+angle][5],
                        ]
                    )
            if end_beads > 0 and (blocked == False or block_pos == "second")and repeats >1:
                writer.writerow(
                    [
                        str(
                            i+adj
                            + a_block_repeats
                            + int(reader[end_beads + num_bonds + 1][0])
                        ),
                        str(
                            i+adj
                            + a_block_repeats
                            + int(reader[end_beads + num_bonds + 1][1])
                        ),
                        str(
                            i+adj
                            + a_block_repeats
                            + int(reader[end_beads + num_bonds + 1][2])
                        ),
                        "",
                        reader[num_bonds + 1][4],
                        reader[num_bonds + 1][5],
                    ]
                )

            # write dihedrals
            if dihedrals_per_repeat != 0:
                for monomer in range(repeats - 3):
                    i = monomer + 1
                    j = monomer + 2
                    k = monomer + 3
                    l = monomer + 4
                    writer.writerow(
                        [
                            str(i + a_block_repeats),
                            str(j + a_block_repeats),
                            str(k + a_block_repeats),
                            str(l + a_block_repeats),
                            reader[5][4],
                            reader[5][5],
                            reader[5][6],
                            reader[5][7],
                        ]
                    )
    return bondlengths, bondangles


def update_curr_pos(curr_pos, bondlength, bondangle, counter, pitch=0.25):
    """
    Places molecules in a helix. Returns the position of the bead just placed

    Args:
      curr_pos (list of floats): [x,y,z]
      bondlength (float): length of equilibrium bond... (doesn't really matter so i typically right .4)
      bondangle (float): determines the angle between beads (maybe)
      counter (int): determines the number of radians around the circle that have been traveled
      pitch (float): pitch of the helix

    Returns:
      new_curr_pos (list of floats) position of bead just placed
    """
    x = curr_pos[0]
    y = curr_pos[1]
    z = curr_pos[2]
    newx = (
        x
        + bondlength * np.sin(bondangle * counter / 2 * np.pi / 180)
        - pitch / 2 * bondlength * np.sin(counter * 2 * np.pi * pitch)
    )

    newy = (
        y
        + bondlength * np.cos(bondangle * counter / 2 * np.pi / 180)
        - pitch / 2 * bondlength * np.cos(counter * 2 * np.pi * pitch)
    )
    newz = z + bondlength * pitch
    new_curr_pos = [newx, newy, newz]

    return new_curr_pos


def get_bond_angle(bond_angles, index):
    """
    Either return the ith bond angle from the _bonds.csv or return 45

    Args:
      bond_angles (list of ints): list of bond angles
      index (int): index

    Returns:
      angle (int): angle
    """
    # molecyule can exist without an angle but cant without a bond
    try:
        angle = bond_angles[index]
    except:
        angle = 45
    return angle


def make_beads(
    name,
    csv_name_beads,
    bondlengths,
    bondangles,
    repeats,
    block=(False, None),
    a_block_repeats=0,
    fin_pos=[0.0, 0.0, 0.0],
    bond_angle=45,
):
    """
    Creates (or appends to) file (named csv_name_bead) in martini3/molecules/ 
    which contains a polymer with N repeats and assigns beads and bead positions according 
    to the name (name) of the molecule which is being added. If adding to a BCP,
    indexing starts after a_block_repeats. Only tested for PCL, PEO, PBD. Previous 
    version worked for PDMAEMA but untested.

    Args:
      name (string): name of polymer added (PEO,PCL,PBD)
      csv_name_beads (string): name of file to save to
      bondlengths: DEPRECATED
      bondangles: DEPRECATED
      repeats (int): N
      block (tuple): (bool: is_blocked, string: position ("first","middle", "second"))
      a_block_repeats (int): If blocked, number of repeats in the a-block
      fin_pos (list of int): position of last bead placed if continuing a bcp
      bond_angle: determines the angle of the helix bonds are placed in

    Returns:
      None
    """
    blocked = block[0]
    block_pos = block[1]
    with open(csv_name_beads, "a", newline="") as write_file:
        writer = csv.writer(write_file)
        with open(
            root + "/martini3/molecules/" + name + "_bead" + ".csv",
            newline="",
        ) as f:
            reader = csv.reader(f)
            reader = list(reader)
            curr_pos = fin_pos
            num_beads = count_beads(reader)

            # write header
            if not blocked or block_pos == "first":
                writer.writerow(reader[0])

            # write beads
            counter = 1

            # for end groups
            if name == "PEO":
                end_beads = 1
                if blocked == False or block_pos == "first":
                    writer.writerow(reader[1])
                    bond_angle_i = 45
                    curr_pos = update_curr_pos(
                        curr_pos, .4, bond_angle_i, counter
                    )
                    counter = counter + 1
            elif name == "PBD" or name == "PDMAEMA" or "PCL" in name:
                end_beads = 0
            beads_per_repeat = num_beads - end_beads

            #Handle repeat unit of 1 bead/1 bead w side chains
            if name == "PDMAEMA" or name =="PBD" or name == "PEO":
                for i in range(repeats):
                    side_pos = curr_pos.copy()
                    side_counter = counter
                    for bead in range(beads_per_repeat):
                        writer.writerow(
                            [
                                str(
                                    i * beads_per_repeat
                                    + bead
                                    + 1
                                    + end_beads
                                    + a_block_repeats
                                ),
                                reader[end_beads + 1 + bead][1],
                                str(side_pos[0]),
                                str(side_pos[1]),
                                str(side_pos[2]),
                            ]
                        )
                        side_pos = update_curr_pos(
                            side_pos, 0.4, -135 + 45 * side_counter, counter=1, pitch=0
                        )
                        side_counter = side_counter + 1
                    counter = counter + 1
                    if bond_angle == None:
                        bond_angle = 45

                    curr_pos = update_curr_pos(
                        curr_pos,
                        .4,
                        bond_angle,
                        counter,
                        pitch=0.25,
                    )
            elif "PCL" in name:
                #PCL has 2 beads per repeat unit
                for i in range(2*repeats):
                    side_pos = curr_pos.copy()
                    writer.writerow(
                        [
                            str(i 
                                + 1
                                + end_beads
                                + a_block_repeats
                            ),
                            reader[1+ i % 2][1],
                            str(side_pos[0]),
                            str(side_pos[1]),
                            str(side_pos[2]),
                        ]
                    )
                    counter = counter + 1
                    if bond_angle == None:
                        bond_angle = 45
                    curr_pos = update_curr_pos(
                        curr_pos,
                        .4,
                        bond_angle,
                        counter,
                        pitch=0.25,
                    )
            if end_beads > 0 and (blocked == False or block_pos == "second"):
                writer.writerow(
                    [
                        str(repeats * beads_per_repeat + 2 + a_block_repeats),
                        reader[end_beads][1],
                        str(curr_pos[0]),
                        str(curr_pos[1]),
                        str(curr_pos[2]+.1),
                    ]
                )


def make_polym(name, repeats):
    """
    Makes a homopolymer. Saves information in name+str(repeats).

    Args:
      name (string): name of polymer added (PEO,PCL,PBD)
      repeats (int): N

    Returns:
      csv_names (tuple): (name of _bead.csv file, name of _bonds.csv file)
    """
    # If blocked = true, the head group of this polymer will not be included on the starting end
    block = (False, None)
    csv_names = (
        root + "/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv",
        root
        + "/martini3/molecules/"
        + name
        + str(repeats)
        + "_bonds"
        + ".csv",
    )
    bondlengths, bondangles = make_bonds(name, csv_names[1], repeats, block)
    make_beads(name, csv_names[0], bondlengths, bondangles, repeats, block)
    return csv_names


def make_bonds_slab(name, write_name):
    """
    Make an empty bonds file with just a header

    Args:
      name (string): name of file to read from (any real file with a header)
      write_name (string): name of file to save to

    Returns:
      None
    """
    with open(
        root + "/martini3/molecules/" + name + "_bonds" + ".csv",
        newline="",
    ) as f:
        reader = csv.reader(f)
        bonds = list(reader)
    with open(write_name, "w", newline="") as write_file:
        writer = csv.writer(write_file)
        writer.writerow(bonds[0])


def make_beads_slab(
    name, write_name, repeats, box_size, charge_frac=0, amorphous=False, dipole=False
):
    """
    Creates a beads file with a slab of thickness repeats  with beads of space ~.47 apart 
    and potentially an amorphous upper and bottom layer and charg_frac percent of charged beads

    Args:
      name (String): name of initial bead file to read from (i.e. Si_bead)
      write_name (string): name of file to write to 
      repeats (int): number of layers of bead
      box_size (int): x/y dimension of box (assumes tetrahedral)
      charge_frac (int): fraction of surface beads to be replaced with negatively charged surface beads
      amorphous (bool): if true, make the upper layer randomly amorphous
      dipole (bool): if true, make the top and bottom surfaces have opposite charges

    Returns:
      num_charged: number of charged beads placed in slab
    """
    beads_per_layer = int(box_size // 0.47)
    rng = np.random.default_rng()
    with open(
        root + "/martini3/molecules/" + name + "_bead" + ".csv",
        newline="",
    ) as f:
        reader = csv.reader(f)
        bead = list(reader)
    bead_num = 1
    num_charged = 0
    with open(write_name, "w", newline="") as write_file:
        writer = csv.writer(write_file)
        writer.writerow(bead[0])

        for iter_num in range(repeats):
            i = iter_num - repeats // 2
            z = -0.47 * (i)
            bounds = np.linspace(
                -box_size / 2 + 0.24,
                box_size / 2 - 0.24,
                beads_per_layer,
            )
            for j in bounds:
                for k in bounds:
                    if i**2 >= (repeats // 2) ** 2:
                        if amorphous:
                            offset = 0.235 * (rng.random(size=(3)) * 2 - 1)
                        else:
                            offset = np.array([0.0, 0.0, 0.0])

                        if np.random.randint(0, 100) < charge_frac:
                            if dipole == True and i < 0:
                                charge = +1
                                bead_type = "Q5"
                            else:
                                charge = -1
                                bead_type = "Q5"
                                num_charged = num_charged + 1

                        else:
                            charge = 0
                            bead_type = bead[1][1]
                        writer.writerow(
                            [
                                str(bead_num),
                                bead_type,
                                str(j + offset[0]),
                                str(k + offset[1]),
                                str(z + offset[2]),
                                charge,
                            ]
                        )
                    else:
                        writer.writerow(
                            [str(bead_num), bead[1][1], str(j), str(k), str(z)]
                        )
                    bead_num = bead_num + 1

    return num_charged


def make_slab(name, repeats, box_size, charge_frac=0, amorphous=False, dipole=False):
    """
    Creates a beads and bonds file with a slab of thickness repeats  with beads of space ~.47 apart 
    and potentially an amorphous upper and bottom layer and charg_frac percent of charged beads

    Args:
      name (String): name of initial bead file to read from (i.e. Si_bead)
      repeats (int): number of layers of bead
      box_size (int): x/y dimension of box (assumes tetrahedral)
      charge_frac (int): fraction of surface beads to be replaced with negatively charged surface beads
      amorphous (bool): if true, make the upper layer randomly amorphous
      dipole (bool): if true, make the top and bottom surfaces have opposite charges

    Returns:
      num_charged: number of charged beads placed in slab
    """
    csv_names = (
        root + "/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv",
        root
        + "/martini3/molecules/"
        + name
        + str(repeats)
        + "_bonds"
        + ".csv",
    )
    num_charged = make_beads_slab(
        name,
        csv_names[0],
        repeats,
        box_size,
        charge_frac=charge_frac,
        amorphous=amorphous,
        dipole=dipole,
    )
    make_bonds_slab(name, csv_names[1])
    return csv_names, num_charged


def get_finpos(name):
    """
    Get the position of the bottommost bead written to a _bead file. add .48 to the fin_pos_z so no overlap

    Args:
      name (string): name of _bead file to get the positions from.

    Returns:
      fin_pos (list of floats): [x,y,z]
    """
    with open(name) as f:
        reader = csv.reader(f)
        reader = list(reader)
        x_pos, y_pos, z_pos = reader[-1][2:5]
    return [
        np.float64(x_pos),
        np.float64(y_pos),
        np.float64(z_pos) + 0.48,
    ]  # adjust z so nothing is ontop of eachother


def block_bonds(name, file_path, a_block_repeats, bonds_per_repeat=1, end_beads=0,angles_off = False):
    """
    Appends to file (named file_path) in martini3/molecules/ with the bond/angles that connect two blocks of a block-copolymer.
     Only tested for PCL, PEO, PBD. Previous 
    version worked for PDMAEMA but untested.

    Args:
      name (string): name of block_bond (i.e. PCLbPEO,PCLbPBD)
      file_path (string): name of file to save to
      a_block_repeats (int): If blocked, number of repeats in the a-block
      bonds_per_repeat (int): DEPRECATED
      end_beads (int): DEPRECATED
      angles_off (bool): if true, do not attempt to place any angles. Typically used in conjunction with block_angles

    Returns:
      None
    """
    with open(file_path, "a", newline="") as write_file:
        writer = csv.writer(write_file)
        read_path = root + "/martini3/molecules/" + name + "_bonds" + ".csv"
        with open(read_path) as f:
            reader = csv.reader(f)
            reader = list(reader)[1::]
            for row in reader:
                if row[2] == "":
                    writer.writerow(
                        [
                            str(int(row[0]) + a_block_repeats - 1),
                            str(int(row[1]) + a_block_repeats - 1),
                            "",
                            "",
                            row[4],
                            row[5],
                        ]
                    )
                elif row[2] != "":
                    if (angles_off == False):
                        writer.writerow(
                            [
                                str(int(row[0]) + a_block_repeats - 1),
                                str(int(row[1]) + a_block_repeats - 1),
                                str(int(row[2]) + a_block_repeats - 1),
                                "",
                                row[4],
                                row[5],
                            ]
                        )

def block_angles(file_path, a_block_repeats, a2_on = True):
    """
    Appends to file (named file_path) in martini3/molecules/ with the angles that connect two blocks of a block-copolymer.
    Places 2 angles with angle = 135 and angle strength = 120. 

    Args:
      file_path (string): name of file to save to
      a_block_repeats (int): If blocked, number of repeats in the a-block
      a2_on (bool): if false, do not add the second angle (i.e. at the end of a polymer)

    Returns:
      None
    """
    #Note. block_bonds will insert angles if angle_type ==True. This is if you want seperate generic angles
    with open(file_path, "a", newline="") as write_file:
        writer = csv.writer(write_file)
        index1 = a_block_repeats-2
        index2 = a_block_repeats-1
        index3 = a_block_repeats-0
        index4 = a_block_repeats+1
        if index1>-1:
            writer.writerow(
                [
                    index1,
                    index2,
                    index3,
                    "",
                    135,
                    20,
                ]
            )
        if a2_on:
            writer.writerow(
                [
                    index2,
                    index3,
                    index4,
                    "",
                    135,
                    20,
                ]
            )

def make_block_polym(name, a_block_repeats, b_block_repeats, bond_angle=None):
    """
    Creates a block co-polymer _bead and _bonds files.

    Args:
      name (string): name of block co polymer (form "PBDbPEO")
      a_block_repeats (int): number of repeats in the a-block
      b_block_repeats (int): number of repeats in the b-block
      bond_angle (int): angle which determines shape of helix

    Returns:
      csv_names (tuple of strings)
    """
    name_a, name_b = name.split("b")
    parsed_name = name_a + str(a_block_repeats) + name_b + str(b_block_repeats)
    csv_names = (
        root + "/martini3/molecules/" + parsed_name + "_bead" + ".csv",
        root + "/martini3/molecules/" + parsed_name + "_bonds" + ".csv",
    )
    bondlengths, bondangles = make_bonds(
        name_a, csv_names[1], a_block_repeats, block=(True, "first")
    )
    if name_a == "PEO":
        end_beads_a = 1  # if more complex polymers with end beads that dont add just 1 bond and 1 angle, will have to change this
    elif name_a == "PBD" or name_a=="PDMAEMA"or "PCL" in name_a:
        end_beads_a = 0
    else:
        raise Exception("This pol.ymer isnt implemented yet")
    if name_b == "PEO":
        end_beads_b = 1  # if more complex polymers with end beads that dont add just 1 bond and 1 angle, will have to change this
    elif name_b == "PBD" or "PDMAEMA"or "PCL" in name_b:
        end_beads_b = 0
    else:
        raise Exception("This pol.ymer isnt implemented yet")

    bonds_per_repeat_a = len(bondlengths) - end_beads_a
    make_beads(
        name_a,
        csv_names[0],
        bondlengths,
        bondangles,
        a_block_repeats,
        block=(True, "first"),
        bond_angle=bond_angle,
    )
    fin_pos = get_finpos(csv_names[0])
    # a_block_repeats has a -1 to account for the missing head group on PEO
    bondlengths, bondangles = make_bonds(
        name_b,
        csv_names[1],
        b_block_repeats,
        block=(True, "second"),
        a_block_repeats=a_block_repeats * bonds_per_repeat_a - end_beads_b,
    )
    make_beads(
        name_b,
        csv_names[0],
        bondlengths,
        bondangles,
        b_block_repeats,
        block=(True, "second"),
        a_block_repeats=a_block_repeats * bonds_per_repeat_a - end_beads_b,
        fin_pos=fin_pos,
        bond_angle=bond_angle,
    )
    block_bonds(name, csv_names[1], a_block_repeats * bonds_per_repeat_a - end_beads_b)
    return csv_names

def name(a_or_b,name_a,name_b):
    """
    Return name a if a_or_b = A otherwise return B

    Args:
      a_or_b (string): Either "A" or "B"
      name_a (string): name_a
      name_b (string): name_
      
    Returns:
      name_a if "A" else name_b
    """
    return name_a if a_or_b =="A" else name_b

def len_beads(path):
    """
    Number of beads in _bead

    Args:
      path (string): filename/path of _bead
      
    Returns:
      beads_in_csv (int)
    """
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return len(list(reader))-1
    
def beads_per_monomer(name):
    """
    Hard-coded beads per monomer

    Args:
      name (string): name of monomer
      
    Returns:
      beads_per_monomer (int)
    """
    if name =="PEO" or name == "PBD":
        return 1
    elif "PCL" in name:
        return 2
    elif name =="PDMAEMA":
        return 3
    
def make_sequence(name_a,name_b,sequence, id,bond_angle=None):
    """
    Creates a sequence defined macromolecular (PCL,PEO,PBD are tested).

    Args:
      name_a (string): name of monomer a (form "PBD")
      name_b (string): name of monomer b
      sequence (string): sequence (form "AAA-BB-AA-B")
      id (string): id of sequence to be added
      bond_angle (int): angle which determines shape of helix

    Returns:
      csv_names (tuple of strings)
    """
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    csv_names = (
        root + "/martini3/molecules/" + parsed_name + "_bead" + ".csv",
        root + "/martini3/molecules/" + parsed_name + "_bonds" + ".csv",
    )
    path_a_bead = root + "/martini3/molecules/" + name_a + "_bead" + ".csv"
    path_b_bead = root + "/martini3/molecules/" + name_b + "_bead" + ".csv"
    #my initial thought looking at this is we can break any sequence down to N homopolymers. Then. All we have to do is add the linking bonds
    if os.path.exists(root+ "/martini3/molecules/" +name_a + "b" + name_b + "_bonds.csv"):
        name_block = name_a + "b" + name_b
    else:
        name_block = name_b + "b" + name_a
    split_seq = sequence.split("-")      
    
    a_or_b = name(split_seq[0][0],name_a,name_b)
    bondlengths, bondangles = make_bonds(
        a_or_b, csv_names[1], len(split_seq[0]), block=(True, "first")
    )
    make_beads(
        a_or_b,
        csv_names[0],
        bondlengths,
        bondangles,
        len(split_seq[0]),
        block=(True, "first"),
        bond_angle=bond_angle,
    )

    beads_places = len_beads(csv_names[0])
    i = 1

    while len(split_seq[i:]) > 0:
        pos = "middle" if len(split_seq[i:])>1 else "second"
        seq = split_seq[i]
        number_beads_in_csv = len_beads(path_a_bead) if seq[0]=="A" else len_beads(path_b_bead)
        num_beads_per_monomer= beads_per_monomer(name(seq[0],name_a,name_b))
        
        name_s = name(seq[0],name_a,name_b)
        fin_pos = get_finpos(csv_names[0])
        if number_beads_in_csv-num_beads_per_monomer>0 and pos == "second" and len(seq)==1:
            adj = num_beads_per_monomer
        else:
            adj = 0
        # a_block_repeats has a -1 to account for the missing head group on PEO
        if number_beads_in_csv-num_beads_per_monomer>0 and pos == "second" and len(split_seq[i:]) ==1:
            adj_2 = num_beads_per_monomer
        else:
            adj_2 = 0

        block_bonds(name_block, csv_names[1], beads_places-1,angles_off=True)
        a2_on = False if len(split_seq[i:]) ==1 and len(seq)==1 and number_beads_in_csv == 1 else True
        block_angles(csv_names[1],beads_places,a2_on)
        bondlengths, bondangles = make_bonds(
            name_s,
            csv_names[1],
            len(seq),
            block=(True, pos),
            a_block_repeats=beads_places-number_beads_in_csv+num_beads_per_monomer-adj,
        )
        make_beads(
            name_s,
            csv_names[0],
            bondlengths,
            bondangles,
            len(seq),
            block=(True, pos),
            a_block_repeats=beads_places-number_beads_in_csv+num_beads_per_monomer, # this one and many that follow need changed PCL
            fin_pos=fin_pos,
            bond_angle=bond_angle,
        )
        
        beads_places += len(seq)*num_beads_per_monomer
        i+=1
    # block_bonds(name_block, csv_names[1], beads_places-1,angles_off=True)
    return csv_names

def add_charges(name, charge_frac, replace_bead, charged_bead, charge):
    """
    Given an existing polymer, replace charge_frac of replace_bead with charged_bead

    Args:
      name (string): name of block co polymer (form "PBDbPEO")
      charge_frac (int): fraction of beads to replace
      replace_bead (string): bead to replace
      charged bead (string): bead to replace replace_bead with
      charge (int): charge of charged bead 

    Returns:
      write_name (string): csv of new bonds
      num_charged: number of charged beads replaced 
    """
    read_name = root + "/martini3/molecules/" + name + "_bead" + ".csv"

    write_name = (
        root
        + "/martini3/molecules/"
        + name
        + "c"
        + str(charge_frac)
        + "_bead"
        + ".csv"
    )
    with open(read_name, "r", newline="") as read_input, open(
        write_name, "w", newline=""
    ) as write_input:
        reader = csv.reader(read_input)
        writer = csv.writer(write_input)
        num_charged = 0
        for row in reader:
            if replace_bead in row and np.random.randint(0, 101) < charge_frac:
                writer.writerow(
                    [row[0], charged_bead, row[2], row[3], row[4], str(charge)]
                )
                num_charged = num_charged + 1

            else:
                writer.writerow(row)
    return write_name, num_charged
