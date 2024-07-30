import pytest
import sys
import os
import csv
import numpy as np

sys.path.append("../")

from martini3 import polymers


# What do I want to test here since its big
# Now I have enough written that I feel like I can implement PDMAEMA
# Bt its late so whatever
@pytest.fixture()
def root():
    script_path = os.path.abspath(__file__)
    root = script_path.split("/projects")[0]
    return root


@pytest.fixture
def make_PEO_bonds(root):
    repeats = 10
    name = "PEO"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    bond_angles, bond_lengths = polymers.make_bonds(name, path, repeats)

    return bond_angles, bond_lengths


@pytest.fixture
def make_PBD_bonds(root):
    repeats = 10
    name = "PBD"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    bond_angles, bond_lengths = polymers.make_bonds(name, path, repeats)

    return bond_angles, bond_lengths


@pytest.fixture
def make_PBD_beads(root, make_PBD_bonds):
    repeats = 10
    name = "PBD"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    polymers.make_beads(name, path, make_PBD_bonds[0], make_PBD_bonds[1], repeats)


@pytest.fixture
def make_PEO_beads(root, make_PEO_bonds):
    repeats = 10
    name = "PEO"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    polymers.make_beads(name, path, make_PEO_bonds[0], make_PEO_bonds[1], repeats)


@pytest.fixture
def make_PBDbPEO(root):
    name = "PBDbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_block_polym(name,repeats_a,repeats_b)

@pytest.fixture
def make_PCLbPEO(root):
    name = "PCLbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_block_polym(name,repeats_a,repeats_b)

@pytest.fixture
def reader_PBDbPEO_bead(root, make_PBDbPEO):
    name = "PBDbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_PBDbPEO
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_PBDbPEO_bonds(root, make_PBDbPEO):
    name = "PBDbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_PBDbPEO
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_PCLbPEO_bead(root, make_PCLbPEO):
    name = "PCLbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_PCLbPEO
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_PCLbPEO_bonds(root, make_PCLbPEO):
    name = "PCLbPEO"
    name_a, name_b = name.split("b")
    repeats_a = 10
    repeats_b = 10
    parsed_name = name_a + str(repeats_a) + name_b + str(repeats_b)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_PCLbPEO
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


@pytest.fixture
def reader_PEO(root, make_PEO_beads):
    name = "PEO"
    repeats = 10
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    make_PEO_beads
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)

@pytest.fixture
def reader_PEO_bonds(root, make_PEO_bonds):
    name = "PEO"
    repeats = 10
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )
    make_PEO_bonds
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)

@pytest.fixture
def reader_PBD(root, make_PBD_beads):
    name = "PBD"
    repeats = 10
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    make_PBD_beads
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)

@pytest.fixture
def reader_PBD_bonds(root, make_PBD_bonds):
    name = "PBD"
    repeats = 10
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )
    make_PBD_bonds
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
    
@pytest.fixture
def make_PDMAEMA_bonds(root):
    repeats = 10
    name = "PDMAEMA"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    bond_angles, bond_lengths = polymers.make_bonds(name, path, repeats)

    return bond_angles, bond_lengths


@pytest.fixture
def make_PDMAEMA_beads(root, make_PDMAEMA_bonds):
    repeats = 10
    name = "PDMAEMA"
    path = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    if os.path.exists(path):
        os.remove(path)
    polymers.make_beads(
        name, path, make_PDMAEMA_bonds[0], make_PDMAEMA_bonds[1], repeats
    )


def test_make_PEO(root):
    name = "PBD"
    repeats = 10
    path_beads = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bead" + ".csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + name + str(repeats) + "_bonds" + ".csv"
    )

    csv_names = polymers.make_polym(name, repeats)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    if os.path.exists(path_beads):
        os.remove(path_beads)
    assert csv_names == (path_beads, path_bonds)


def test_bond_angles_PEO(make_PEO_bonds):
    assert len(make_PEO_bonds[1]) == 2


def test_bond_lengths_PEO(make_PEO_bonds):
    assert len(make_PEO_bonds[0]) == 2


def test_bond_angles_PBD(make_PBD_bonds):
    assert len(make_PBD_bonds[1]) == 1


def test_bond_lengths_PBD(make_PBD_bonds):
    assert len(make_PBD_bonds[0]) == 1


def test_len_beads_PEO(reader_PEO):
    repeats = 10
    end_beads = 2
    assert len(reader_PEO) - 1 == repeats + end_beads


def test_len_beads_PBD(reader_PBD):
    repeats = 10
    end_beads = 0
    assert len(reader_PBD) - 1 == repeats + end_beads


def test_len_beads_PBD(reader_PBD):
    repeats = 10
    end_beads = 0
    assert len(reader_PBD) - 1 == repeats + end_beads

def test_len_bonds_PEO(reader_PEO_bonds):
    repeats = 10
    end_beads = 2
    bonds = repeats + end_beads - 1
    angles = repeats + end_beads - 2
    dihedrals = repeats - 3
    assert len(reader_PEO_bonds) - 1 == bonds + angles + dihedrals

def test_len_bonds_PBD(reader_PBD_bonds):
    repeats = 10
    end_beads = 0
    bonds = repeats + end_beads - 1
    angles = repeats + end_beads - 2
    # dihedrals = repeats - 3
    assert len(reader_PBD_bonds) - 1 == bonds + angles

def test_len_bead_PBDbPEO(reader_PBDbPEO_bead):

    repeats_a = 10
    repeats_b = 10
    end_beads = 1
    assert len(reader_PBDbPEO_bead) - 1 == repeats_a + repeats_b+end_beads

def test_len_bead_PBDbPEO_continuous(reader_PBDbPEO_bead):
    for i,row in enumerate(reader_PBDbPEO_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i-1:
            check == False
    assert check

def test_len_bonds_PBDbPEO(reader_PBDbPEO_bonds):
    repeats_a = 10
    bonds_a = repeats_a-1
    angles_a = repeats_a-2
    bonds_block = 1
    angles_block = 2
    repeats_b = 10
    end_beads_b = 1
    bonds_b = repeats_b-1 + end_beads_b
    angles_b = repeats_b-2 + end_beads_b
    dihedrals_b = repeats_b-3
    tot_b = bonds_b+ angles_b + dihedrals_b
    tot_a = bonds_a + angles_a
    tot_block = bonds_block + angles_block
    assert len(reader_PBDbPEO_bonds) - 1 == tot_a+tot_b+tot_block

def test_len_bead_PCLbPEO(reader_PCLbPEO_bead):

    repeats_a = 10*2
    repeats_b = 10
    end_beads = 1
    assert len(reader_PCLbPEO_bead) - 1 == repeats_a + repeats_b+end_beads

def test_bond_continuity_PBDbPEO(reader_PBDbPEO_bonds,reader_PBDbPEO_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_PBDbPEO_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    num_bead = len(reader_PBDbPEO_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False
    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
    assert check_bead_1 and check_bead_2

def test_angle_continuity_PCLbPEO(reader_PCLbPEO_bonds,reader_PCLbPEO_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    bond_list_i3=[]

    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_PCLbPEO_bonds:
        if (row[2] != "" and row[3]==""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
            bond_list_i3.append(row[2])
    num_bead = len(reader_PCLbPEO_bead)-1
    check_bead_1 = 0
    for i in range(num_bead-2):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 +=1
    check_bead_2 = 0
    for i in range(1,num_bead-1):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 +=1
    check_bead_3 = 0
    for i in range(2,num_bead):
        found_index = False
        for bond_i3 in bond_list_i3:
            if int(bond_i3) == i:
                found_index = True
        if found_index==False:
            check_bead_3 +=1
    assert check_bead_1 + check_bead_2 + check_bead_3 ==6

def test_angles_in_bounds_PBDbPEO(reader_PBDbPEO_bonds,reader_PBDbPEO_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_PBDbPEO_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_PBDbPEO_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead:
            check_bead_1=False
    assert check_bead_1 


@pytest.fixture
def make_seq_def_1(root):
    name_a = "PEO"
    name_b = "PBD"
    id = "test"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "AA-BB-AAA-BBBB-AAA-BBB-AAA-BBB-AAA-BBB-AAAAAAAA"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_1_bead(root, make_seq_def_1):
    name_a = "PEO"
    name_b = "PBD"
    id = "test"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_1
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_1_bonds(root, make_seq_def_1):
    name_a = "PEO"
    name_b = "PBD"
    id = "test"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_1
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq1(reader_seq_def_1_bead):
    seq = "AA-BB-AAA-BBBB-AAA-BBB-AAA-BBB-AAA-BBB-AAAAAAAA"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    end_beads = 2
    assert len(reader_seq_def_1_bead) - 1 == len_seq-num_dashes+end_beads

def test_len_bead_seq1_continuous(reader_seq_def_1_bead):
    for i,row in enumerate(reader_seq_def_1_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq1(reader_seq_def_1_bonds,reader_seq_def_1_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_1_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_1_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq1(reader_seq_def_1_bonds,reader_seq_def_1_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_1_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_1_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
    assert check_bead_1 and check_bead_2


#test second sequence
@pytest.fixture
def make_seq_def_2(root):
    name_a = "PEO"
    name_b = "PBD"
    id = "test2"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "A-BB-AAA-BBBB-A-B-A"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_2_bead(root, make_seq_def_2):
    name_a = "PEO"
    name_b = "PBD"
    id = "test2"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_2
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_2_bonds(root, make_seq_def_2):
    name_a = "PEO"
    name_b = "PBD"
    id = "test2"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_2
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq2(reader_seq_def_2_bead):
    seq = "A-BB-AAA-BBBB-A-B-A"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    end_beads = 2
    assert len(reader_seq_def_2_bead) - 1 == len_seq-num_dashes+end_beads

def test_len_bead_seq1_continuous(reader_seq_def_2_bead):
    for i,row in enumerate(reader_seq_def_2_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq2(reader_seq_def_2_bonds,reader_seq_def_2_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_2_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_2_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq2(reader_seq_def_2_bonds,reader_seq_def_2_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_2_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_2_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
            
    assert check_bead_1 and check_bead_2



#test second sequence
@pytest.fixture
def make_seq_def_3(root):
    name_a = "PBD"
    name_b = "PEO"
    id = "test3"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "A-BB-AAA-BBBB-A-B-A"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_3_bead(root, make_seq_def_3):
    name_a = "PBD"
    name_b = "PEO"
    id = "test3"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_3
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_3_bonds(root, make_seq_def_3):
    name_a = "PBD"
    name_b = "PEO"
    id = "test3"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_3
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq3(reader_seq_def_3_bead):
    seq = "A-BB-AAA-BBBB-A-B-A"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    end_beads = 0
    assert len(reader_seq_def_3_bead) - 1 == len_seq-num_dashes+end_beads

def test_len_bead_seq3_continuous(reader_seq_def_3_bead):
    for i,row in enumerate(reader_seq_def_3_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq3(reader_seq_def_3_bonds,reader_seq_def_3_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_3_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_3_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq3(reader_seq_def_3_bonds,reader_seq_def_3_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_3_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_3_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
            
    assert check_bead_1 and check_bead_2

#test PCL b PEO
@pytest.fixture
def make_seq_def_4(root):
    name_a = "PEO"
    name_b = "PCL"
    id = "test4"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "A-BB-AAA-BBBB-A-B-A"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_4_bead(root, make_seq_def_4):
    name_a = "PEO"
    name_b = "PCL"
    id = "test4"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_4
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_4_bonds(root, make_seq_def_4):
    name_a = "PEO"
    name_b = "PCL"
    id = "test4"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_4
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq4(reader_seq_def_4_bead):
    seq = "A-BB-AAA-BBBB-A-B-A"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    num_PCL = 2*(2+4+1)
    num_PEO = 1+3+1+1
    end_beads = 2
    assert len(reader_seq_def_4_bead) - 1 == num_PCL+num_PEO+end_beads

def test_len_bead_seq4_continuous(reader_seq_def_4_bead):
    for i,row in enumerate(reader_seq_def_4_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq4(reader_seq_def_4_bonds,reader_seq_def_4_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_4_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_4_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq4(reader_seq_def_4_bonds,reader_seq_def_4_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_4_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_4_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
            
    assert check_bead_1 and check_bead_2


    #test PCL b PBD
@pytest.fixture
def make_seq_def_5(root):
    name_a = "PBD"
    name_b = "PCL"
    id = "test5"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "A-BB-AAA-BBBB-A-B-A"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_5_bead(root, make_seq_def_5):
    name_a = "PBD"
    name_b = "PCL"
    id = "test5"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_5
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_5_bonds(root, make_seq_def_5):
    name_a = "PBD"
    name_b = "PCL"
    id = "test5"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_5
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq5(reader_seq_def_5_bead):
    seq = "A-BB-AAA-BBBB-A-B-A"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    num_PCL = 2*(2+4+1)
    num_PEO = 1+3+1+1
    end_beads = 0
    assert len(reader_seq_def_5_bead) - 1 == num_PCL+num_PEO+end_beads

def test_len_bead_seq5_continuous(reader_seq_def_5_bead):
    for i,row in enumerate(reader_seq_def_5_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq5(reader_seq_def_5_bonds,reader_seq_def_5_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_5_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_5_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq5(reader_seq_def_5_bonds,reader_seq_def_5_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_5_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_5_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
            
    assert check_bead_1 and check_bead_2

def test_angle_continuity_seq5(reader_seq_def_5_bonds,reader_seq_def_5_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    bond_list_i3=[]

    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_5_bonds:
        if (row[2] != "" and row[3]==""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
            bond_list_i3.append(row[2])
    num_bead = len(reader_seq_def_5_bead)-1
    check_bead_1 = True
    for i in range(num_bead-2):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False
    check_bead_2 = True
    for i in range(1,num_bead-1):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
    check_bead_3 = True
    for i in range(2,num_bead):
        found_index = False
        for bond_i3 in bond_list_i3:
            if int(bond_i3) == i:
                found_index = True
        if found_index==False:
            check_bead_3 = False
    assert check_bead_1 and check_bead_2 and check_bead_3
#test PCL b PEO -> terminal bead = PCL
@pytest.fixture
def make_seq_def_6(root):
    name_a = "PEO"
    name_b = "PCL"
    id = "test6"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)

    path_bead = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bead.csv"
    )
    path_bonds = (
        root + "/projects/martini3/molecules/" + parsed_name+ "_bonds.csv"
    )
    seq = "A-BB-AAA-BBBB-AA-B-A-B"
    if os.path.exists(path_bead):
        os.remove(path_bead)
    if os.path.exists(path_bonds):
        os.remove(path_bonds)
    csv_names = polymers.make_sequence(name_a,name_b,seq,id)

@pytest.fixture
def reader_seq_def_6_bead(root, make_seq_def_6):
    name_a = "PEO"
    name_b = "PCL"
    id = "test6"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bead" + ".csv"
    )
    make_seq_def_6
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   
@pytest.fixture
def reader_seq_def_6_bonds(root, make_seq_def_6):
    name_a = "PEO"
    name_b = "PCL"
    id = "test6"
    parsed_name = "seq_" + name_a + "_"+ name_b + "_" + str(id)
    path = (
        root + "/projects/martini3/molecules/" + parsed_name + "_bonds" + ".csv"
    )
    make_seq_def_6
    with open(path, "r", newline="") as read_file:
        reader = csv.reader(read_file)
        return list(reader)
   


def test_len_bead_seq6(reader_seq_def_6_bead):
    seq = "A-BB-AAA-BBBB-AA-B-A-B"
    len_seq = len(seq)
    num_dashes = len(seq.split("-"))-1
    num_PCL = 2*(2+4+1+1)
    num_PEO = 1+3+1+2
    end_beads = 1
    assert len(reader_seq_def_6_bead) - 1 == num_PCL+num_PEO+end_beads

def test_len_bead_seq6_continuous(reader_seq_def_6_bead):
    for i,row in enumerate(reader_seq_def_6_bead):
        if i ==0:
            check = True
        elif int(row[0]) != i:
            check = False
    assert check

def test_angles_in_bounds_seq6(reader_seq_def_6_bonds,reader_seq_def_6_bead):
    angle_list_i3 = []
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_6_bonds:
        if (row[2] != "" and row[3] == ""):
            angle_list_i3.append(row[1])
    num_bead = len(reader_seq_def_6_bead)-1
    check_bead_1 = True
    for i in angle_list_i3:
        if int(i)>num_bead-1:
            check_bead_1=False
    assert check_bead_1 


def test_bond_continuity_seq6(reader_seq_def_6_bonds,reader_seq_def_6_bead):
    bond_list_i1 = []
    bond_list_i2=[]
    #if every number shows up in both the first index and 2nd index (besides end bead),should have connectivity
    for row in reader_seq_def_6_bonds:
        if (row[2] == ""):
            bond_list_i1.append(row[0])
            bond_list_i2.append(row[1])
    print(bond_list_i1)
    print(bond_list_i2)
    num_bead = len(reader_seq_def_6_bead)-1
    check_bead_1 = True
    for i in range(num_bead-1):
        found_index = False
        for bond_i1 in bond_list_i1:
            if int(bond_i1) == i:
                found_index = True
        if found_index==False:
            check_bead_1 = False

    check_bead_2 = True
    for i in range(1,num_bead):
        found_index = False
        for bond_i2 in bond_list_i2:
            if int(bond_i2) == i:
                found_index = True
        if found_index==False:
            check_bead_2 = False
            
    assert check_bead_1 and check_bead_2
