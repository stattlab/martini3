import pytest
import sys
import os
sys.path.append('../')

from martini3 import molecules as m

# Test Bond
@pytest.fixture
def bond_1():
  return m.Bond(0,[0,1],.3,5000)

@pytest.fixture
def bond_diff_idx():
  return m.Bond(0,[1,2],.3,5000)

@pytest.fixture
def bond_diff_force():
  return m.Bond(1,[1,2],.3,3000)

@pytest.fixture
def bond_4():
  return m.Bond(2,[1,2],.5,3500)

@pytest.fixture
def bond_list(bond_1,bond_diff_idx,bond_diff_force):
  return [bond_1,bond_diff_idx,bond_diff_force]

def test_bond_init(bond_1):
  assert bond_1.spatial == .3 and bond_1.bead_indices == [0,1] and bond_1.force ==5000 and bond_1.idx ==0
  
def test_bond_update(bond_1):
  bond_1.update_idx(1)
  assert bond_1.spatial == .3 and bond_1.bead_indices == [0,1] and bond_1.force ==5000 and bond_1.idx ==1

def test_bond_equal(bond_1,bond_diff_idx):
  assert bond_1==bond_diff_idx

def test_bond_not_equal(bond_1,bond_diff_force):
  assert bond_1 != bond_diff_force

def test_bond_hash_not_equal(bond_1,bond_diff_force):
  assert (hash(bond_1)) != hash(bond_diff_force)

def test_bond_hash_equal(bond_1,bond_diff_idx):
  assert (hash(bond_1)) == hash(bond_diff_idx)

def test_bond_hash_formula(bond_1,bond_diff_idx):
  assert (hash(bond_1)) == hash((.3,5000))

def test_find_existing_bond(bond_1,bond_list):
  assert bond_1 == bond_1.find_existing_bond(bond_list)

def test_find_existing_bond_none(bond_4,bond_list):
  assert None == bond_4.find_existing_bond(bond_list)

def test_count_unique_bonds(bond_list):
  assert m.count_unique_bonds(bond_list) == 2


# Test Angle
@pytest.fixture
def angle_1():
  return m.Angle(0,[0,1,2],.3,5000)

@pytest.fixture
def angle_diff_idx():
  return m.Angle(0,[1,2,3],.3,5000)

@pytest.fixture
def angle_diff_force():
  return m.Angle(1,[1,2,3],.3,3000)

@pytest.fixture
def angle_4():
  return m.Angle(2,[1,2,3],.5,3500)

@pytest.fixture
def angle_list(angle_1,angle_diff_idx,angle_diff_force):
  return [angle_1,angle_diff_idx,angle_diff_force]

def test_angle_init(angle_1):
  assert angle_1.spatial == .3 and angle_1.bead_indices == [0,1,2] and angle_1.force ==5000 and angle_1.idx ==0
  
def test_angle_update(angle_1):
  angle_1.update_idx(1)
  assert angle_1.spatial == .3 and angle_1.bead_indices == [0,1,2] and angle_1.force ==5000 and angle_1.idx ==1

def test_angle_equal(angle_1,angle_diff_idx):
  assert angle_1==angle_diff_idx

def test_angle_not_equal(angle_1,angle_diff_force):
  assert angle_1 != angle_diff_force

def test_angle_hash_not_equal(angle_1,angle_diff_force):
  assert (hash(angle_1)) != hash(angle_diff_force)

def test_angle_hash_equal(angle_1,angle_diff_idx):
  assert (hash(angle_1)) == hash(angle_diff_idx)

def test_angle_hash_formula(angle_1,angle_diff_idx):
  assert (hash(angle_1)) == hash((.3,5000))

def test_find_existing_angle(angle_1,angle_list):
  assert angle_1 == angle_1.find_existing_angle(angle_list)

def test_find_existing_angle_none(angle_4,angle_list):
  assert None == angle_4.find_existing_angle(angle_list)

def test_count_unique_angles(angle_list):
  assert m.count_unique_angles(angle_list) == 2

# Test Dihedrals
@pytest.fixture
def dihedral_1():
  return m.Dihedral(0,[0,1,2,3],.3,5000,3,4)

@pytest.fixture
def dihedral_diff_idx():
  return m.Dihedral(0,[1,2,3,4],.3,5000,3,4)

@pytest.fixture
def dihedral_diff_force():
  return m.Dihedral(1,[1,2,3,4],.3,3000,3,4)

@pytest.fixture
def dihedral_4():
  return m.Dihedral(2,[1,2,3,4],.5,3500,4,6)

@pytest.fixture
def dihedral_list(dihedral_1,dihedral_diff_idx,dihedral_diff_force):
  return [dihedral_1,dihedral_diff_idx,dihedral_diff_force]

def test_dihedral_init(dihedral_1):
  assert dihedral_1.k1 == .3 and dihedral_1.bead_indices == [0,1,2,3] and dihedral_1.k2 ==5000 and dihedral_1.idx ==0  and dihedral_1.k3 ==3 and dihedral_1.k4==4
  
def test_dihedral_update(dihedral_1):
  dihedral_1.update_idx(1)
  assert dihedral_1.k1 == .3 and dihedral_1.bead_indices == [0,1,2,3] and dihedral_1.k2 ==5000 and dihedral_1.idx ==1 and dihedral_1.k3 ==3 and dihedral_1.k4==4

def test_dihedral_equal(dihedral_1,dihedral_diff_idx):
  assert dihedral_1==dihedral_diff_idx

def test_dihedral_not_equal(dihedral_1,dihedral_diff_force):
  assert dihedral_1 != dihedral_diff_force

def test_dihedral_hash_not_equal(dihedral_1,dihedral_diff_force):
  assert (hash(dihedral_1)) != hash(dihedral_diff_force)

def test_dihedral_hash_equal(dihedral_1,dihedral_diff_idx):
  assert (hash(dihedral_1)) == hash(dihedral_diff_idx)

def test_dihedral_hash_formula(dihedral_1):
  assert (hash(dihedral_1)) == hash((.3,5000,3,4))

def test_find_existing_dihedral(dihedral_1,dihedral_list):
  assert dihedral_1 == dihedral_1.find_existing_dihedral(dihedral_list)

def test_find_existing_dihedral_none(dihedral_4,dihedral_list):
  assert None == dihedral_4.find_existing_dihedral(dihedral_list)

def test_count_unique_dihedrals(dihedral_list):
  assert m.count_unique_dihedrals(dihedral_list) == 2


# Test Contents
@pytest.fixture
def contents():
  return m.Contents()

def test_init_contents(contents):
  assert len(contents.contents)==0

# Test Molecules
@pytest.fixture
def water(contents):
  w = m.Molecule('W','W','W',contents)
  return w 

def test_water_types(water):
  assert water.types == ["W"]

def test_water_index(water):
  assert water.id ==[0]

def test_water_name(water):
  assert water.name == "W"

def test_water_position(water):
  assert water.position == [[0.0,0.0,0.0]]

def test_water_shiftposition(water):
  water.shift_positions(0,0,1)
  print(water.position)
  assert water.position == [[0.0,0.0,1.0]]

def test_water_invertposition(water):
  water.shift_positions(0,0,1)
  water.invert_positions(True)
  assert water.position == [[0.0,0.0,-1.0]]

def test_water_invertposition_false(water):
  water.shift_positions(0,0,1)
  water.invert_positions(False)
  assert water.position == [[0.0,0.0,1.0]]

def test_water_expandposition(water):
  water.shift_positions(0,0,1)
  water.expand(5)
  assert water.position == [[0.0,0.0,5.0]]

def test_water_bonded(water):
  a = len(water.bonds)
  b = len(water.angles)
  c = len(water.dihedrals)
  assert a+b+c == 0

def test_water_charges(water):
  assert water.charges ==[None]

@pytest.fixture
def tm(contents):
  bead_path, bond_path = m.path_to_beads("test")
  tm = m.Molecule("test", bead_path, bond_path, contents)
  return tm


def test_tm_types(tm):
  assert tm.types == ["Q1",'Q5',"SN4a","TN4a","C1",'C4h']

def test_tm_index(tm):
  assert tm.id == ['1','2','3','4','5','6']

def test_tm_name(tm):
  assert tm.name == "test"

def test_tm_position(tm):
  assert tm.position == [[0.0,0.0,0.0],[0.0,0.0,1.0],[1.0,0.0,1.0],[2.0,0.0,1.0],[2.0,1.0,1.0],[2.0,2.0,1.0]]

def test_tm_charges(tm):
  assert tm.charges ==['1','-1',None,None,None,None]


  #Next, test bonds. 

  #Test Path_to_beads 
  # tests contents.add molecule
# test extract molecules