import numpy as np
import sys
from martini3 import molecules
from martini3 import init_cell
import hoomd
import os


def gen_grid(box, num_spaces=100):
    grid = []
    x_arr = np.linspace(-box[0]/2, box[0]/2, num_spaces)
    y_arr = np.linspace(-box[1]/2, box[1]/2, num_spaces)

    for x in x_arr:
        for y in y_arr:
            grid.append((x, y))
    return grid


def make_periodic(positions,box):
    min_image = np.where(positions > 0.5*box, positions - box, np.where(positions < -0.5*box, positions + box , positions))

    return np.array(min_image)

def check_dist(molecule_position, positions,box):
    threshold = .355
    distance_matrix = np.linalg.norm(make_periodic(molecule_position[:, np.newaxis, :] - positions,box), axis=2)    
    return np.all(distance_matrix> threshold)

def main():
    path = "data/" +"DPPC_membrane/" 
    if not os.path.exists(path):
        os.makedirs(path)
    # Initialize simulation box
    x_box_size = 10
    y_box_size= 10
    z_box_size = 10
    box = np.array([x_box_size,y_box_size,z_box_size])


    # Compute total number of lipids to be placed -- can replace total lipids with a specific number)
    initial_density = 0.69  # since we cant realistically fill 100%
    total_h_beads = (
        2 * 9 * x_box_size*y_box_size * 1.4 * initial_density
    )  # this is using the fact that the density of lipids is ~1.4 per nm**2 at RT, and theres ab 2*9 beads in each layer
    total_lipids = int(total_h_beads / 9)
    lipid_placed = 0
    positions = np.empty((0, 3))
    np.random.seed()
    grid_lipid = gen_grid(box, num_spaces=x_box_size * 40)
    is_inverted=0

    #init water grid (when placing lipids, remember lowest and highest places you have placed molecules. Do not place water in these regions)
    num_water_x = int(x_box_size*2.1)
    num_water_y = int(y_box_size*2.1)
    grid_upper = np.zeros((num_water_x,num_water_y)) -20
    grid_lower = np.zeros((num_water_x,num_water_y)) + 20

    x_range = np.linspace(-x_box_size/2,x_box_size/2,num_water_x)
    y_range = np.linspace(-y_box_size/2,y_box_size/2,num_water_y)
    contents = molecules.Contents()

    # Place lipids
    while lipid_placed < total_lipids:
        rand_pos_lipid = int(np.random.randint(0, len(grid_lipid)))
        x_shift, y_shift = grid_lipid[rand_pos_lipid]
        theta = np.random.uniform(0, np.pi * 2)
        z_shift = -0.4
        
        molecule = molecules.make_DPPC(
            contents,
            x_shift=x_shift,
            y_shift=y_shift,
            z_shift=z_shift,
            theta=theta,
            is_inverted=is_inverted,
        )
        del grid_lipid[rand_pos_lipid]
    
        molecule.position = make_periodic(molecule.position,box)
        molecule_position = np.array(molecule.position)
        if len(positions) == 0 or (check_dist(molecule_position,positions,box)
        ):
            contents.add_molecule(molecule)
            positions = np.append(positions, molecule_position, axis=0)
            for pos in positions:
                x_pos = pos[0]
                y_pos = pos[1]
                z_pos = pos[2]
                closest_x = np.abs(x_range-x_pos).argmin()
                closest_y = np.abs(y_range-y_pos).argmin()
                if z_pos > grid_upper[closest_x,closest_y]:
                    grid_upper[closest_x,closest_y] = z_pos
                elif z_pos < grid_lower[closest_x,closest_y]:
                    grid_lower[closest_x,closest_y] = z_pos
            lipid_placed = lipid_placed + 1
            if is_inverted ==0 :
                is_inverted = 1
            else:
                is_inverted=0

    # Place water

    x_place_water = np.linspace(
        -x_box_size / 2 , x_box_size / 2 , int(x_box_size * 2.1)
    )
    y_place_water = np.linspace(
        -y_box_size / 2, y_box_size / 2, int(y_box_size * 2.1)
    )
    z_place_water = np.linspace(
        -z_box_size / 2 + 0.2, z_box_size / 2 - 0.25, int(z_box_size * 2.2)
    )

    for x_iter2, x2 in enumerate(x_place_water):
        for y_iter2, y2 in enumerate(y_place_water):
            for z2 in z_place_water:
                dont_place = False
                if x_iter2 > 1 and x_iter2 < len(x_place_water)-2 and y_iter2 > 1 and y_iter2 < len(y_place_water)-2:
                    for l in range(5):
                        for k in range(5):
                            if z2 > grid_lower[x_iter2+l-2,y_iter2+k-2]-.3 and z2<grid_upper[x_iter2+l-2,y_iter2+k-2]+.3:
                                dont_place = True
                else:
                    dont_place = True
                if not dont_place:
                    contents = molecules.add_water(
                        contents, x_shift=x2+np.random.randint()*.05, y_shift=y2, z_shift=z2)
    
    # Make init.gsd
    lj, coulomb, bond_harmonic, angle_forces, _,_,_ = init_cell.init_cell(
        contents, path, box_size=[x_box_size, y_box_size, z_box_size], pair_on=False
    )


if __name__ == "__main__":
    if len(sys.argv) != 1:
        print("Usage: python3 init_small.py")
        sys.exit(1)
    main()
