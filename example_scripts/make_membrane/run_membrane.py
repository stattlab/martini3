import numpy as np
import sys
from martini3 import particles
from martini3 import force_fields
from martini3 import molecules
from martini3 import init_cell
import datetime
import hoomd
import os

#Writes outputs to terminal about progress/remaining time
class Status:
    def __init__(self, simulation):
        self.simulation = simulation

    @property
    def seconds_remaining(self):
        try:
            return (
                self.simulation.final_timestep - self.simulation.timestep
            ) / self.simulation.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))


def main():
    path = "data/" + "DPPC_membrane/"
    if not os.path.exists(path):
        raise Exception("no data present")
    
    #init all force fields
    lj, coulomb, bond_harmonic, angle_forces, dihedrals,impropers,rigid = force_fields.forces_from_gsd(
        path, "init.gsd"
    )
    # create sim with GPU if you have it otherwise CPU
    try:
        sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=int(50) * 500)
    except:
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=int(50) * 500)

    #Load init.gsd
    name = "init.gsd"
    sim.create_state_from_gsd(filename=path + name)
    # thermostat and barostat correctly
    velocity_rescaling = 1
    pressure_rescaling = 12
    pressure = 0.061019  # 1 unit is 16.388 atm
    mttk = hoomd.md.methods.thermostats.MTTK(kT=2.47, tau=velocity_rescaling)
    rigid_centers_and_free = hoomd.filter.All()
    cp = hoomd.md.methods.ConstantPressure(
        filter=rigid_centers_and_free,
        tauS=pressure_rescaling,
        S=[0, 0, pressure, 0, 0, 0],
        couple="xy",box_dof=[True,True,True,False,False,False],
        thermostat=mttk,
    )

    #assign to integrator
    integrator = hoomd.md.Integrator(
        dt=0.02,
        methods=[cp],
        forces=[lj, bond_harmonic, angle_forces, coulomb],
    )
    #Write to terminal
    status = Status(sim)
    logger = hoomd.logging.Logger(categories=["scalar", "string"])
    logger.add(sim, quantities=["timestep", "tps"])
    logger[("Status", "etr")] = (status, "etr", "string")
    table = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(period=50000), logger=logger
    )
    sim.operations.writers.append(table)

    #write to gsd
    gsd_writer = hoomd.write.GSD(
        filename=path + "traj.gsd",
        trigger=hoomd.trigger.Periodic(1000),
        filter=hoomd.filter.All(),
        mode="wb",
    )
    sim.operations.integrator = integrator
    sim.operations.writers.append(gsd_writer)

    sim.run(100000)



if __name__ == "__main__":
    if len(sys.argv) != 1:
        print("Usage: python3 run_membrane.py")
        sys.exit(1)
    main()
