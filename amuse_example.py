#!/usr/bin/env python
'''
A simple example/test of the AMUSE interface for falcON:
evolve a 10^5 particle Plummer sphere for a short time
'''
import time
from amuse.units import nbody_system
from amuse.community.falcon.interface import Falcon
#from amuse.community.bhtree.interface import BHTree
from amuse.ic.plummer import new_plummer_model

particles = new_plummer_model(100000)
nbody = Falcon(channel_type='sockets')
#nbody = BHTree(channel_type='sockets')  # this is 10-20x slower than falcon
nbody.particles.add_particles(particles)
nbody.parameters.epsilon_squared = (0.02 | nbody_system.length)**2
nbody.parameters.timestep = 0.0625 | nbody_system.time
wallclock_t0=time.time()
for i in range(33):
    Ekin = nbody.get_kinetic_energy  ().value_in(nbody_system.mass * nbody_system.length**2 * nbody_system.time**-2)
    Epot = nbody.get_potential_energy().value_in(nbody_system.mass * nbody_system.length**2 * nbody_system.time**-2)
    print("T=%6.4g, Epot=%10.7f, Ekin=%10.7f, Etot=%10.7g, wallclock time=%5.2f s" %
        (nbody.get_time().value_in(nbody_system.time), Epot, Ekin, Epot+Ekin, time.time()-wallclock_t0))
    if i==32: break
    nbody.evolve_model((i+1) * nbody.parameters.timestep)
nbody.stop()

