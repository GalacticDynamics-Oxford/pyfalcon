'''
Test of pyfalcon community module for AMUSE
'''
from amuse.test.amusetest import TestWithMPI
from amuse.support import exceptions
from amuse.units import nbody_system, units
from amuse.community.falcon.interface import Falcon
from amuse.datamodel import Particles

class TestFalcon(TestWithMPI):

    def test0(self):
        instance = Falcon()
        particles = Particles(2)
        particles.radius = 0.5 | nbody_system.length
        particles[0].mass= 1.0 | nbody_system.mass
        particles[1].mass= 2.0 | nbody_system.mass
        particles[0].x   =-2.0 | nbody_system.length
        particles[1].x   = 1.0 | nbody_system.length
        particles[0].vy  = 0.5 | nbody_system.speed
        particles[1].vy  =-0.25| nbody_system.speed
        particles.y  = particles.z  = 0 | nbody_system.length
        particles.vx = particles.vz = 0 | nbody_system.speed
        instance.particles.add_particles(particles)
        instance.parameters.individual_epsilon = True
        # evolve_model should fail because timestep has not been set
        self.assertRaises(exceptions.AmuseException, lambda: instance.evolve_model(1 | nbody_system.time))
        # now set timestep and try again
        instance.parameters.timestep = 0.25 | nbody_system.time
        E0 = instance.get_kinetic_energy() + instance.get_potential_energy()
        instance.evolve_model(10 | nbody_system.time)
        E1 = instance.get_kinetic_energy() + instance.get_potential_energy()
        self.assertLess(abs(1-E1/E0), 2e-5)
        # without individual softening lengths, need to set the global one, otherwise it should fail
        instance.parameters.individual_epsilon = False
        self.assertRaises(exceptions.AmuseException, lambda: instance.evolve_model(1 | nbody_system.time))
        instance.parameters.epsilon = particles.radius  # set the global softening length
        instance.evolve_model(1 | nbody_system.time)    # now should work
        instance.stop()


if __name__ == '__main__':
    try:
        test = TestFalcon('test0')
        test.test0()
        print('PASS')
    except Exception as e:
        print('FAIL: %s %s' % (type(e), e))
