#!/usr/bin/env python
info = '''Python interface for the falcON N-body simulation code from NEMO.
Command-line arguments (same meaning as in gyrfalcON):
  in     name of the input snapshot file in the NEMO format
  out    name of the output snapshot file
  eps    softening length
  kmax   timestep parameter: dt = 0.5^kmax
  step   interval between writing output snapshots
  tstop  total simulation time
'''
import numpy, pyfalcon, os, sys, time

class Snapshot:
    '''Container for position, velocity, mass, softening length, acceleration'''

    def __init__(self, pos, vel, mass, eps, time=0.0):
        self.pos = pos  # positions  (Nx3)
        self.vel = vel  # velocities (Nx3)
        self.mass= mass # masses (N or a single number)
        self.eps = eps  # softening lengths (N or a single number)
        # initialize the accelerations (Nx3) and potential (N, not used in any computations)
        self.acc, self.pot = pyfalcon.gravity(self.pos, self.mass, self.eps)
        self.time = time

    def leapfrog(self, dt):
        '''Perform one step of a basic kick-drift-kick leapfrog method'''
        # kick for half-step, using accelerations computed at the end of the previous step
        self.vel += self.acc * (dt/2)
        # drift for full step
        self.pos += self.vel * dt
        # recompute accelerations
        self.acc, self.pot = pyfalcon.gravity(self.pos, self.mass, self.eps)
        # kick again for half-step
        self.vel += self.acc * (dt/2)
        self.time += dt


### the following lengthy piece of code implements input/output of N-body snapshots in the NEMO format
try:
    import unsio
except ImportError:
    try: import py_unsio as unsio  # old name
    except ImportError: unsio = None

if unsio is not None:
    def readSnapshot(filename, eps):
        '''Read a NEMO snaphot using unsio'''
        if not os.path.isfile(filename):
            raise RuntimeError('Input file %s is not found' % filename)
        snapfile = unsio.CunsIn(filename, 'all', 'all')
        snapfile.nextFrame('')
        pos = snapfile.getArrayF('all', 'pos')[1].reshape(-1,3).copy()
        vel = snapfile.getArrayF('all', 'vel')[1].reshape(-1,3).copy()
        mass= snapfile.getArrayF('all', 'mass')[1].copy()
        time= snapfile.getValueF('time')[1]
        return Snapshot(pos, vel, mass, eps, time)

    def writeSnapshot(filename, snapshot):
        '''Write a NEMO snapshot using unsio'''
        if os.path.exists(filename):  # UNSIO, like NEMO, will refuse to overwrite an existing file
            raise RuntimeError('Output file %s exists' % filename)
        snapfile = unsio.CunsOut(filename, 'nemo')
        snapfile.setValueF('time', snapshot.time)
        snapfile.setArrayF('all', 'pos', snapshot.pos. astype(numpy.float32).reshape(-1))
        snapfile.setArrayF('all', 'vel', snapshot.vel. astype(numpy.float32).reshape(-1))
        snapfile.setArrayF('all', 'mass',snapshot.mass.astype(numpy.float32))
        snapfile.setArrayF('all', 'acc', snapshot.acc. astype(numpy.float32).reshape(-1))
        snapfile.setArrayF('all', 'pot', snapshot.pot. astype(numpy.float32))
        snapfile.save()

else:
    print("UNSIO is not available; using a poor man's substitute for I/O (*much* slower)")
    import subprocess

    def readSnapshot(filename, eps):
        '''Read NEMO snapshot using NEMO snapprint program'''
        data = numpy.loadtxt(subprocess.Popen('snapprint %s x,y,z,vx,vy,vz,m' % filename,
            shell=True, stdout=subprocess.PIPE).communicate()[0].decode().splitlines()
            ).reshape(-1,7).astype(numpy.float32)
        return Snapshot(data[:,0:3], data[:,3:6], data[:,6], eps)

    def writeSnapshot(filename, snapshot):
        '''Write a NEMO snapshot using NEMO a2s program'''
        if os.path.exists(filename):  # any NEMO program will refuse to overwrite an existing file
            raise RuntimeError('Output file %s exists' % filename)
        data = numpy.column_stack((snapshot.pos, snapshot.vel, snapshot.mass, snapshot.acc, snapshot.pot))
        stdin = '\n'.join([' '.join(x.astype(str)) for x in data])
        subprocess.Popen('a2s - %s read=x,v,m,a,p N=%i time=%g' % (filename, len(data), snap.time),
            shell=True, stdin=subprocess.PIPE).communicate(stdin.encode())


# finally comes the main program
if __name__ == '__main__':
    # parse input parameters
    arglist = []
    for iarg, arg in enumerate(sys.argv[1:]):
        nameval = arg.split('=')
        if len(nameval)!=2:
            # first two arguments are 'in' and 'out', and their names need not be specified explicitly
            if   iarg==0: nameval=('in', arg)
            elif iarg==1: nameval=('out',arg)
            else: raise ValueError('Command-line arguments should be in the form  name=value')
        arglist.append([nameval[0], nameval[1]])
    args = dict(arglist)
    infile  = args.get('in')    # input filename (required)
    outfile = args.get('out')   # output filename (required)
    eps     = args.get('eps')   # softening length (required)
    kmax    = args.get('kmax')  # timestep parameter (required)
    if infile is None or outfile is None or eps is None or kmax is None:
        print(info)
        exit()
    eps   = float(eps)
    dt    = 0.5**int(kmax)
    tstop = float(args.get('tstop', numpy.inf))  # total integration time
    step  = float(args.get('step',  numpy.inf))  # interval between outputs
    snap  = readSnapshot(infile, eps)

    # run the integration loop
    wallclock_t0 = time.time()
    print("#   sim_time   total_energy   pot_energy   kin_energy   virial   cpu_time")
    while snap.time <= tstop:
        if snap.time % step == 0:  # write a snapshot to a file, appending time to the filename
            writeSnapshot(outfile + '_t%g' % snap.time, snap)
        # print out diagnostic info
        Epot = numpy.sum(snap.mass[:,None] * snap.pos * snap.acc)
        Ekin = numpy.sum(snap.mass[:,None] * snap.vel**2) * 0.5
        wallclock_t = time.time() - wallclock_t0
        print("%12.6f %14.8g %12.6g %12.6g %8.4f %2i:%02i:%04.1f" %
            (snap.time, Epot+Ekin, Epot, Ekin, -2*Ekin/Epot,
            wallclock_t/3600, wallclock_t/60%60, wallclock_t%60))
        if snap.time >= tstop: break
        # perform one step of simulation
        snap.leapfrog(dt)
