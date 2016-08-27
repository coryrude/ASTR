# RUNORBIT.PY: run orbit in 3-D logarithmic potential
# Written by Scott Griffiths, 17 Apr 2014
# Based on RUNORBIT.C by Josh Barnes, 12 Feb 2005

import argparse
import numpy as np
import matplotlib.pyplot as plt
from os import makedirs, path
    

class OrbitIntegrator:
    def __init__(self, **kwargs):
        # set default parameters if keyword not passed
        pos = kwargs.pop('pos', [1.0, 0.0, 0.0])
        vel = kwargs.pop('vel', [0.0, 0.4, 0.0])
        nstep = kwargs.pop('nstep', 25000)
        self.dtime = kwargs.pop('dtime', 0.01)
        detol = kwargs.pop('detol', 1.0e-3)
        self.Rc = kwargs.pop('Rc', 0.2)
        self.b =  kwargs.pop('b', 0.9)
        self.c =  kwargs.pop('c', 0.8)
        self.vy = vel[1]
        
        # compute initial energy and save initial position
        pos = np.array(pos); vel = np.array(vel)
        Einit = self._energy(pos, vel)
        self.orbit = [np.append(pos, 0.0)]
        
        # loop nstep times...
        for n in xrange(1, nstep+1):                      
            pos,vel = self._leapstep(pos, vel)          # take a leapfrog step     
            tnow = self.dtime * n                       # advance value of time    
            Enow = self._energy(pos, vel)               # compute current energy   
            if np.abs(Einit - Enow) > detol:            # test energy conservation 
                print "Energy conservation violated:"
                print "  E(0) = %f,  E(%f) = %f,  Delta E = %f" % (Einit, tnow, Enow, Einit - Enow)
                break
            self.orbit.append(np.append(pos, tnow))     # save position and time
        print "[ run completed : E(0) = %f,  E(%f) = %f ]" % (Einit, tnow, Enow)

    def plot(self, plane='XY'):
        plt.ion()
        plt.clf()
        orbit = np.array(self.orbit)
        x,y,z,t = np.transpose(orbit)
        { 'XY' : plt.plot(x, y, 'k-'),
          'YZ' : plt.plot(y, z, 'k-'),
          'XZ' : plt.plot(x, z, 'k-') }.get(plane)
        plt.xlabel(plane[0])
        plt.ylabel(plane[1])
        plt.title("Initial y velocity = {0}".format(self.vy))
    
    def write(self, file):
        with open(file, 'w') as f:
            f.write("#{0:>11s} {1:>11s} {2:>11s} {3:>11s}\n".format("x", "y", "z", "time"))
            for coords in self.orbit:
                f.write(" {0:11.6f} {1:11.6f} {2:11.6f} {3:11.6f}\n".format(*coords))
        
    def _energy(self, pos, vel):
        """ENERGY: compute binding energy from phase-space coordinates."""
        # return sum of KE and PE
        return (0.5*np.sum(np.square(vel)) + self._potential(pos))

    def _leapstep(self, pos, vel):
        """LEAPSTEP: advance coordinates using leap-frog integrator."""
        acc = self._acceleration(pos)       # compute starting acceleration
        vel = vel + 0.5*self.dtime*acc      # advance velocities a half-step
        pos = pos + self.dtime*vel          # advance positions a whole step
        acc = self._acceleration(pos)       # compute ending acceleration
        vel = vel + 0.5*self.dtime*acc      # complete velocity step
        return pos, vel

    def _potential(self, pos):
        """POTENTIAL: compute gravitational potential of logarithmic model."""
        x,y,z = pos
        mu2 = np.sum(np.square([self.Rc, x, y/self.b, z/self.c]))
        return 0.5*np.log(mu2)

    def _acceleration(self, pos):
        """ACCELERATION: compute gravitational acceleration of logarithmic model."""
        x,y,z = pos
        mu2 = np.sum(np.square([self.Rc, x, y/self.b, z/self.c]))
        acc = -np.array([x, y/(self.b**2), z/(self.c**2)])/mu2
        return acc

    
if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description='RUNORBIT.PY: run orbit in 3-D logarithmic potential')
    parser.add_argument('--pos', type=float, nargs=3, default=[1.0, 0.0, 0.0], help='Initial position vector')
    parser.add_argument('--vel', type=float, nargs=3, default=[0.0, 0.4, 0.0], help='Initial velocity vector')
    parser.add_argument('--nstep', type=int, default=25000, help='Number of time-steps')
    parser.add_argument('--dtime', type=float, default=0.01, help='Value of time-step')
    parser.add_argument('--detol', type=float, default=1.0e-3, help='Tolerance for energy error')
    parser.add_argument('-R', type=float, default=0.2, help='Potential core radius')
    parser.add_argument('-b', type=float, default=0.9, help='Y:X axis ratio of potential')
    parser.add_argument('-c', type=float, default=0.8, help='Z:X axis ratio of potential')
    args = parser.parse_args()
    orbit = OrbitIntegrator(**vars(args))
    orbit.plot()
    
    # Save results
    if not path.exists('RunOrbit'):
        makedirs('RunOrbit')
    fname = "RunOrbit/vy_{0}".format(str(args.vel[1]).replace('.',','))
    plt.savefig(fname+'.png')
    orbit.write(fname+'.dat')
