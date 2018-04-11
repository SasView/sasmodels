#!/usr/bin/env python
"""
Evaluate the scattering for a shape using spherical integration over the entire
surface in theta-phi polar coordinates, and compare it to the 1D scattering
pattern for that shape.

Parameters are the same as sascomp.  Unlike the -sphere option in sascomp,
the evaluation is restricted to a single radial line for performance reasons,
with angle set by -angle=alpha in the qx-qy plane.
"""

from __future__ import print_function, division

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import numpy as np

from sasmodels import compare, data

def main(angle=0, steps=76):
    # Parse the options using the parser in sasmodels.compare.  -angle
    # is an additional parameter that is parsed separately., with -sphere=n and -angle=a as additional arguments
    # Pull out angle ans
    argv = []
    for arg in sys.argv[1:]:
        if arg.startswith('-sphere='):
            steps = int(arg[8:])
        elif arg.startswith('-angle='):
            angle = float(arg[7:])
        elif arg != "-2d":
            argv.append(arg)
    opts = compare.parse_opts(argv)

    # Create a 2D radial slice
    qr = opts['data'].x
    qx, qy = qr*np.cos(np.radians(angle)), qr*np.sin(np.radians(angle))
    radial_data = data.Data2D(x=qx, y=qy)
    radial_data.radial = True  # let plotter know it is actual 1D

    # Create an engine to evaluate it
    comp = compare.make_engine(opts['info'][0], radial_data,
                               opts['engine'][0], opts['cutoff'][0])
    opts['engines'] = [opts['engines'][0], comp]

    # Set the integration parameters to the half sphere
    compare.set_spherical_integration_parameters(opts, steps)

    # Set the random seed
    if opts['seed'] > -1:
        print("Randomize using -random=%i"%opts['seed'])
        np.random.seed(opts['seed'])

    # Run the comparison
    compare.compare(opts, maxdim=np.inf)

if __name__ == "__main__":
    main(angle=0)
