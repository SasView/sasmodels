
from numpy import pi, inf

name = "CylinderModel"
title = "Cylinder with uniform scattering length density"
description = """\
         f(q)= 2*(sldCyl - sldSolv)*V*sin(qLcos(alpha/2))
                /[qLcos(alpha/2)]*J1(qRsin(alpha/2))/[qRsin(alpha)]

                P(q,alpha)= scale/V*f(q)^(2)+background
                V: Volume of the cylinder
                R: Radius of the cylinder
                L: Length of the cylinder
                J1: The bessel function
                alpha: angle betweenthe axis of the
                cylinder and the q-vector for 1D
                :the ouput is P(q)=scale/V*integral
                from pi/2 to zero of...
                f(q)^(2)*sin(alpha)*dalpha+ bkg
"""

parameters = [
    #   [ "name", "units", default, [lower, upper], "type",
    #     "description" ],
    [ "sldCyl", "1/Ang^2", 4e-6, [-inf,inf], "",
      "Cylinder scattering length density" ],
    [ "sldSolv", "1/Ang^2", 1e-6, [-inf,inf], "",
      "Solvent scattering length density" ],
    [ "radius", "Ang",  20, [0, inf], "volume",
      "Cylinder radius" ],
    [ "length", "Ang",  400, [0, inf], "volume",
      "Cylinder length" ],
    [ "cyl_theta", "degrees", 60, [-inf, inf], "orientation",
      "In plane angle" ],
    [ "cyl_phi", "degrees", 60, [-inf, inf], "orientation",
      "Out of plane angle" ],
    ]
source = [ "lib/J1.c", "lib/gauss76.c", "cylinder_clone.c" ]

def ER(radius, length):
    ddd = 0.75*radius*(2*radius*length + (length+radius)*(length+pi*radius))
    return 0.5 * (ddd)**(1./3.)

