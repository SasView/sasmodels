r"""
Calculates the macroscopic scattering intensity for a multi-component
homogeneous mixture of polymers using the Random Phase Approximation.
This general formalism contains 10 specific cases

Case 0: C/D binary mixture of homopolymers

Case 1: C-D diblock copolymer

Case 2: B/C/D ternary mixture of homopolymers

Case 3: C/C-D mixture of a homopolymer B and a diblock copolymer C-D

Case 4: B-C-D triblock copolymer

Case 5: A/B/C/D quaternary mixture of homopolymers

Case 6: A/B/C-D mixture of two homopolymers A/B and a diblock C-D

Case 7: A/B-C-D mixture of a homopolymer A and a triblock B-C-D

Case 8: A-B/C-D mixture of two diblock copolymers A-B and C-D

Case 9: A-B-C-D tetra-block copolymer

**NB: these case numbers are different from those in the NIST SANS package!**

Only one case can be used at any one time.

The RPA (mean field) formalism only applies only when the multicomponent
polymer mixture is in the homogeneous mixed-phase region.

**Component D is assumed to be the "background" component (ie, all contrasts
are calculated with respect to component D).** So the scattering contrast
for a C/D blend = [SLD(component C) - SLD(component D)]\ :sup:`2`.

Depending on which case is being used, the number of fitting parameters - the
segment lengths (ba, bb, etc) and $\chi$ parameters (Kab, Kac, etc) - vary.
The *scale* parameter should be held equal to unity.

The input parameters are the degrees of polymerization, the volume fractions,
the specific volumes, and the neutron scattering length densities for each
component.


References
----------

A Z Akcasu, R Klein and B Hammouda, *Macromolecules*, 26 (1993) 4136
"""

from numpy import inf

name = "rpa"
title = "Random Phase Approximation - unfinished work in progress"
description = """
This formalism applies to multicomponent polymer mixtures in the
homogeneous (mixed) phase region only.
Case 0: C/D binary mixture of homopolymers
Case 1: C-D diblock copolymer
Case 2: B/C/D ternary mixture of homopolymers
Case 3: B/C-D mixture of homopolymer b and diblock copolymer C-D
Case 4: B-C-D triblock copolymer
Case 5: A/B/C/D quaternary mixture of homopolymers
Case 6: A/B/C-D mixture of two homopolymers A/B and a diblock C-D
Case 7: A/B-C-D mixture of a homopolymer A and a triblock B-C-D
Case 8: A-B/C-D mixture of two diblock copolymers A-B and C-D
Case 9: A-B-C-D four-block copolymer
See details in the model function help
"""
category = ""

CASES = [
    "C+D binary mixture",
    "C:D diblock copolymer",
    "B+C+D ternary mixture",
    "B+C:D binary mixture",
    "B:C:D triblock copolymer",
    "A+B+C+D quaternary mixture",
    "A+B+C:D ternary mixture",
    "A+B:C:D binary mixture",
    "A:B+C:D binary mixture",
    "A:B:C:D quadblock copolymer",
]

#   ["name", "units", default, [lower, upper], "type","description"],
parameters = [
    ["case_num", "", 1, [CASES], "", "Component organization"],

    ["N[4]", "", 1000.0, [1, inf], "", "Degree of polymerization"],
    ["Phi[4]", "", 0.25, [0, 1], "", "volume fraction"],
    ["v[4]", "mL/mol", 100.0, [0, inf], "", "specific volume"],
    ["L[4]", "fm", 10.0, [-inf, inf], "", "scattering length"],
    ["b[4]", "Ang", 5.0, [0, inf], "", "segment length"],

    ["Kab", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
    ["Kac", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
    ["Kad", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
    ["Kbc", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
    ["Kbd", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
    ["Kcd", "", -0.0004, [-inf, inf], "", "Interaction parameter"],
]

category = "shape-independent"

source = ["rpa.c"]

control = "case_num"
HIDE_NONE = set()
HIDE_A = set("Na Phia va La Kab Kac Kad".split())
HIDE_AB = set("Nb Phib vb Lb Kbc Kbd".split()).union(HIDE_A)
def hidden(pars):
    case_num = pars.get("case_num", parameters[0][2])
    if case_num < 2:
        return HIDE_AB
    elif case_num < 5:
        return HIDE_A
    else:
        return HIDE_NONE

