#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from cylcode import GpuCylinder
from lamellarcode import GpuLamellar
from ellipsecode import GpuEllipse
from coreshellcylcode import GpuCoreShellCylinder
from sasmodel import SasModel, load_data, set_beam_stop


radial_data = load_data('JUN03289.DAT')
set_beam_stop(radial_data, 0.004)
trans_data = load_data('JUN03305.DAT')
set_beam_stop(trans_data, 0.004)




dtype='float32'
radial = SasModel(radial_data,
                  GpuCylinder, dtype=dtype,
                  scale=1, radius=64.1, length=266.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
                  cyl_theta=0, cyl_phi=0, radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,length_pd=0.1,
                  length_pd_n=5, length_pd_nsigma=3, cyl_theta_pd=0.1, cyl_theta_pd_n=5, cyl_theta_pd_nsigma=3,
                  cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3)
trans = SasModel(trans_data,
                  GpuCylinder, dtype=dtype,
                  **radial.parameters())

radial.radius.range(0,100)
radial.length.range(0, 1000)
radial.cyl_theta.range(0,90)
radial.cyl_phi.range(0,90)
radial.scale.range(0,10)
trans.cyl_theta = radial.cyl_theta + 90.


problem = FitProblem([radial,trans])



