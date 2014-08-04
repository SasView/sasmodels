#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from Models.code_coreshellcyl import GpuCoreShellCylinder
from sasmodel import SasModel, load_data, set_beam_stop, set_top


radial_data = load_data('December/DEC07267.DAT')
set_beam_stop(radial_data, 0.00669, outer=0.025)
set_top(radial_data, -.0185)

tan_data = load_data('December/DEC07266.DAT')
set_beam_stop(tan_data, 0.00669, outer=0.025)
set_top(tan_data, -.0185)




dtype='float32'
radial = SasModel(radial_data,
                 GpuCoreShellCylinder,
                 scale= 3.75e-7, radius=378, thickness=30, length=1806,
                 core_sld=7.105e-6, shell_sld=.291e-6, solvent_sld=7.105e-6,
                 background=0.2, axis_theta=0, axis_phi=90,

                 radius_pd=0.26, radius_pd_n=20, radius_pd_nsigma=3,
                 length_pd=0.26, length_pd_n=20, length_pd_nsigma=3,
                 thickness_pd=1, thickness_pd_n=1, thickness_pd_nsigma=0,
                 axis_theta_pd=1, axis_theta_pd_n=10, axis_theta_pd_nsigma=3,
                 axis_phi_pd=0.1, axis_phi_pd_n=1, axis_phi_pd_nsigma=0,
                 dtype='float')
tan = SasModel(tan_data,
                  GpuCoreShellCylinder, dtype=dtype,
                  **radial.parameters())

radial.radius.range(15, 1000)
radial.length.range(0, 2500)
#radial.thickness.range(18, 38)
#radial.thickness_pd.range(0, 1)
#radial.axis_phi.range(0, 90)
#radial.radius_pd.range(0, 1)
#radial.length_pd.range(0, 1)
#radial.axis_theta_pd.range(0, 360)
#radial.background.range(0,5)
#radial.scale.range(0, 1)
radial.axis_phi = tan.axis_phi + 90


problem = FitProblem([radial,tan])



