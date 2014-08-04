#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from Models.code_lamellar import GpuLamellar
from Models.code_ellipse import GpuEllipse
from Models.code_cylinder import GpuCylinder
from multisasmodels import SasModel, load_data, set_beam_stop, set_half
import numpy as np

data = load_data('December/DEC07238.DAT')
set_beam_stop(data, 0.0052)#, outer=0.025)
#set_half(data, 'left')

truth = SasModel(data, GpuCylinder,
scale=0.08, radius=46, length=719,
sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
cyl_theta=90, cyl_phi=0,
cyl_theta_pd=23, cyl_theta_pd_n=40, cyl_theta_pd_nsigma=3,
radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
length_pd=0.1, length_pd_n=10, length_pd_nsigma=3,
cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3,
dtype='float')


#model.radius_a.range(15, 1000)
#model.radius_b.range(15, 1000)
#model.axis_theta_pd.range(0, 360)
#model.background.range(0,1000)
#truth.scale.range(0, 1)

arrayone = truth.theory()


lies = SasModel(data, GpuCylinder,
scale=0.08, radius=46, length=719,
sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
cyl_theta=90, cyl_phi=0,
cyl_theta_pd=23, cyl_theta_pd_n=40, cyl_theta_pd_nsigma=3,
radius_pd=0.1, radius_pd_n=10, radius_pd_nsigma=3,
length_pd=0.1, length_pd_n=10, length_pd_nsigma=3,
cyl_phi_pd=0.1, cyl_phi_pd_n=10, cyl_phi_pd_nsigma=3,
dtype='float')



#modeltwo.bi_thick.range(0, 1000)
#modeltwo.scale.range(0, 1)
#model.bi_thick_pd.range(0, 1000)
#model.background.range(0, 1000)


arraytwo = lies.theory()


a = np.add(np.multiply(arrayone, .08), np.multiply(arraytwo, .08))
truth.set_result(a)


problem = FitProblem(truth)

