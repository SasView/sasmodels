#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bumps.names import *
from code_lamellar import GpuLamellar
from code_ellipse import GpuEllipse
from sasmodel import SasModel, load_data, set_beam_stop, set_half
import numpy as np

data = load_data('DEC07269.DAT')
set_beam_stop(data, 0.0052)#, outer=0.025)
#set_half(data, 'left')

truth = SasModel(data, GpuEllipse,
                scale=.5, radius_a=45.265, radius_b=600.8, sldEll=.291e-6, sldSolv=7.105e-6,
                 background=8.30161, axis_theta=0, axis_phi=0,
                 radius_a_pd=0.222296, radius_a_pd_n=1, radius_a_pd_nsigma=0,
                 radius_b_pd=.000128, radius_b_pd_n=1, radius_b_pd_nsigma=0,
                 axis_theta_pd=20, axis_theta_pd_n=40, axis_theta_pd_nsigma=3,
                 axis_phi_pd=2.63698e-05, axis_phi_pd_n=20, axis_phi_pd_nsigma=0,
                 dtype='float')

#model.radius_a.range(15, 1000)
#model.radius_b.range(15, 1000)
#model.axis_theta_pd.range(0, 360)
#model.background.range(0,1000)
#truth.scale.range(0, 1)

arrayone = truth.theory()


lies = SasModel(data, GpuLamellar, scale=0, bi_thick=5, sld_bi=.291e-6, sld_sol=5.77e-6, background=85.23,
                 bi_thick_pd= 0.0013, bi_thick_pd_n=5, bi_thick_pd_nsigma=3, dtype='float')


#modeltwo.bi_thick.range(0, 1000)
#modeltwo.scale.range(0, 1)
#model.bi_thick_pd.range(0, 1000)
#model.background.range(0, 1000)


arraytwo = lies.theory()


a = np.add(np.multiply(arrayone, .5), np.multiply(arraytwo, 0))
truth.set_result(a)


problem = FitProblem(truth)

