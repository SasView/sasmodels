#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sans.dataloader.loader import Loader
from sans.dataloader.manipulations import Ringcut
from Kernel import fit


def cylinder(data): ##################################################################################
    from sans.models.CylinderModel import CylinderModel
    model1 = CylinderModel()
    model1.setParam("scale", 10.0)
    model1.setParam("radius",18)
    model1.setParam("length", 397)
    model1.setParam("sldCyl",3e-006 )
    model1.setParam("sldSolv",0.0 )
    model1.setParam("background", 0.0)

    # Dispersion parameters
    model1.dispersion['radius']['width'] = 0.25
    model1.dispersion['radius']['npts'] = 50

    theory = model1.evalDistribution([data.qx_data, data.qy_data])
    return theory

def load_data(filename):
    loader = Loader()
    f = loader.load(filename)
    return f

def set_beam_stop(data, radius):
    data.mask = Ringcut(0, radius)(data)

def plot_data(data):
    from numpy.ma import masked_array
    import matplotlib.pyplot as plt
    img = masked_array(data.data, data.mask)
    xmin, xmax = min(data.qx_data), max(data.qx_data)
    ymin, ymax = min(data.qy_data), max(data.qy_data)
    plt.imshow(img.reshape(128,128),
               interpolation='nearest', aspect=1, origin='upper',
               extent=[xmin, xmax, ymin, ymax])

def demo():

    data = load_data('NOV07090.DAT')
    """
    print data
    print type(data)
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()
    """
    import matplotlib.pyplot as plt
    import numpy as np
    sasviewcyl = cylinder(data)
    gpucyl = fit.Fit(data)
    diff = gpucyl - sasviewcyl
    np.linalg.norm(diff, order=None)
    np.linalg.norm(diff, order=np.inf)

    plot_data(data, diff)
    plt.show()




if __name__ == "__main__":
    demo()

