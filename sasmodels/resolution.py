from scipy.special import erf
from numpy import sqrt
import numpy as np

SLIT_SMEAR_POINTS = 500

class MatrixResolution:
    def apply(self, Iq):
        return np.dot(Iq, self.resolution_matrix)

class Pinhole1D(MatrixResolution):
    def __init__(self, q, q_width):
        self.q, self.q_width = q, q_width
        self.q_calc = pinhole_extend_q(q,q_width)
        self.resolution_matrix = \
            pinhole_resolution(self.q_calc, self.q, self.q_width)

class Slit1D(MatrixResolution):
    def __init__(self, q, qx_width, qy_width):
        if np.isscalar(qx_width):
            qx_width = qx_width*np.ones(len(q))
        if np.isscalar(qy_width):
            qy_width = qy_width*np.ones(len(q))
        self.q, self.qx_width, self.qy_width = q, qx_width, qy_width
        self.q_calc = slit_extend_q(q, qx_width, qy_width)
        self.resolution_matrix = \
            slit_resolution(self.q_calc, self.q, self.qx_width, self.qy_width)

def pinhole_resolution(q_calc, q, q_width):
    """
    Compute the convolution matrix *W* for pinhole resolution 1-D data.

    Each row *W[i]* determines the normalized weight that the corresponding
    points *q_calc* contribute to the resolution smeared point *q[i]*.  Given
    *W*, the resolution smearing can be computed using *dot(W,q)*.

    *q_calc* must be increasing.
    """
    edges = bin_edges(q_calc)
    edges[edges<0.] = 0. # clip edges below zero
    G = erf( (edges[:,None] - q[None,:]) / (sqrt(2.0)*q_width)[None,:] )
    weights = G[1:] - G[:-1]
    weights /= np.sum(weights, axis=1)
    return weights


def pinhole_extend_q(q, q_width):
    return q


def slit_resolution(q_calc, q, qx_width, qy_width):
    edges = bin_edges(q_calc) # Note: requires q > 0
    edges[edges<0.] = 0.0 # clip edges below zero
    qy_min, qy_max = 0.0, edges[-1]

    # Make q_calc into a row vector, and q, qx_width, qy_width into columns
    # Make weights match [ q_calc X q ]
    weights = np.zeros((len(q),len(q_calc)),'d')
    q_calc = q_calc[None,:]
    q, qx_width, qy_width, edges = [
        v[:,None] for v in (q, qx_width, qy_width, edges)]

    # Loop for width (height is analytical).
    # Condition: height >>> width, otherwise, below is not accurate enough.
    # Smear weight numerical iteration for width>0 when height>0.
    # When width = 0, the numerical iteration will be skipped.
    # The resolution calculation for the height is done by direct integration,
    # assuming the I(q'=sqrt(q_j^2-(q+shift_w)^2)) is constant within
    # a q' bin, [q_high, q_low].
    # In general, this weight numerical iteration for width>0 might be a rough
    # approximation, but it must be good enough when height >>> width.
    E_sq = edges**2
    y_points = SLIT_SMEAR_POINTS if np.any(qy_width>0) else 1
    qy_step = 0 if y_points == 1 else qy_width/(y_points-1)
    for k in range(-y_points+1,y_points):
        qy = np.clip(q + qy_step*k, qy_min, qy_max)
        qx_low = qy
        qx_high = sqrt(qx_low**2 + qx_width**2)
        in_x = (q_calc>=qx_low)*(q_calc<=qx_high)
        qy_sq = qy**2
        weights += (sqrt(E_sq[1:]-qy_sq) - sqrt(qy_sq - E_sq[:-1]))*in_x
    weights /= np.sum(weights, axis=1)
    return weights


def slit_extend_q(q, qx_width, qy_width):
    return q


def bin_edges(x):
    if len(x) < 2 or (np.diff(x)<0).any():
        raise ValueError("Expected bins to be an increasing set")
    edges = np.hstack([
        x[0]  - 0.5*(x[1]  - x[0]),  # first point minus half first interval
        0.5*(x[1:] + x[:-1]),        # mid points of all central intervals
        x[-1] + 0.5*(x[-1] - x[-2]), # last point plus half last interval
        ])
    return edges


############################################################################
# usage demo
############################################################################

def _eval_demo_1d(resolution, title):
    from sasmodels import core
    from sasmodels.models import cylinder
    ## or alternatively:
    # cylinder = core.load_model_definition('cylinder')
    model = core.load_model(cylinder)

    kernel = core.make_kernel(model, [resolution.q_calc])
    Iq_calc = core.call_kernel(kernel, {'length':210, 'radius':500})
    Iq = resolution.apply(Iq_calc)

    import matplotlib.pyplot as plt
    plt.loglog(resolution.q_calc, Iq_calc, label='unsmeared')
    plt.loglog(resolution.q, Iq, label='smeared', hold=True)
    plt.legend()
    plt.title(title)
    plt.xlabel("Q (1/Ang)")
    plt.ylabel("I(Q) (1/cm)")

def demo_pinhole_1d():
    q = np.logspace(-3,-1,400)
    dq = 0.1*q
    resolution = Pinhole1D(q, dq)
    _eval_demo_1d(resolution, title="10% dQ/Q Pinhole Resolution")

def demo_slit_1d():
    q = np.logspace(-3,-1,400)
    qx_width = 0.005
    qy_width = 0.0
    resolution = Slit1D(q, qx_width, qy_width)
    _eval_demo_1d(resolution, title="0.005 Qx Slit Resolution")

def demo():
    import matplotlib.pyplot as plt
    plt.subplot(121)
    demo_pinhole_1d()
    plt.subplot(122)
    demo_slit_1d()
    plt.show()


if __name__ == "__main__":
    demo()


