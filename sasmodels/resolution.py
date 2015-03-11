from scipy.special import erf
from numpy import sqrt
import numpy as np

SLIT_SMEAR_POINTS = 500

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
    weights /= sum(weights, axis=1)
    return weights

def slit_resolution(q_calc, q, qx_width, qy_width):
    edges = bin_edges(q_calc) # Note: requires q > 0
    edges[edges<0.] = 0.0 # clip edges below zero
    qy_min, qy_max = 0.0, edges[-1]

    weights = np.zeros((len(q),len(q_calc)),'d')
    # Loop for width (height is analytical).
    # Condition: height >>> width, otherwise, below is not accurate enough.
    # Smear weight numerical iteration for width>0 when height>0.
    # When width = 0, the numerical iteration will be skipped.
    # The resolution calculation for the height is done by direct integration,
    # assuming the I(q'=sqrt(q_j^2-(q+shift_w)^2)) is constant within
    # a q' bin, [q_high, q_low].
    # In general, this weight numerical iteration for width>0 might be a rough
    # approximation, but it must be good enough when height >>> width.
    E_sq = edges**2[:,None]
    y_points = SLIT_SMEAR_POINTS if np.any(qy_width>0) else 1
    qy_step = 0 if y_points == 1 else qy_width/(y_points-1)
    for k in range(-y_points+1,y_points):
        qy = np.clip(q + qy_step*k, qy_min, qy_max)
        qx_low = qy
        qx_high = sqrt(qx_low**2 + qx_width**2)
        in_x = (q_calc[:,None]>=qx_low[None,:])*(q_calc[:,None]<=qx_high[None,:])
        qy_sq = qy**2[None,:]
        weights += (sqrt(E_sq[1:]-qy_sq) - sqrt(E_sq[:-1]-qy_sq))*in_x
    weights /= sum(weights, axis=1)
    return weights

def bin_edges(x):
    if len(x) < 2 or (np.diff(x)<0).any():
        raise ValueError("Expected bins to be an increasing set")
    edges = np.hstack([
        x[0]  - 0.5*(x[1]  - x[0]),  # first point minus half first interval
        0.5*(x[1:] + x[:-1]),        # mid points of all central intervals
        x[-1] + 0.5*(x[-1] - x[-2]), # last point plus half last interval
        ])
    return edges
