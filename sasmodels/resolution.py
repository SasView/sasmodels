from scipy.special import erf
from numpy import sqrt
import numpy as np

def pinhole_resolution(q_calc, q, dq):
    """
    Compute the convolution matrix *W* for pinhole resolution 1-D data.

    Each row *W[i]* determines the normalized weight that the corresponding
    points *q_calc* contribute to the resolution smeared point *q[i]*.  Given
    *W*, the resolution smearing can be computed using *dot(W,q)*.

    *q_calc* must be increasing.
    """
    edges = bin_edges(q_calc)
    edges[edges<0.] = 0. # clip edges below zero
    G = erf( (edges[:,None]-q[None,:]) / (sqrt(2.0)*dq)[None,:] )
    weights = G[1:,:] - G[:-1,:]
    weights /= sum(weights, axis=1)
    return weights

def slit_resolution(q_calc, q, qx_width, qy_width):
    edges = bin_edges(q_calc) # Note: requires q > 0
    edges[edges<0.] = 0. # clip edges below zero

    weights = np.zeros((len(q),len(q_calc)),'d')
    # Loop for width (;Height is analytical.)
    # Condition: height >>> width, otherwise, below is not accurate enough.
    # Smear weight numerical iteration for width >0 when the height (>0) presents.
    # When width = 0, the numerical iteration will be skipped.
    # The resolution calculation for the height is done by direct integration,
    # assuming the I(q'=sqrt(q_j^2-(q+shift_w)^2)) is constant within a q' bin, [q_high, q_low].
    # In general, this weight numerical iteration for width >0 might be a rough approximation,
    # but it must be good enough when height >>> width.
    E_sq = edges**2[:,None]
    y_pts = 500 if np.any(qy_width>0) else 1
    for k in range(-y_pts+1,y_pts):
        qy = q if y_pts == 1 else q + qy_width/(y_pts-1)*k
        qy = np.clip(qy, 0.0, edges[-1])
        qx_low = qy
        qx_high = sqrt(qx_low**2 + qx_width**2)
        in_x = (q_calc[:,None]>=qx_low[None,:])*(q_calc[:,None]<=qx_high[None,:])
        weights += (sqrt(E_sq[1:]-qy[None,:]**2)-sqrt(E_sq[:-1]-qy[None,:]**2))*in_x
    weights /= sum(weights, axis=1)
            # Condition: zero slit smear.
            if (npts_w == 1 and npts_h == 1):
                if(q_j == q) :
                    weights[i,j] = 1.0
            #Condition:Smear weight integration for width >0 when the height (=0) does not present.
            #Or height << width.
            elif (npts_w!=1 and npts_h==1)or(npts_w!=1 and npts_h != 1 and width/height > 100.0):
                shift_w = width
                #del_w = width/((double)npts_w-1.0);
                q_shifted_low = q - shift_w
                # High limit of the resolution range
                q_shifted_high = q + shift_w
                # Go through all the q_js for weighting those points
                if(q_j >= q_shifted_low and q_j <= q_shifted_high):
                    # The weighting factor comes,
                    # Give some weight (delq_bin) for the q_j within the resolution range
                    # Weight should be same for all qs except
                    # for the q bin size at j.
                    # Note that the division by q_0 is only due to the precision problem
                    # where q_high - q_low gets to very small.
                    # Later, it will be normalized again.
                    weights[i,j] += (q_high - q_low)/q_0
            else:
                # Loop for width (;Height is analytical.)
                # Condition: height >>> width, otherwise, below is not accurate enough.
                # Smear weight numerical iteration for width >0 when the height (>0) presents.
                # When width = 0, the numerical iteration will be skipped.
                # The resolution calculation for the height is done by direct integration,
                # assuming the I(q'=sqrt(q_j^2-(q+shift_w)^2)) is constant within a q' bin, [q_high, q_low].
                # In general, this weight numerical iteration for width >0 might be a rough approximation,
                # but it must be good enough when height >>> width.
                for k in range(-npts_w + 1,npts_w+1):
                    if(npts_w!=1):
                        shift_w = width/(npts_w-1.0)*k
                    # For each q-value, compute the weight of each other q-bin
                    # in the I(q) array
                    # Low limit of the resolution range
                    q_shift = q + shift_w
                    if (q_shift < 0.0):
                        q_shift = 0.0
                    q_shifted_low = q_shift
                    # High limit of the resolution range
                    q_shifted_high = sqrt(q_shift * q_shift + shift_h * shift_h)


                    # Go through all the q_js for weighting those points
                    if(q_j >= q_shifted_low and q_j <= q_shifted_high) :
                        # The weighting factor comes,
                        # Give some weight (delq_bin) for the q_j within the resolution range
                        # Weight should be same for all qs except
                        # for the q bin size at j.
                        # Note that the division by q_0 is only due to the precision problem
                        # where q_high - q_low gets to very small.
                        # Later, it will be normalized again.

                        # The fabs below are not necessary but in case: the weight should never be imaginary.
                        # At the edge of each sub_width. weight += u(at q_high bin) - u(0), where u(0) = 0,
                        # and weighted by (2.0* npts_w -1.0)once for each q.
                        #if (q == q_j) {
                        if (q_low <= q_shift and q_high > q_shift) :
                            #if (k==0)
                            weights[i,j] += (sqrt(abs((q_high)*(q_high)-q_shift * q_shift)))/q_0# * (2.0*double(npts_w)-1.0);
                        # For the rest of sub_width. weight += u(at q_high bin) - u(at q_low bin)
                        else:# if (u > 0.0){
                            weights[i,j] += (sqrt(abs((q_high)*(q_high)- q_shift * q_shift))-sqrt(abs((q_low)*(q_low)- q_shift * q_shift)))/q_0


def bin_edges(x):
    if len(x) < 2 or (np.diff(x)<0).any():
        raise ValueError("Expected bins to be an increasing set")
    edges = np.hstack([
        x[0]  - 0.5*(x[1]  - x[0]),  # first point minus half first interval
        0.5*(x[1:] + x[:-1]),        # mid points of all central intervals
        x[-1] + 0.5*(x[-1] - x[-2]), # last point plus half last interval
        ])
    return edges
