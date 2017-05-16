"""
Application to explore the difference between sasview 3.x orientation
dispersity and possible replacement algorithms.
"""
import mpl_toolkits.mplot3d   # Adds projection='3d' option to subplot
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
from matplotlib import cm

import numpy as np
from numpy import pi, cos, sin, sqrt, exp, degrees, radians

def draw_sphere(ax, radius=10., steps=100):
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)

    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')

def draw_mesh_current(ax, theta, dtheta, phi, dphi, radius=10., dist='gauss'):
    theta = radians(theta)
    phi = radians(phi)
    dtheta = radians(dtheta)
    dphi = radians(dphi)

    # 10 point 3-sigma gaussian weights
    t = np.linspace(-3., 3., 11)
    if dist == 'gauss':
        weights = exp(-0.5*t**2)
    elif dist == 'rect':
        weights = np.ones_like(t)
    else:
        raise ValueError("expected dist to be 'gauss' or 'rect'")
    theta = theta + dtheta*t
    phi = phi + dphi*t

    x = radius * np.outer(cos(phi), cos(theta))
    y = radius * np.outer(sin(phi), cos(theta))
    z = radius * np.outer(np.ones_like(phi), sin(theta))
    w = np.outer(weights, weights*abs(cos(theta)))

    x,y,z,w = [v.flatten() for v in (x,y,z,w)]

    ax.scatter(x, y, z, c=w, marker='o', vmin=0., vmax=1.0)

def draw_mesh_new(ax, theta, dtheta, phi, dphi, flow, radius=10., dist='gauss'):
    theta_center = radians(90-theta)
    phi_center = radians(phi)
    flow_center = radians(flow)
    dtheta = radians(dtheta)
    dphi = radians(dphi)

    # 10 point 3-sigma gaussian weights
    t = np.linspace(-3., 3., 11)
    if dist == 'gauss':
        weights = exp(-0.5*t**2)
    elif dist == 'rect':
        weights = np.ones_like(t)
    else:
        raise ValueError("expected dist to be 'gauss' or 'rect'")
    theta = dtheta*t
    phi = dphi*t

    x = radius * np.outer(cos(phi), cos(theta))
    y = radius * np.outer(sin(phi), cos(theta))
    z = radius * np.outer(np.ones_like(phi), sin(theta))
    #w = np.outer(weights, weights*abs(cos(dtheta*t)))
    w = np.outer(weights, weights*abs(cos(theta)))

    x, y, z, w = [v.flatten() for v in (x,y,z,w)]
    x, y, z = rotate(x, y, z, phi_center, theta_center, flow_center)

    ax.scatter(x, y, z, c=w, marker='o', vmin=0., vmax=1.)

def rotate(x, y, z, phi, theta, psi):
    R = rotation_matrix(psi, theta, phi)
    p = np.vstack([x,y,z])
    q = np.dot(R,p)
    return q

def rotation_matrix(xa,ya,za):
    Rz = [[cos(za), -sin(za), 0.],
          [sin(za),  cos(za), 0.],
          [0., 0., 1.]]
    Ry = [[cos(ya), 0., -sin(ya)],
          [0., 1., 0.],
          [sin(ya), 0.,  cos(ya)]]
    Rx = [[1., 0., 0.],
          [0.,  cos(xa), sin(xa)],
          [0., -sin(xa), cos(xa)]]
    R = np.dot(np.dot(Rz, Ry), Rx)
    return R

def main():
    plt.hold(True)
    plt.set_cmap('gist_earth')
    plt.clf()
    #gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    #ax = plt.subplot(gs[0], projection='3d')
    ax = plt.axes([0.0, 0.2, 1.0, 0.8], projection='3d')

    phi, dphi = -45., 3.
    theta, dtheta = 70., 10.
    flow = 0.
    #dist = 'rect'
    dist = 'gauss'

    axcolor = 'lightgoldenrodyellow'
    axphi = plt.axes([0.1, 0.1, 0.45, 0.04], axisbg=axcolor)
    axtheta  = plt.axes([0.1, 0.15, 0.45, 0.04], axisbg=axcolor)
    sphi = Slider(axphi, 'Phi', -180, 180, valinit=phi)
    stheta = Slider(axtheta, 'Theta', -180, 180, valinit=theta)
    axdphi = plt.axes([0.75, 0.1, 0.15, 0.04], axisbg=axcolor)
    axdtheta  = plt.axes([0.75, 0.15, 0.15, 0.04], axisbg=axcolor)
    sdphi = Slider(axdphi, 'dPhi', 0, 30, valinit=dphi)
    sdtheta = Slider(axdtheta, 'dTheta', 0, 30, valinit=dtheta)

    axflow = plt.axes([0.1, 0.05, 0.45, 0.04], axisbg=axcolor)
    sflow = Slider(axflow, 'Flow', -180, 180, valinit=flow)
    axusenew= plt.axes([0.75, 0.05, 0.15, 0.04], axisbg=axcolor)
    susenew = CheckButtons(axusenew, ['New'], [True])

    def update(val, axis=None):
        phi, theta = sphi.val, stheta.val
        dphi, dtheta = sdphi.val, sdtheta.val
        flow = sflow.val
        use_new = susenew.lines[0][0].get_visible()
        ax.cla()
        draw_sphere(ax)
        if use_new:
            draw_mesh_new(ax, theta=theta, dtheta=dtheta, phi=phi, dphi=dphi,
                          flow=flow, radius=11., dist=dist)
        else:
            draw_mesh_current(ax, theta=theta, dtheta=dtheta, phi=phi, dphi=dphi,
                              radius=11., dist=dist)
        if not axis.startswith('d'):
            ax.view_init(elev=90-theta if use_new else theta, azim=phi)
        plt.gcf().canvas.draw()

    stheta.on_changed(lambda v: update(v,'theta'))
    sphi.on_changed(lambda v: update(v, 'phi'))
    sdtheta.on_changed(lambda v: update(v, 'dtheta'))
    sdphi.on_changed(lambda v: update(v, 'dphi'))
    sflow.on_changed(lambda v: update(v, 'dflow'))
    susenew.on_clicked(lambda v: update(v, 'use_new'))

    update(None, 'phi')

    plt.show()

if __name__ == "__main__":
    main()