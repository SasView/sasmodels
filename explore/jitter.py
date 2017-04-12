"""
Application to explore the difference between sasview 3.x orientation
dispersity and possible replacement algorithms.
"""
import sys

import mpl_toolkits.mplot3d   # Adds projection='3d' option to subplot
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
from matplotlib import cm
import numpy as np
from numpy import pi, cos, sin, sqrt, exp, degrees, radians

def draw_beam(ax, view=(0, 0)):
    #ax.plot([0,0],[0,0],[1,-1])
    #ax.scatter([0]*100,[0]*100,np.linspace(1, -1, 100), alpha=0.8)

    steps = 25
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(-1, 1, steps)

    r = 0.02
    x = r*np.outer(np.cos(u), np.ones_like(v))
    y = r*np.outer(np.sin(u), np.ones_like(v))
    z = 1.3*np.outer(np.ones_like(u), v)

    theta, phi = view
    shape = x.shape
    points = np.matrix([x.flatten(), y.flatten(), z.flatten()])
    points = Rz(phi)*Ry(theta)*points
    x, y, z = [v.reshape(shape) for v in points]

    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='y', alpha=0.5)

def draw_jitter(ax, view, jitter):
    size = [0.1, 0.4, 1.0]
    draw_shape = draw_parallelepiped
    #draw_shape = draw_ellipsoid

    #np.random.seed(10)
    #cloud = np.random.randn(10,3)
    cloud = [
        [-1, -1, -1],
        [-1, -1,  0],
        [-1, -1,  1],
        [-1,  0, -1],
        [-1,  0,  0],
        [-1,  0,  1],
        [-1,  1, -1],
        [-1,  1,  0],
        [-1,  1,  1],
        [ 0, -1, -1],
        [ 0, -1,  0],
        [ 0, -1,  1],
        [ 0,  0, -1],
        [ 0,  0,  0],
        [ 0,  0,  1],
        [ 0,  1, -1],
        [ 0,  1,  0],
        [ 0,  1,  1],
        [ 1, -1, -1],
        [ 1, -1,  0],
        [ 1, -1,  1],
        [ 1,  0, -1],
        [ 1,  0,  0],
        [ 1,  0,  1],
        [ 1,  1, -1],
        [ 1,  1,  0],
        [ 1,  1,  1],
    ]
    dtheta, dphi, dpsi = jitter
    if dtheta == 0:
        cloud = [v for v in cloud if v[0] == 0]
    if dphi == 0:
        cloud = [v for v in cloud if v[1] == 0]
    if dpsi == 0:
        cloud = [v for v in cloud if v[2] == 0]
    draw_shape(ax, size, view, [0, 0, 0], steps=100, alpha=0.8)
    for point in cloud:
        delta = [dtheta*point[0], dphi*point[1], dpsi*point[2]]
        draw_shape(ax, size, view, delta, alpha=0.8)
    for v in 'xyz':
        a, b, c = size
        lim = np.sqrt(a**2+b**2+c**2)
        getattr(ax, 'set_'+v+'lim')([-lim, lim])
        getattr(ax, v+'axis').label.set_text(v)

def draw_ellipsoid(ax, size, view, jitter, steps=25, alpha=1):
    a,b,c = size
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)
    x = a*np.outer(np.cos(u), np.sin(v))
    y = b*np.outer(np.sin(u), np.sin(v))
    z = c*np.outer(np.ones_like(u), np.cos(v))
    x, y, z = transform_xyz(view, jitter, x, y, z)

    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w', alpha=alpha)

    draw_labels(ax, view, jitter, [
         ('c+', [ 0, 0, c], [ 1, 0, 0]),
         ('c-', [ 0, 0,-c], [ 0, 0,-1]),
         ('a+', [ a, 0, 0], [ 0, 0, 1]),
         ('a-', [-a, 0, 0], [ 0, 0,-1]),
         ('b+', [ 0, b, 0], [-1, 0, 0]),
         ('b-', [ 0,-b, 0], [-1, 0, 0]),
    ])

def draw_parallelepiped(ax, size, view, jitter, steps=None, alpha=1):
    a,b,c = size
    x = a*np.array([ 1,-1, 1,-1, 1,-1, 1,-1])
    y = b*np.array([ 1, 1,-1,-1, 1, 1,-1,-1])
    z = c*np.array([ 1, 1, 1, 1,-1,-1,-1,-1])
    tri = np.array([
        # counter clockwise triangles
        # z: up/down, x: right/left, y: front/back
        [0,1,2], [3,2,1], # top face
        [6,5,4], [5,6,7], # bottom face
        [0,2,6], [6,4,0], # right face
        [1,5,7], [7,3,1], # left face
        [2,3,6], [7,6,3], # front face
        [4,1,0], [5,1,4], # back face
    ])

    x, y, z = transform_xyz(view, jitter, x, y, z)
    ax.plot_trisurf(x, y, triangles=tri, Z=z, color='w', alpha=alpha)

    draw_labels(ax, view, jitter, [
         ('c+', [ 0, 0, c], [ 1, 0, 0]),
         ('c-', [ 0, 0,-c], [ 0, 0,-1]),
         ('a+', [ a, 0, 0], [ 0, 0, 1]),
         ('a-', [-a, 0, 0], [ 0, 0,-1]),
         ('b+', [ 0, b, 0], [-1, 0, 0]),
         ('b-', [ 0,-b, 0], [-1, 0, 0]),
    ])

def draw_mesh(ax, view, jitter, radius=1.2, n=11, dist='gauss'):
    theta, phi, psi = view
    dtheta, dphi, dpsi = jitter
    if dist == 'gauss':
        t = np.linspace(-3, 3, n)
        weights = exp(-0.5*t**2)
    elif dist == 'rect':
        t = np.linspace(0, 1, n)
        weights = np.ones_like(t)
    else:
        raise ValueError("expected dist to be 'gauss' or 'rect'")

    # mesh in theta, phi formed by rotating z
    z = np.matrix([[0], [0], [radius]])
    points = np.hstack([Rx(phi_i)*Ry(theta_i)*z
                        for theta_i in dtheta*t
                        for phi_i in dphi*t])
    # rotate relative to beam
    points = orient_relative_to_beam(view, points)

    w = np.outer(weights, weights)

    x, y, z = [np.array(v).flatten() for v in points]
    ax.scatter(x, y, z, c=w.flatten(), marker='o', vmin=0., vmax=1.)

def Rx(angle):
    a = radians(angle)
    R = [[1., 0., 0.],
         [0.,  cos(a), sin(a)],
         [0., -sin(a), cos(a)]]
    return np.matrix(R)

def Ry(angle):
    a = radians(angle)
    R = [[cos(a), 0., -sin(a)],
         [0., 1., 0.],
         [sin(a), 0.,  cos(a)]]
    return np.matrix(R)

def Rz(angle):
    a = radians(angle)
    R = [[cos(a), -sin(a), 0.],
         [sin(a),  cos(a), 0.],
         [0., 0., 1.]]
    return np.matrix(R)

def transform_xyz(view, jitter, x, y, z):
    x, y, z = [np.asarray(v) for v in (x, y, z)]
    shape = x.shape
    points = np.matrix([x.flatten(),y.flatten(),z.flatten()])
    points = apply_jitter(jitter, points)
    points = orient_relative_to_beam(view, points)
    x, y, z = [np.array(v).reshape(shape) for v in points]
    return x, y, z

def apply_jitter(jitter, points):
    dtheta, dphi, dpsi = jitter
    points = Rx(dphi)*Ry(dtheta)*Rz(dpsi)*points
    return points

def orient_relative_to_beam(view, points):
    theta, phi, psi = view
    points = Rz(phi)*Ry(theta)*Rz(psi)*points
    return points

def draw_labels(ax, view, jitter, text):
    labels, locations, orientations = zip(*text)
    px, py, pz = zip(*locations)
    dx, dy, dz = zip(*orientations)

    px, py, pz = transform_xyz(view, jitter, px, py, pz)
    dx, dy, dz = transform_xyz(view, jitter, dx, dy, dz)

    for label, p, zdir in zip(labels, zip(px, py, pz), zip(dx, dy, dz)):
        zdir = np.asarray(zdir).flatten()
        ax.text(p[0], p[1], p[2], label, zdir=zdir)

def draw_sphere(ax, radius=10., steps=100):
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)

    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')

def main():
    #plt.hold(True)
    plt.set_cmap('gist_earth')
    plt.clf()
    #gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    #ax = plt.subplot(gs[0], projection='3d')
    ax = plt.axes([0.0, 0.2, 1.0, 0.8], projection='3d')

    theta, dtheta = 70., 10.
    phi, dphi = -45., 3.
    psi, dpsi = -45., 3.
    theta, phi, psi = 0, 0, 0
    dtheta, dphi, dpsi = 0, 0, 0
    #dist = 'rect'
    dist = 'gauss'

    axcolor = 'lightgoldenrodyellow'

    axtheta  = plt.axes([0.1, 0.15, 0.45, 0.04], axisbg=axcolor)
    axphi = plt.axes([0.1, 0.1, 0.45, 0.04], axisbg=axcolor)
    axpsi = plt.axes([0.1, 0.05, 0.45, 0.04], axisbg=axcolor)
    stheta = Slider(axtheta, 'Theta', -90, 90, valinit=theta)
    sphi = Slider(axphi, 'Phi', -180, 180, valinit=phi)
    spsi = Slider(axpsi, 'Psi', -180, 180, valinit=psi)

    axdtheta  = plt.axes([0.75, 0.15, 0.15, 0.04], axisbg=axcolor)
    axdphi = plt.axes([0.75, 0.1, 0.15, 0.04], axisbg=axcolor)
    axdpsi= plt.axes([0.75, 0.05, 0.15, 0.04], axisbg=axcolor)
    sdtheta = Slider(axdtheta, 'dTheta', 0, 30, valinit=dtheta)
    sdphi = Slider(axdphi, 'dPhi', 0, 30, valinit=dphi)
    sdpsi = Slider(axdpsi, 'dPsi', 0, 30, valinit=dpsi)

    def update(val, axis=None):
        view = stheta.val, sphi.val, spsi.val
        jitter = sdtheta.val, sdphi.val, sdpsi.val
        ax.cla()
        draw_beam(ax, (0, 0))
        draw_jitter(ax, view, jitter)
        #draw_jitter(ax, view, (0,0,0))
        draw_mesh(ax, view, jitter)
        plt.gcf().canvas.draw()

    stheta.on_changed(lambda v: update(v,'theta'))
    sphi.on_changed(lambda v: update(v, 'phi'))
    spsi.on_changed(lambda v: update(v, 'psi'))
    sdtheta.on_changed(lambda v: update(v, 'dtheta'))
    sdphi.on_changed(lambda v: update(v, 'dphi'))
    sdpsi.on_changed(lambda v: update(v, 'dpsi'))

    update(None, 'phi')

    plt.show()

if __name__ == "__main__":
    main()