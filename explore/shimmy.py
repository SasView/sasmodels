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

def draw_beam(ax):
    #ax.plot([0,0],[0,0],[1,-1])
    #ax.scatter([0]*100,[0]*100,np.linspace(1, -1, 100), alpha=0.8)

    steps = 25
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(-1, 1, steps)

    r = 0.02
    x = r*np.outer(np.cos(u), np.ones_like(v))
    y = r*np.outer(np.sin(u), np.ones_like(v))
    z = np.outer(np.ones_like(u), v)

    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='y', alpha=0.5)
    
def draw_shimmy(ax, theta, phi, psi, dtheta, dphi, dpsi):
    size=[0.1, 0.4, 1.0] 
    view=[theta, phi, psi]
    shimmy=[0,0,0]
    #draw_shape = draw_parallelepiped
    draw_shape = draw_ellipsoid
    
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
    if dtheta == 0:
        cloud = [v for v in cloud if v[0] == 0]
    if dphi == 0:
        cloud = [v for v in cloud if v[1] == 0]
    if dpsi == 0:
        cloud = [v for v in cloud if v[2] == 0]
    draw_shape(ax, size, view, shimmy, steps=100, alpha=0.8)
    for point in cloud:
        shimmy=[dtheta*point[0], dphi*point[1], dpsi*point[2]]
        draw_shape(ax, size, view, shimmy, alpha=0.8)
    for v in 'xyz':
        a, b, c = size
        lim = np.sqrt(a**2+b**2+c**2)
        getattr(ax, 'set_'+v+'lim')([-lim, lim])
        getattr(ax, v+'axis').label.set_text(v)

def draw_ellipsoid(ax, size, view, shimmy, steps=25, alpha=1):
    a,b,c = size
    theta, phi, psi = view
    dtheta, dphi, dpsi = shimmy

    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)
    x = a*np.outer(np.cos(u), np.sin(v))
    y = b*np.outer(np.sin(u), np.sin(v))
    z = c*np.outer(np.ones_like(u), np.cos(v))

    shape = x.shape
    points = np.matrix([x.flatten(),y.flatten(),z.flatten()])
    points = Rz(dpsi)*Ry(dtheta)*Rx(dphi)*points
    points = Rz(phi)*Ry(theta)*Rz(psi)*points
    x,y,z = [v.reshape(shape) for v in points]

    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w', alpha=alpha)

def draw_parallelepiped(ax, size, view, shimmy, alpha=1):
    a,b,c = size
    theta, phi, psi = view
    dtheta, dphi, dpsi = shimmy

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

    points = np.matrix([x,y,z])
    points = Rz(dpsi)*Ry(dtheta)*Rx(dphi)*points
    points = Rz(phi)*Ry(theta)*Rz(psi)*points
	
    x,y,z = [np.array(v).flatten() for v in points]
    ax.plot_trisurf(x, y, triangles=tri, Z=z, color='w', alpha=alpha)

def draw_sphere(ax, radius=10., steps=100):
    u = np.linspace(0, 2 * np.pi, steps)
    v = np.linspace(0, np.pi, steps)

    x = radius * np.outer(np.cos(u), np.sin(v))
    y = radius * np.outer(np.sin(u), np.sin(v))
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')

def draw_mesh_new(ax, theta, dtheta, phi, dphi, flow, radius=10., dist='gauss'):
    theta_center = radians(theta)
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
    R = Rz(phi)*Ry(theta)*Rz(psi)
    p = np.vstack([x,y,z])
    return R*p

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
    sdpsi = Slider(axdpsi, 'dPsi', 0, 30, valinit=dphi)

    def update(val, axis=None):
        theta, phi, psi = stheta.val, sphi.val, spsi.val
        dtheta, dphi, dpsi = sdtheta.val, sdphi.val, sdpsi.val
        ax.cla()
        draw_beam(ax)
        draw_shimmy(ax, theta, phi, psi, dtheta, dphi, dpsi)
        #if not axis.startswith('d'):
        #    ax.view_init(elev=theta, azim=phi)
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