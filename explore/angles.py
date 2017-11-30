#!/usr/bin/env python
"""
Generate code for orientation transforms using symbolic algebra.

To make it easier to generate correct transforms for oriented shapes, we
use the sympy symbolic alegbra package to do the matrix multiplication.
The transforms are displayed both using an ascii math notation, and as
C or python code which can be pasted directly into the kernel driver.

If ever we decide to change conventions, we simply need to adjust the
order and parameters to the rotation matrices.  For display we want to
use forward transforms for the mesh describing the shape, first applying
jitter, then adjusting the view.  For calculation we know the effective q
so we instead need to first unwind the view, using the inverse rotation,
then undo the jitter to get the q to calculate for the shape in its
canonical orientation.

Set *OUTPUT* to the type of code you want to see: ccode, python, math
or any combination.
"""

from __future__ import print_function

import codecs
import sys
import re

import sympy as sp
from sympy import pi, sqrt, sin, cos, Matrix, Eq

# Select output
OUTPUT = ""
OUTPUT = OUTPUT + "ccode"
#OUTPUT = OUTPUT + "python "
OUTPUT = OUTPUT + "math "
REUSE_SINCOS = True
QC_ONLY = True # show only what is needed for dqc in the symmetric case

# include unicode symbols in output, even if piping to a pager
if sys.version_info[0] < 3:
    sys.stdout = codecs.getwriter('utf8')(sys.stdout)
sp.init_printing(use_unicode=True)

def subs(s):
    """
    Transform sympy generated code to follow sasmodels naming conventions.
    """
    if REUSE_SINCOS:
        s = re.sub(r'(phi|psi|theta)\^\+', r'\1', s)  # jitter rep:  V^+ => V
    s = re.sub(r'([a-z]*)\^\+', r'd\1', s)  # jitter rep:  V^+ => dV
    s = re.sub(r'(cos|sin)\(([a-z]*)\)', r'\1_\2', s)  # cos(V) => cos_V
    s = re.sub(r'pow\(([a-z]*), 2\)', r'\1*\1', s)  # pow(V, 2) => V*V
    return s

def comment(s):
    r"""
    Add a comment to the generated code.  Use '\n' to separate lines.
    """
    if 'ccode' in OUTPUT:
        for line in s.split("\n"):
            print("// " + line if line else "")
    if 'python' in OUTPUT:
        for line in s.split("\n"):
            print("    ## " + line if line else "")

def vprint(var, vec, comment=None, post=None):
    """
    Generate assignment statements.

    *var* could be a single sympy symbol or a 1xN vector of symbols.

    *vec* could be a single sympy expression or a 1xN vector of expressions
    such as results from a matrix-vector multiplication.

    *comment* if present is added to the start of the block as documentation.
    """
    #for v, row in zip(var, vec): sp.pprint(Eq(v, row))
    desc = sp.pretty(Eq(var, vec), wrap_line=False)
    if not isinstance(var, Matrix):
        var, vec = [var], [vec]
    if 'ccode' in OUTPUT:
        if 'math' in OUTPUT:
            print("\n// " + comment if comment else "")
            print("/*")
            for line in desc.split("\n"):
                print(" * "+line)
            print(" *\n */")
        else:
            print("\n    // " + comment if comment else "")
        if post:
            print("    // " + post)
        for v, row in zip(var, vec):
            print(subs("    const double " + sp.ccode(row, assign_to=v)))
    if 'python' in OUTPUT:
        if comment:
            print("\n    ## " + comment)
        if 'math' in OUTPUT:
            for line in desc.split("\n"):
                print("    # " + line)
        if post:
            print("    ## " + post)
        for v, row in zip(var, vec):
            print(subs("    " + sp.ccode(row, assign_to=v)[:-1]))

    if OUTPUT == 'math ':
        print("\n// " + comment if comment else "")
        if post: print("// " + post)
        print(desc)

def mprint(var, mat, comment=None, post=None):
    """
    Generate assignment statements for matrix elements.
    """
    n = sp.prod(var.shape)
    vprint(var.reshape(n, 1), mat.reshape(n, 1), comment=comment, post=post)

# From wikipedia:
#    https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
def Rx(a):
    """Rotate y and z about x"""
    R = [[1, 0, 0],
         [0, +cos(a), -sin(a)],
         [0, +sin(a), +cos(a)]]
    return Matrix(R)

def Ry(a):
    """Rotate x and z about y"""
    R = [[+cos(a), 0, +sin(a)],
         [0, 1, 0],
         [-sin(a), 0, +cos(a)]]
    return Matrix(R)

def Rz(a):
    """Rotate x and y about z"""
    R = [[+cos(a), -sin(a), 0],
         [+sin(a), +cos(a), 0],
         [0, 0, 1]]
    return Matrix(R)


## ===============  Describe the transforms ====================

# Define symbols used.  Note that if you change the symbols for the jitter
# angles, you will need to update the subs() function accordingly.
dphi, dpsi, dtheta = sp.var("phi^+ psi^+ theta^+")
phi, psi, theta = sp.var("phi psi theta")
#dphi, dpsi, dtheta = sp.var("beta^+ gamma^+ alpha^+")
#phi, psi, theta = sp.var("beta gamma alpha")
x, y, z = sp.var("x y z")
q = sp.var("q")
qx, qy, qz = sp.var("qx qy qz")
dqx, dqy, dqz = sp.var("qx^+ qy^+ qz^+")
qa, qb, qc = sp.var("qa qb qc")
dqa, dqb, dqc = sp.var("qa^+ qb^+ qc^+")
qab = sp.var("qab")

# 3x3 matrix M
J = Matrix([sp.var("J(1:4)(1:4)")]).reshape(3,3)
V = Matrix([sp.var("V(1:4)(1:4)")]).reshape(3,3)
R = Matrix([sp.var("R(1:4)(1:4)")]).reshape(3,3)

# various vectors
xyz = Matrix([[x], [y], [z]])
x_hat = Matrix([[x], [0], [0]])
y_hat = Matrix([[0], [y], [0]])
z_hat = Matrix([[0], [0], [z]])
q_xy = Matrix([[qx], [qy], [0]])
q_abc = Matrix([[qa], [qb], [qc]])
q_xyz = Matrix([[qx], [qy], [qz]])
dq_abc = Matrix([[dqa], [dqb], [dqc]])
dq_xyz = Matrix([[dqx], [dqy], [dqz]])

def print_steps(jitter, jitter_inv, view, view_inv, qc_only):
    """
    Show the forward/reverse transform code for view and jitter.
    """
    if 0:  # forward calculations
        vprint(q_xyz, jitter*q_abc, "apply jitter")
        #vprint(xyz, jitter*z_hat, "r")
        #mprint(J, jitter, "forward jitter")
        vprint(dq_xyz, view*q_xyz, "apply view after jitter")
        #mprint(V, view, "forward view")

        #vprint(dq_xyz, view*jitter*q_abc, "combine view and jitter")
        mprint(R, view*jitter, "forward matrix")

    if 1:  # reverse calculations
        pre_view = "set angles from view" if REUSE_SINCOS else None
        pre_jitter = "set angles from jitter" if REUSE_SINCOS else None
        index = slice(2,3) if qc_only else slice(None,None)

        comment("\n**** direct ****")
        vprint(q_abc, view_inv*q_xy, "reverse view", post=pre_view)
        vprint(dq_abc[index,:], (jitter_inv*q_abc)[index,:],
               "reverse jitter after view", post=pre_jitter)

        comment("\n\n**** precalc ****")
        #vprint(q_abc, jitter_inv*view_inv*q_xy, "combine jitter and view reverse")
        mprint(V[:,:2], view_inv[:,:2], "reverse view matrix", post=pre_view)
        mprint(J[index,:], jitter_inv[index,:], "reverse jitter matrix", post=pre_jitter)
        mprint(R[index,:2], (J*V)[index,:2], "reverse matrix")
        comment("\n**** per point ****")
        mprint(q_abc[index,:], (R*q_xy)[index,:], "applied reverse matrix")
        #mprint(q_abc, J*V*q_xy, "applied reverse matrix")
        #mprint(R[index,:2], jitter_inv*view_inv, "reverse matrix direct")

        #vprint(q_abc, M*q_xy, "matrix application")

if 1:
    comment("==== asymmetric ====")
    print_steps(
        jitter=Rx(dphi)*Ry(dtheta)*Rz(dpsi),
        jitter_inv=Rz(-dpsi)*Ry(-dtheta)*Rx(-dphi),
        view=Rz(phi)*Ry(theta)*Rz(psi),
        view_inv=Rz(-psi)*Ry(-theta)*Rz(-phi),
        qc_only=False,
    )

if 1:
    comment("\n\n==== symmetric ====")
    print_steps(
        jitter=Rx(dphi)*Ry(dtheta),
        jitter_inv=Ry(-dtheta)*Rx(-dphi),
        view=Rz(phi)*Ry(theta),
        view_inv=Ry(-theta)*Rz(-phi),
        qc_only=QC_ONLY,
    )

    comment("\n**** qab from qc ****")
    # The indirect calculation of qab is better than directly c
    # alculating qab^2 = qa^2 + qb^2 since qc can be computed
    # as qc = M31*qx + M32*qy, thus requiring only two elements
    # of the rotation matrix.
    #vprint(qab, sqrt(qa**2 + qb**2), "Direct calculation of qab")
    vprint(dqa, sqrt((qx**2+qy**2) - dqc**2),
        "Indirect calculation of qab, from qab^2 = |q|^2 - qc^2")

if 0:
    comment("==== asymmetric (3.x) ====")
    view_inv = Rz(-psi)*Rx(theta)*Ry(-(pi/2 - phi))
    vprint(q_abc, view_inv*q_xy, "reverse view")
    print("""  existing code
    cos_alpha = cos_theta*cos_phi*qxhat + sin_theta*qyhat;
    cos_mu = (-sin_theta*cos_psi*cos_phi - sin_psi*sin_phi)*qxhat + cos_theta*cos_psi*qyhat;
    cos_nu = (-cos_phi*sin_psi*sin_theta + sin_phi*cos_psi)*qxhat + sin_psi*cos_theta*qyhat;
    """)
