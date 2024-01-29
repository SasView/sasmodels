from collections import namedtuple
import numpy as np
from numpy import sqrt, exp, expm1

AVOGADRO = 6.022e23

Polymer = namedtuple("Polymer", "n phi v a b".split())
def E(q, poly):
    qrsq = (q*Rg(poly))**2
    retval = exp(-qrsq)
    return retval

def F(q, poly):
    qrsq = (q*Rg(poly))**2
    retval = -expm1(-qrsq)/qrsq
    return retval

def P_ii(q, poly):
    qrsq = (q*Rg(poly))**2
    retval = 2 * (expm1(-qrsq) + qrsq)/qrsq**2
    return retval

def P_ij(q, poly_list):
    i, j = poly_list[0], poly_list[-1]
    retval = F(q, i) * np.prod([E(q,p) for p in poly_list[1:-1]]) * F(q, j)
    return retval

def Rg(poly):
    return sqrt(poly.n/6.)*poly.a

def S0_ii(q, poly):
    retval = poly.n*poly.phi*poly.v*P_ii(q, poly)
    return retval

def S0_ij(q, poly_list):
    i,j = poly_list[0], poly_list[-1]
    retval = sqrt(i.n*i.phi*i.v*j.n*j.phi*j.v) * P_ij(q, poly_list)
    return retval

def drho(poly, base):
    return (poly.b/poly.v -  base.b/base.v)*1e-13*sqrt(AVOGADRO)

def binary_blend(q, C, D, Kcd):
    """
    de Gennes, 1979
    """
    S0cc = S0_ii(q, C)
    S0dd = S0_ii(q, D)
    vcc = 1/S0dd - 2*Kcd #/v0
    Scc = S0cc/(1 + vcc*S0cc)
    rhocd = drho(C,D)
    Iq = rhocd**2 * Scc
    return Iq

def ternary_blend(q, B, C, D, Kbc, Kbd, Kcd):
    S0bb = S0_ii(q, B)
    S0cc = S0_ii(q, C)
    S0dd = S0_ii(q, D)
    vbb = 1/S0dd - 2*Kbd
    vcc = 1/S0dd - 2*Kcd
    vbc = 1/S0dd + Kbd - Kbc - Kcd
    rhobd = drho(B,D)
    rhocd = drho(C,D)
    det = (1+vbb*S0bb)*(1+vcc*S0cc) - vbc**2*S0bb*S0cc
    Sbb = S0bb*(1+vcc*S0cc)/det
    Scc = S0cc*(1+vbb*S0bb)/det
    Sbc = -S0bb*vbc*S0cc/det
    Iq = rhobd**2*Sbb + rhocd**2*Scc + 2*rhobd*rhocd*Sbc
    return Iq

def diblock_copolymer(q, C, D, Kcd):
    """
    Leibler, 1980
    """
    S0cc = S0_ii(q, C)
    S0dd = S0_ii(q, D)
    S0cd = S0_ij(q, [C, D])
    Scc = (S0cc*S0dd - S0cd**2)/((S0cc+S0dd + 2*S0cd)-2*Kcd*(S0cc+S0dd-2*S0cd))
    rhocd = drho(C,D)
    Iq = rhocd**2 * Scc
    return Iq

def matrix_form(q, case_num, polys,
                Kab=0., Kac=0., Kad=0., Kbc=0., Kbd=0., Kcd=0.):
    if case_num < 2:
        C, D = polys
        A = B = D
    elif case_num < 5:
        B, C, D = polys
        A = D
    else:
        A, B, C, D = polys

    rho = np.array([[drho(p, D)] for p in (A,B,C)])
    S0aa = S0_ii(q, A)
    S0bb = S0_ii(q, B)
    S0cc = S0_ii(q, C)
    S0ab = S0_ij(q, [A, B])
    S0ac = S0_ij(q, [A, B, C])
    S0bc = S0_ij(q, [B, C])
    if case_num == 4: # No a-c interaction in triblock copolymer
        S0ac *= 0.0
    elif case_num == 9: # No a-c or a-d interaction in tetrablock copolymer
        S0ac *= 0.0

    S0 = np.array([[S0aa, S0ab, S0ac], [S0ab, S0bb, S0bc], [S0ac, S0bc, S0cc]])
    S0dd = S0_ii(q, D)
    vaa = 1./S0dd - 2*Kad
    vbb = 1./S0dd - 2*Kbd
    vcc = 1./S0dd - 2*Kcd
    vab = 1./S0dd + Kab - Kad - Kbd
    vac = 1./S0dd + Kac - Kad - Kcd
    vbc = 1./S0dd + Kbc - Kbd - Kcd
    v = np.array([[vaa, vab, vac], [vab, vbb, vbc], [vac, vbc, vcc]])

    Iq = np.empty_like(q)
    for k, qk in enumerate(q):
        S0_k = S0[:,:,k].M
        v_k = v[:,:,k].M
        S_k = np.linalg.inv(1 + S0_k*v_k)*S0_k
        Iq[k] = rho.T @ S_k @ rho

def build_pars(case_num, polys, **interactions):
    def set_poly(x, poly):
        pars["N"+x] = poly.n
        pars["Phi"+x] = poly.phi
        pars["v"+x] = poly.v
        pars["b"+x] = poly.a
        pars["L"+x] = poly.b
    pars = interactions.copy()
    pars["case_num"] = case_num
    polys = list(reversed(polys))
    if len(polys) >= 4: set_poly("a",polys[3])
    if len(polys) >= 3: set_poly("b",polys[2])
    if len(polys) >= 2: set_poly("c",polys[1])
    if len(polys) >= 1: set_poly("d",polys[0])
    return pars

def sasmodels_rpa(q, pars):
    from sasmodels.models import rpa
    from sasmodels.core import load_model
    from sasmodels.direct_model import DirectModel
    from sasmodels.data import empty_data1D
    data = empty_data1D(q, resolution=0.0)
    model = load_model(rpa, dtype="double", platform="dll")
    #model = load_model(rpa, dtype="single", platform="ocl")
    M = DirectModel(data, model)
    return M(**pars)

def sasview_rpa(q, pars):
    from sasmodels.models import rpa
    from sasmodels.compare import eval_sasview
    from sasmodels.data import empty_data1D
    data = empty_data1D(q, resolution=0.0)
    M = eval_sasview(rpa, data)
    return M(**pars)

def demo():
    import sys
    case_num = 0 if len(sys.argv) < 2 else int(sys.argv[1])

    B = Polymer(n=525,phi=0.57,v=97.5,b=-4.99,a=8)
    C = Polymer(n=525,phi=0.57,v=97.5,b=-4.99,a=8)
    D = Polymer(n=1105,phi=0.43,v=81.9,b=53.1,a=2)
    q = np.logspace(-4,-1,400)
    #q = np.array([0.1])
    Kab=Kac=Kad=0.0
    Kcd = 0.0106*0.0035 - 1.84e-5
    Kbd = Kcd + 2e-4
    Kbc = (Kcd + 1e-4)*0.5
    K = dict(Kab=Kab,Kac=Kac,Kad=Kad,Kbc=Kbc,Kbd=Kbd,Kcd=Kcd)
    if case_num == 0:
        Iq = binary_blend(q, C, D, Kcd)
    elif case_num == 1:
        Iq = diblock_copolymer(q, C, D, Kcd)
    elif case_num == 2:
        Iq = ternary_blend(q, B, C, D, Kbc, Kbd, Kcd)
    else:
        raise ValueError("Case %d not implmented"%case_num)

    pars = build_pars(case_num, [B, C, D], **K)
    print "eval sasmodels"
    Iq_sasmodels = sasmodels_rpa(q, pars)
    print "eval sasview"
    Iq_sasview = sasview_rpa(q, pars)
    print 1./Iq[0], 1./Iq_sasmodels[0], 1./Iq_sasview[0]

    #return
    import pylab
    pylab.subplot(121)
    pylab.loglog(q, Iq, label='direct')
    pylab.loglog(q, Iq_sasmodels, label='sasmodels')
    pylab.loglog(q, Iq_sasview, label='sasview')
    pylab.legend()
    pylab.subplot(122)
    pylab.loglog(q, abs(Iq_sasview - Iq_sasmodels)/Iq_sasmodels, label='sasview-sasmodels')
    pylab.loglog(q, abs(Iq_sasmodels - Iq)/Iq, label='sasmodels-direct')
    pylab.legend()
    pylab.show()

if __name__ == "__main__":
    demo()
