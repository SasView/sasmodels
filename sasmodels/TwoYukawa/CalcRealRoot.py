import numpy as np

from .Epoly import make_poly

def CalRealRoot(coeff):

    Ecoefficient = coeff.Ecoefficient

    # Ick! Code is duplicated in Epoly.py
    gE1d20 = Ecoefficient[0]
    gE1d11 = Ecoefficient[1]
    gE1d21 = Ecoefficient[2]
    gE1d31 = Ecoefficient[3]
    gE1d02 = Ecoefficient[4]
    gE1d12 = Ecoefficient[5]
    gE1d22 = Ecoefficient[22]
    gE1d32 = Ecoefficient[23]
    gE1d42 = Ecoefficient[6]
    gE1d13 = Ecoefficient[7]
    gE1d23 = Ecoefficient[8]
    gE1d33 = Ecoefficient[9]
    gE1d24 = Ecoefficient[10]

    gE2d20 = Ecoefficient[11]
    gE2d11 = Ecoefficient[12]
    gE2d21 = Ecoefficient[13]
    gE2d31 = Ecoefficient[14]
    gE2d02 = Ecoefficient[15]
    gE2d12 = Ecoefficient[16]
    gE2d22 = Ecoefficient[24]
    gE2d32 = Ecoefficient[25]
    gE2d42 = Ecoefficient[17]
    gE2d13 = Ecoefficient[18]
    gE2d23 = Ecoefficient[19]
    gE2d33 = Ecoefficient[20]
    gE2d24 = Ecoefficient[21]

    d2Coeff = make_poly(Ecoefficient)
    # print("Polynomial coefficients:"); print(d2Coeff)
    # print("Poly: "+ " ".join(f"{coeff:+8.1e}" for coeff in d2Coeff))

    N=len(d2Coeff)

    if 1:
        Factor = 1  # 20
        FactorArray = Factor ** np.arange(N)[::-1]
        d2Coeff_factor = d2Coeff*FactorArray
        myd2Root_factor = np.roots(d2Coeff_factor)
        myd2Root = myd2Root_factor*Factor
    else:
        # Very slow high precision root finder
        import mpmath
        myd2Root = mpmath.polyroots(d2Coeff, maxsteps=500, extraprec=1000)

    if 0:
        import mpmath
        myd2Root_mpmath = mpmath.polyroots(d2Coeff, maxsteps=500, extraprec=1000)
        nproots = list(sorted(myd2Root, key=lambda v: (v.real, v.imag)))
        # print(nproots/Factor)
        mproots = list(complex(v) for v in sorted(myd2Root_mpmath, key=lambda v: (v.real,v.imag)))
        for k, (npr, mpr) in enumerate(zip(nproots, mproots)):
            if (mpr.imag == 0 and mpr != 0):
                print(f"{k:02} mp:{float(mpr.real):.2f} np:{npr.real:.2f} {abs((npr-mpr)/(abs(mpr) + (mpr==0))):8.2e}")

    myRd2Root = []
    myRd1Root = []

    # Solve Eq 15 in DOI:10.1063/1.1830433, giving two d1 solutions for each proposed d2 root.
    # Check both d1 against Eq 14, rejecting those that that are not near zero.
    for d2 in myd2Root:
        # Skip imaginary d2 roots
        if np.imag(d2) != 0:
            if abs(d2.imag) < 1e-12*abs(d2.real):
                print("Skipping nearly real {d2}")
            continue

        d2 = float(np.real(d2))

        # Coeff of Equation 14
        ycoe14 = gE1d42 * d2
        ycoe13 = gE1d31 + gE1d32 * d2 + gE1d33 * d2**2
        ycoe12 = gE1d22 * d2
        ycoe11 = gE1d11 + gE1d12 * d2 + gE1d13 * d2**2
        ycoe10 = gE1d02 * d2
        def checkd1(d1):
            """Verify equation 14 is near zero."""
            # TODO: why are we dividing by d1^2?
            return d1 if np.abs(ycoe14 * d1**2 + ycoe13 * d1 + ycoe12 + ycoe11 / d1 + ycoe10 / d1**2) < 1e-5 else np.nan

        # Coeff of Equation 15
        ycoe22 = gE2d31 * d2 + gE2d33 * d2**3
        ycoe21 = gE2d20 + gE2d21 * d2 + gE2d22 * d2**2 + gE2d23 * d2**3 + gE2d24 * d2**4
        ycoe20 = gE2d11 * d2 + gE2d13 * d2**3

        discriminant = ycoe21**2 - 4 * ycoe22 * ycoe20
        if discriminant < 0:
            continue

        if ycoe22 == 0:
            # Eq. 15 is linear. Avoid divide by zero when looking for a pair of quadratic roots.
            # TODO: should we skip this root entirely?
            # TODO: shouldn't the solution to $y2_1 d_1 + y2_0 = 0$ be $d_1 = -y2_0 / y2_1$?
            d1p = 0
            d1m = np.nan
        else:
            d1p = checkd1((-ycoe21 + np.sqrt(discriminant)) / (2 * ycoe22))
            d1m = checkd1((-ycoe21 - np.sqrt(discriminant)) / (2 * ycoe22))
        d1pair = [d1p, d1m] if not np.isnan(d1p) else [d1m, d1p]

        myRd2Root.append(d2)
        myRd1Root.append(d1pair)

    # Check the roots
    if 0:
        print(myRd2Root)
        print(myRd1Root)
        for root in myRd2Root:
            val = np.polyval(d2Coeff, root)
            print(f"{'X ' if abs(val)>1e-10 else '  '} P({root}) = {val}")
    if 0:
        testRoot = np.array(myRd2Root)
        resi = testRoot * 0
        print(f"{testRoot[0]=}")
        for k, coeff in enumerate(d2Coeff[::-1]):
            resi = coeff + resi * testRoot
            # print(f"{len(d2Coeff) - k:2} {coeff=} {resi[0]=}")
        print(resi)

    return np.array(myRd1Root), np.array(myRd2Root)
