import numpy as np
from numpy import pi, inf, mean

from .Ecoefficient import TYCoeff
from .CalcRealRoot import CalRealRoot
from .TInvFourier import TInvFourier

# Supplied Q vector must go out to at least this value to calculate g(r).
Q_UPPER = 700
Q_STEP = 0.04 # Suggested step size for Q

# Suggested parameter bounds. These are not enforced in the code, allowing
# the exploration of limiting cases. For now the limits should be implemented
# by the caller.
K_MIN = 1e-4  # Minimum absolute value for K
Z_MIN = 1e-2  # Minimum value for Z
Z_MIN_DIFF = 1e-2 # Minimum difference between Z1 and Z2

def CalTYSk(Z1, Z2, K1, K2, volF, Q, warnFlag=True, debugFlag=False):
    """
    Python implementation of the MATLAB CalTYSk function

    Parameters:
    Z1, Z2 : float
        Parameters for the two-Yukawa potential
    K1, K2 : float
        Parameters for the two-Yukawa potential
    volF : float
        Volume fraction
    Q : array_like
        Wave vector values

    Returns:
    calSk : array_like
        Calculated structure factor
    rootCounter : int
        Number of valid roots found
    calr : array_like
        Spatial distance variable for pair distribution function
    calGr : array_like
        Pair distribution function g(r)
    errorCode : int
        Error code (0: no good roots but some results, 1: one good root, >1: multiple good roots, -1: error)
    cVar : array_like or 0
        Coefficient variables if choice=1, otherwise 0
    """

    # # Check if maximum Q is sufficient
    # if np.max(Q) < Q_UPPER:
    #     print('Maximum Q is too small, possible error when checking g(r)')
    #     print(f'Please increase the maximum Q at least to {Q_UPPER}')
    #     calSk = np.zeros_like(Q)
    #     rootCounter = -1
    #     calr = 0
    #     calGr = 1
    #     errorCode = -1
    #     cVar = 0
    #     return calSk, rootCounter, calr, calGr, errorCode, cVar

    # # Check Q spacing
    # if np.max(Q) / len(Q) > 0.05:
    #     print('Please make the interval of neighbouring Q smaller')
    #     print('max(Q)/length(Q) > 0.05')
    #     calSk = np.zeros_like(Q)
    #     rootCounter = -1
    #     calr = 0
    #     calGr = 1
    #     errorCode = -1
    #     cVar = 0
    #     return calSk, rootCounter, calr, calGr, errorCode, cVar

    # Note: we are using the opposite sign convention for K1 and K2 in python compared
    # to the matlab code. We are reversing it here rather than updating the rest of the
    # code to minimize the difference between this code and matlab.
    coeff = TYCoeff(Z=[Z1, Z2], K=[-K1, -K2], phi=volF)

    # Calculate the roots for d1 and d2
    Rd1, Rd2 = CalRealRoot(coeff)
    Rd1 = Rd1 * coeff.d1Factor
    Rd2 = Rd2 * coeff.d2Factor

    # Initialize errorCode and rootCounter
    errorCode = 0
    rootCounter = 0
    goodRoot = 0

    # Arrays to store results
    calCoeArray = []
    calrArray = []
    calGrArray = []
    calSkArray = []
    goodRootPos = []

    # For debugging purposes keep track of the roots we are considering and their status.
    scanned_roots = []  # tuple[float,str,int] giving root:flag:position

    # kk is used to compute g(r) for selecting between available solutions
    # Use a fixed Q range for this calculation
    kk = np.arange(0, Q_UPPER + Q_STEP, Q_STEP)
    # kk = Q
    for i in range(len(Rd2)):
        for j in range(2):
            # Check if there is NaN root
            if np.isnan(Rd1[i, j]) or Rd2[i] == 0:
                continue

            # Based on the solution of d1 and d2, calculate other coefficients in Blum's method
            coe = coeff.CxCoef(Rd1[i, j], Rd2[i])

            # Assign the calculated result to different variables
            # coe contains [a00, b00, v1, v2, a, b, c1, d1, c2, d2]
            a00, b00, v1, v2, *rest = coe

            # Calculate C(k)
            eCk = coeff.Ck(a00, b00, v1, v2, kk)

            # Calculate h(k) and S(k)
            ehk = coeff.hk(eCk)
            eSk = coeff.Sk(ehk)

            # Calculate h(r) and g(r) from h(k)
            rh, hc_r = TInvFourier(kk, 4 * pi * ehk * kk, 0.2)
            rh, hc_r = rh[1:], hc_r[1:]
            h_r = -np.imag(hc_r) / rh / 4 / pi / pi
            g_r = h_r + 1

            if debugFlag:
                import matplotlib.pyplot as plt

                # Prepare a new figure window without clobbering the current one
                old_fig = plt.gcf().number
                plt.figure()

                # Calculate C(r)
                r = np.arange(0.001, 10.001, 0.01)
                eCr = coeff.Cr(a00, b00, v1, v2, r)

                # Plot C(r)
                plt.subplot(3, 1, 1)
                plt.plot(r, -eCr, 'g*-')
                index = np.where(r > 1)[0]

                tempMin = np.min(-eCr[index])
                tempMax = np.max(-eCr[index])
                maxPos = r[np.argmax(-eCr)]
                plt.axis([r[index[0]], 4, min(-eCr[index]), max(1, tempMax)])
                plt.grid(True)
                plt.text(0.5 * 4, 0.4 * (tempMax - tempMin) + tempMin,
                         f'Maximum value = {tempMax}, Position={maxPos}')
                plt.xlabel('r/\\sigma')
                plt.ylabel('-c(r)')
                plt.title(f'φ={volF}, z(1)={Z1}, z(2)={Z2}, K(1)={K1}, K(2)={K2}')

                # Plot S(k)
                plt.subplot(3,1,2)
                plt.plot(kk, eSk, 'bo-')
                plt.xlabel('k\\sigma')
                plt.ylabel('S(k)')
                plt.title(f'φ={volF}, z(1)={Z1}, z(2)={Z2}, K(1)={K1}, K(2)={K2}')
                plt.axis([0, 5*pi, 0, max(eSk)])
                plt.grid()

                # Plot g(r)
                plt.subplot(3,1,3)
                plt.plot(rh, g_r, 'bo-')
                plt.axis([0, 8, None, None])
                plt.xlabel('r/σ')
                plt.ylabel('g(r)')
                plt.title(f'Root: ({i}, {j}) = {Rd1[i, j]:.2f}')
                plt.grid()

                # Restore the original figure
                plt.figure(old_fig)

            hardcore_gr = np.mean(abs(g_r[rh < 0.95][1:]))
            if testGr(rh, g_r, False) == 1:
                scanned_roots.append((Rd1[i, j], "X", -1, hardcore_gr))
                continue

            rootCounter += 1

            if v1*coeff.k[0] >= 0 and v2*coeff.k[1] >= 0:
                goodRoot += 1
                goodRootPos.append(rootCounter-1) # python index arrays are 0-origin
                scanned_roots.append((Rd1[i, j], " ", rootCounter-1, hardcore_gr))
            else:
                scanned_roots.append((Rd1[i, j], "?", rootCounter-1, hardcore_gr))

            calCoeArray.append(coe)
            calrArray.append(rh)
            calGrArray.append(g_r)
            calSkArray.append(eSk)

    # No good roots and reasonable solution
    if rootCounter == 0:
        errorCode = -1
        position = -1
        calr = 0
        calGr = 0
        calSk = np.zeros_like(Q)
        cVar = np.zeros(10)

    # There is some result which may be reasonable
    # There is more than one good solution, only one of them is sent back
    elif goodRoot > 1:
        errorCode = goodRoot
        position = findSGr(calrArray, calGrArray, goodRootPos)
        # position = goodRootPos[0]
        calr = calrArray[position]
        calGr = calGrArray[position]
        calSk = calSkArray[position]
        cVar = calCoeArray[position]
        testGr(calr, calGr, warnFlag)

    # There is only one good result
    elif goodRoot == 1:
        errorCode = 1
        position = goodRootPos[0]
        calr = calrArray[position]
        calGr = calGrArray[position]
        calSk = calSkArray[position]
        cVar = calCoeArray[position]
        testGr(calr, calGr, warnFlag)

    # Sorry, no good result, but some of them can be sent back to be checked
    elif goodRoot == 0:
        errorCode = 0
        # position = 0
        position = findSGr(calrArray, calGrArray)
        calr = calrArray[position]
        calGr = calGrArray[position]
        calSk = calSkArray[position]
        cVar = calCoeArray[position]
        testGr(calr, calGr, warnFlag)

    else:
        raise AssertionError('Internal error in CalTYSk root checking logic.')

    if warnFlag:
        print("roots:", " ".join(f"{root: 7.2f}{flag}({g_r:.3f}){'*' if index==position else ' '}" for root, flag, index, g_r in sorted(scanned_roots)))

    # Only show plots if the last root is not the selected root
    if scanned_roots and debugFlag:
        min_index = np.argmin([abs(r) for r, flag, index, g_r in scanned_roots])
        if scanned_roots[min_index][2] != position:
            import matplotlib.pyplot as plt
            plt.show()

    # TODO: Use stable numerics for Ck and hk for small Q and remove the max(Q, 0.01) hack
    # Compute S(Q) at the desired Q
    if errorCode >= 0:
        a00, b00, v1, v2, *rest = cVar
        eCk = coeff.Ck(a00, b00, v1, v2, np.maximum(Q, 0.01))
        ehk = coeff.hk(eCk)
        calSk = coeff.Sk(ehk)
    else:
        calSk = Q*0.0


    return calSk, rootCounter, calr, calGr, errorCode, cVar


def findSGr(calrArray, calGrArray, positionArray=None):
    """
    Find the position with the smallest sum of absolute values in the hardcore region

    Parameters:
    calrArray : list of arrays
        List of r arrays
    calGrArray : list of arrays
        List of g(r) arrays
    positionArray : list, optional
        List of positions to consider

    Returns:
    position : int
        Position with the smallest sum in the hardcore region
    """
    if positionArray is None:
        sumHardcore = 1000
        position = 0
        for i in range(len(calrArray)):
            index = np.where(calrArray[i] < 0.95)[0]
            if len(index) > 1:  # Ensure there are elements after index 0
                tempSum = np.sum(np.abs(calGrArray[i][1:len(index)])) / (len(index) - 1)
                if sumHardcore > tempSum:
                    sumHardcore = tempSum
                    position = i

        return position
    else:
        position = 0
        sumHardcore = 1000

        for i in range(len(positionArray)):
            j = positionArray[i]
            index = np.where(calrArray[j] < 0.95)[0]
            if len(index) > 1:  # Ensure there are elements after index 0
                tempSum = np.sum(np.abs(calGrArray[j][1:len(index)])) / (len(index) - 1)
                if sumHardcore > tempSum:
                    sumHardcore = tempSum
                    position = j

        return position

def testGr(r, g_r, warnFlag=True):
    """
    Test the g(r) function for various conditions and return a flag indicating the result.
    """
    # global debugFlag  # Uncomment if needed

    # TODO: Test against max absolution value of g(r)
    if max(g_r) > 20:
        if warnFlag:
            print('Warning! Maximum value of g(r) > 20')
        # return 1

    # TODO: Average over abs(g_r) rather than g_r.
    average = mean(g_r[r < 1])

    # TODO: Uncomment to rejection really bad solutions.
    # if average < -10:
    #     if warnFlag:
    #         print('Negative warning g(r), average < -10 inside hardcore')
    #     return 1
    if average < -1:
        if warnFlag:
            print('Negative warning g(r), average < -1 inside hardcore')
        return -2
    if average < -0.3:
        if warnFlag:
            print('Negative warning g(r), average < -0.3 inside hardcore')
        return -1

    # TODO: The r=0 element was already removed, right after the fourier transform.
    # Recompute index but skip the first element
    average = mean(g_r[r < 1][1:])

    # if average > 10:
    #     if warnFlag:
    #         print('Positive warning g(r), averat > 10 inside hardcore')
    #     return 1
    if average > 1:
        if warnFlag:
            print('Positive warning g(r), average > 1 inside hardcore')
        return -2
    if average > 0.3:
        if warnFlag:
            print('Positive warning g(r), average > 0.3 inside hardcore')
        return -1

    return 0

def demo():
    import matplotlib.pyplot as plt

    # Potential From: V(r)/(kB*T)=-K1*exp(-Z1*(r-1))/r - K2*exp(-Z2*(r-1))/r
    # kB: Boltzman constant
    # T: Absolute temperature
    # The hardcore diameter is assumed to be one
    # K1, K2 follow the sign convention of the matlab code
    # In Matlab, negative K1 means attraction, positive means repulsion
    # In Python, negative K1 means repulsion, positive means attraction

    Z1=2
    Z2=10
    K1=-1         # Attraction
    K2=6          # Repulsion

    # Volume Fraction
    volF=0.2

    # Use matlab syntax for defining parameters to ease comparison. Note that
    # K1, K2 are negative in matlab and positive in python.
    K1=-6; K2=2.7; volF=0.1; Z1=10; Z2=2.3;

    # Please give your Q range, the maximum Q should be larger than 700 to let
    # the program to run.
    # (The reason is : we have to use g(r) to select the right result. Therefore
    # large Q range can make g(r) better.
    Q = 2*pi*np.arange(0, 15*10, 0.005)

    [calSk,rootCounter,calr,calGr,errorCode, cVar] = CalTYSk(Z1,Z2,-K1,-K2,volF,Q)

    # Plot structure factor S(Q)
    plt.subplot(211)
    plt.plot(Q, calSk, 'b-')
    plt.xlim([0, 20])
    #plt.xscale('log')xs
    #plt.axis([0,20,0,inf])

    # Plot the pair distribution function g(r)
    plt.subplot(212)
    plt.plot(calr, calGr, 'b-')
    plt.xlim([0, 10])
    #plt.xscale('log')
    #plt.axis([0,10,-inf,inf])
    plt.show()

if __name__ == "__main__":
    demo()
