from numpy import exp, pi, arange, linspace, abs, ceil, log2, interp
from numpy.fft import fft

def TFourier(x, y, deltaX):
    """
    Compute the Fourier transform of a function y(x) with sampling interval deltaX.

    Parameters:
    x (array): x coordinates
    y (array): y coordinates
    deltaX (float): sampling interval (must be positive)

    Returns:
    w (array): frequency coordinates
    Yw (array): Fourier transform values
    """
    # deltaX must be positive
    deltaX = abs(deltaX)

    Ntemp = (x[-1] - x[0]) / deltaX

    N = 2**ceil(log2(Ntemp))

    X = linspace(x[0], x[-1], N)
    Y = interp(X, x, y)

    R_x = (X[-1] - X[0]) * N / (N - 1)
    deltaX = R_x / N
    w_min = -1 * pi / deltaX
    # w_max = pi / deltaX
    R_w = 2 * pi / deltaX
    x_min = x[0]

    i = complex(0, 1)
    Y = Y * exp(-i * w_min * arange(N) * R_x / N)
    Yw = fft(Y) * deltaX
    Yw = Yw * exp(-i * w_min * x_min)
    Yw = Yw * exp(-i * x_min * R_w * arange(N) / N)

    w = w_min + R_w * arange(N) / N
    # Yw(k) = sum_1^N x(j)* exp(-2*pi*i*(j-1)(k-1)/N)

    return w, Yw