from numpy import exp, ceil, log2, pi, arange, interp
from numpy.fft import fft

def TInvFourier(x, y, deltaX):
    """
    Inverse Fourier transform implementation

    Parameters:
    x: array of x values
    y: array of y values
    deltaX: step size

    Returns:
    t: array of t values
    Yt: array of transformed values
    """
    deltaX = abs(deltaX)

    Ntemp = (x[-1] - x[0]) / deltaX

    N = 2**ceil(log2(Ntemp))

    X = arange(0, N) * (x[-1] - x[0]) / N + x[0]
    Y = interp(X, x, y)

    R_x = (X[-1] - X[0]) * N / (N - 1)
    deltaX = R_x / N
    x_min = x[0]

    t_min = 0  # -pi/deltaX
    # t_max = 2 * pi / deltaX
    R_t = 2 * pi / deltaX

    i = 1j
    Y = Y * exp(-i * t_min * (arange(1, N+1) - 1) * R_x / N)
    Yt = fft(Y) * deltaX
    Yt = Yt * exp(-i * t_min * x_min)
    Yt = Yt * exp(-i * x_min * R_t * (arange(1, N+1) - 1) / N) / (2 * pi)

    t = t_min + R_t * (arange(1, N+1) - 1) / N

    # Commented out in original:
    # t = np.flip(t)
    # Yw(k) = sum_1^N x(j)* exp(-2*pi*i*(j-1)(k-1)/N)

    return t, Yt