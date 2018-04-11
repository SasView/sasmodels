from numpy import inf
parameters = [
    ["n", "", 1, [1,5], "", "number of coefficients (or degree+1)"],
    ["c[n]", "", 0, [-inf, inf], "", "coefficients to c_n x^n"],
]

Iq = r"""
    int int_n = (int)n;
    double result = c[int_n-1];
    for (int k=int_n-2; k >= 0; k--) { result = result*q + c[k]; }
    return result;
    """