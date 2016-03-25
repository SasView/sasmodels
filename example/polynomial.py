from numpy import inf
parameters = [
    ["n", "", 0, [0,5], "", "polynomial degree"],
    ["c[n]", "", 0, [-inf, inf], "", "coefficients to c_n x^n"],
]

Iq = """
    int int_n = (int)n;
    double result = c[int_n];
    for (int k=int_n-1; k >= 0; k--) { result = result*q + c[k]; }
    return result;
    """