// To force a fixed rather than adaptive integration scheme, replace [..., "lib/adaptive.c", ...]
// with [..., "lib/gauss<n>.c", "lib/nonadaptive.c", ...] in your source lists.

#define GAUSS_76 0.0 // Doesn't matter since qr is ignored
#define GAUSS_500 0.0 // Doesn't matter since qr is ignored

static int
gauss_weights(double qr, int n_outer, constant double **w, constant double **z)
{
    *w = GAUSS_W; *z = GAUSS_Z; return GAUSS_N;
}