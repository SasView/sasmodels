// To force a fixed rather than adaptive integration scheme, replace [..., "lib/adaptive.c", ...]
// with [..., "lib/gauss<n>.c", "lib/nonadaptive.c", ...] in your source lists.

// Hack for barbell and capped cylinder keeps the outer integral to 76 points or fewer
// by calling gauss_weights with n_outer=ADAPTIVE_MAX_76
#define ADAPTIVE_MAX_76 1 // Doesn't matter since n_outer is ignored
#define ADAPTIVE_MAX_OUTER 1 // Doesn't matter since n_outer is ignored

static int
gauss_weights(double qr, int n_outer, constant double **w, constant double **z)
{
    *w = GAUSS_W; *z = GAUSS_Z; return GAUSS_N;
}