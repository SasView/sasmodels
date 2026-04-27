//#include <math.h>
//#include <stdio.h>

//truncated octahedron volume 
// NOTE: needs to be called form_volume() for a shape category
static double
form_volume(double length_a, double b2a_ratio, double c2a_ratio, double truncation)
{
// length_a is the half height along the a axis of the octahedron without truncature
// length_b is the half height along the b axis of the octahedron without truncature
// length_c is the half height along the c axis of the octahedron without truncature
// b2a_ratio is length_b divided by Length_a
// c2a_ratio is Length_c divided by Length_a
// truncation varies from 0 (octahedron) to 0.5 (cuboctahedron)
// updated convention with truncation=0 for on truncation
// truncation and tinv are exchanged starting from the previous version of the code
    const double tinv = 1.0 - truncation;
    return (4./3.) * cube(length_a) * b2a_ratio * c2a_ratio *(1.-3*cube(truncation));
}


// Fq() is called because option "have_Fq = True" is set to True in the Python file
static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double truncation)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;
    const double tinv = 1.0 - truncation;
    // Because we are integrating over the entire surface, the orientation is arbitrary. Sort
    // the lengths Lc > Lb > La so that the integration will go faster.
    const double maybe_min = fmin(length_a, length_b);
    const double maybe_max = fmax(length_a, length_b);
    const double maybe_mid = fmax(maybe_min, length_c);
    const double La = fmin(maybe_min, length_c);
    const double Lb = fmin(maybe_max, maybe_mid);
    const double Lc = fmax(maybe_max, maybe_mid);

    // Find the circumradius by truncating the line from Lc to Lb. This chops of
    // the similar triangle with sides (truncation)Lc, (truncation)Lb, leaving coordinate ((1-truncation)Lc, truncationLb).
    // The distance to the origin then follows.
    const double qr_max_outer = q*sqrt(square(tinv*Lc) + square(truncation*Lb));
    constant double *z_outer, *w_outer;
    int n_outer = gauss_weights(qr_max_outer, &w_outer, &z_outer);

    const double qr_max_inner = q*sqrt(square(tinv*Lb) + square(truncation*La));
    constant double *z_inner, *w_inner;
    int n_inner = gauss_weights(qr_max_inner, &w_inner, &z_inner);

    //printf("La=%g Lb=%g Lc=%g npoints = %d x %d = %d\n", La, Lb, Lc, n_outer, n_inner, n_outer*n_inner);

    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits

    double outer_sum_F1 = 0.0;
    double outer_sum_F2 = 0.0;

    for(int i=0; i<n_outer; i++) {
        // Do a u=cos(theta) substitution for the outer loop
        const double cos_theta = 0.5*(z_outer[i] + 1.0); // [-1, 1] => [0, 1]
        const double sin_theta = sqrt(1.0 - square(cos_theta)); // = sin(acos(cos_theta))
        const double qc = q * cos_theta;
        const double qz = qc * Lc;
        const double qz2 = square(qz);

        double inner_sum_F1 = 0.0;
        double inner_sum_F2 = 0.0;
        for(int j=0; j<n_inner; j++) {
            // 8-fold symmetry, so we only need one quadrant.
            const double phi = M_PI_4*(z_inner[j] + 1.0); // [-1, 1] => [0, pi/2]
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            //HERE: Octahedron formula
            // q is the modulus of the scattering vector in [A-1]
            // NOTE: capital QX QY QZ are the three components in [A-1] of the scattering vector
            // NOTE: qx qy qz are rescaled components (no unit) for computing AA, BB and CC terms
            const double qa = q * sin_theta * cos_phi;
            const double qb = q * sin_theta * sin_phi;
            const double qx = qa * La;
            const double qy = qb * Lb;

            // TODO: test for q=0 and return the limiting value.
            // From the equations, lim q -> 0 seems to be O(1/q^2), so it diverges. And indeed,
            // for q < 1e-8 the function begins to rise.
            // TODO: calculation is unstable for small q.
            // PAK: reordered the equations and moved factor of 1/2 to normalization.
            // PAK: verified that t=1 matches Wei-Ren Chen et al. for a regular octahedron.
            const double qx2 = square(qx);
            const double qy2 = square(qy);
            const double AA =
                ((qy-qx)*sin(qy*truncation-qx*tinv) + (qy+qx)*sin(qy*truncation+qx*tinv)) / ((qy2-qz2)*(qy2-qx2)) +
                ((qz-qx)*sin(qz*truncation-qx*tinv) + (qz+qx)*sin(qz*truncation+qx*tinv)) / ((qz2-qx2)*(qz2-qy2));

            const double BB =
                ((qz-qy)*sin(qz*truncation-qy*tinv) + (qz+qy)*sin(qz*truncation+qy*tinv)) / ((qz2-qx2)*(qz2-qy2)) +
                ((qx-qy)*sin(qx*truncation-qy*tinv) + (qx+qy)*sin(qx*truncation+qy*tinv)) / ((qx2-qy2)*(qx2-qz2));

            const double CC =
                ((qx-qz)*sin(qx*truncation-qz*tinv) + (qx+qz)*sin(qx*truncation+qz*tinv)) / ((qx2-qy2)*(qx2-qz2)) +
                ((qy-qz)*sin(qy*truncation-qz*tinv) + (qy+qz)*sin(qy*truncation+qz*tinv)) / ((qy2-qz2)*(qy2-qx2));

            // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
            const double AP = 3./(1. - 3.*cube(truncation)) * (AA+BB+CC);

            inner_sum_F1 += w_inner[j] * AP;
            inner_sum_F2 += w_inner[j] * AP * AP;

        }
        // Note: moving the pi/2 scaling outside the sum. No sin(theta) term because of u=cos(theta)
        // TODO: With the u=cos(theta) substitution should the F1 integral be negative?
        // TODO: check for this in all models.
        outer_sum_F1 += w_outer[i] * inner_sum_F1;
        outer_sum_F2 += w_outer[i] * inner_sum_F2;
    }
    // sum(w) = 2, so scaling
    outer_sum_F1 *= 0.5*M_PI_4;
    outer_sum_F2 *= 0.5*M_PI_4;

    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    outer_sum_F1 /= M_PI_2;
    outer_sum_F2 /= M_PI_2;

    // Multiply by contrast and volume
    const double s = (sld-solvent_sld) * form_volume(length_a, b2a_ratio,c2a_ratio, truncation);

    // Convert from [1e-12 A-1] to [cm-1]
    *F1 = 1e-2 * s * outer_sum_F1;
    *F2 = 1e-4 * square(s) * outer_sum_F2;

    // TODO: inf and NaN shouldn't occur, but if they do they shouldn't be replaced by zero.
    if (isnan(*F1) || isinf(*F1)) {
        *F1 = 0.0;
    }
    if (isnan(*F2) || isinf(*F2)) {
        *F2 = 0.0;
    }
}


static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double truncation)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;
    const double tinv = 1.0 - truncation;

    //HERE: Octahedron formula
    // NOTE: qa qb qc are the three components in 1/Ang of the scattering vector
    // NOTE: qx qy qz are rescaled components (no unit) for computing AA, BB and CC terms
    const double qx = qa * length_a;
    const double qy = qb * length_b;
    const double qz = qc * length_c;

    // TODO: calculation is unstable for small q.
    // PAK: reordered the equations and moved factor of 1/2 to normalization.
    const double qx2 = square(qx);
    const double qy2 = square(qy);
    const double qz2 = square(qz);
    const double AA =
        ((qy-qx)*sin(qy*truncation-qx*tinv) + (qy+qx)*sin(qy*truncation+qx*tinv)) / ((qy2-qz2)*(qy2-qx2)) +
        ((qz-qx)*sin(qz*truncation-qx*tinv) + (qz+qx)*sin(qz*truncation+qx*tinv)) / ((qz2-qx2)*(qz2-qy2));

    const double BB =
        ((qz-qy)*sin(qz*truncation-qy*tinv) + (qz+qy)*sin(qz*truncation+qy*tinv)) / ((qz2-qx2)*(qz2-qy2)) +
        ((qx-qy)*sin(qx*truncation-qy*tinv) + (qx+qy)*sin(qx*truncation+qy*tinv)) / ((qx2-qy2)*(qx2-qz2));

    const double CC =
        ((qx-qz)*sin(qx*truncation-qz*tinv) + (qx+qz)*sin(qx*truncation+qz*tinv)) / ((qx2-qy2)*(qx2-qz2)) +
        ((qy-qz)*sin(qy*truncation-qz*tinv) + (qy+qz)*sin(qy*truncation+qz*tinv)) / ((qy2-qz2)*(qy2-qx2));

    // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
    const double AP = 3./(1. - 3.*cube(truncation)) * (AA+BB+CC);


    // Multiply by contrast and volume
    // contrast
    const double s = (sld-solvent_sld);
    // volume
    // s *= form_volume(length_a, b2a_ratio,c2a_ratio, truncation);

    // Convert from [1e-12 A-1] to [cm-1]
    double answer = 1.0e-4 * square(s * form_volume(length_a, b2a_ratio,c2a_ratio, truncation) * AP);
    if (isnan(answer) || isinf(answer)) {
        return 0.0;
    }

    return answer;
}



