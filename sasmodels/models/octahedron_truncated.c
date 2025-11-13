#include <math.h>
#include <stdio.h>
static double
form_volume(double length_a, double b2a_ratio, double c2a_ratio, double t)
{
    /* octahedron volume formula */
    /* length_a is the half height along the a axis of the octahedron */
    return (4./3.) * length_a * (length_a * b2a_ratio) * (length_a * c2a_ratio) * (1. - 3 * (1. - t) * (1. - t) * (1. - t));
}

static double
Iq(double q,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;


   //Integration limits to use in Gaussian quadrature
    const double v1a = 0.0;
    const double v1b = M_PI_2;  //theta integration limits
    const double v2a = 0.0;
    const double v2b = M_PI_2;  //phi integration limits

    double outer_sum = 0.0;
    for (int i = 0; i < GAUSS_N; i++) {
        const double theta = 0.5 * (GAUSS_Z[i] * (v1b - v1a) + v1a + v1b);
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);

        double inner_sum = 0.0;
        for (int j = 0; j < GAUSS_N; j++) {
            double phi = 0.5 * (GAUSS_Z[j] * (v2b - v2a) + v2a + v2b);
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            //HERE: Octahedron formula
            const double Qx = q * sin_theta * cos_phi;
    	    const double Qy = q * sin_theta * sin_phi;
    	    const double Qz = q * cos_theta;
    	    const double qx = Qx * length_a;
    	    const double qy = Qy * length_b;
    	    const double qz = Qz * length_c;

            const double AA = 1./((qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qx)*sin(qy*(1.-t)-qx*t)+(qy+qx)*sin(qy*(1.-t)+qx*t))+
                                1./((qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qx)*sin(qz*(1.-t)-qx*t)+(qz+qx)*sin(qz*(1.-t)+qx*t));

            const double BB = 1./((qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qy)*sin(qz*(1.-t)-qy*t)+(qz+qy)*sin(qz*(1.-t)+qy*t))+
                                1./((qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qy)*sin(qx*(1.-t)-qy*t)+(qx+qy)*sin(qx*(1.-t)+qy*t));

            const double CC = 1./((qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t))+
                                1./((qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));


	    // normalisation to 1. of AP at q = 0. Division by a Factor 4/3.
            const double AP = 6./(1.-3*(1.-t)*(1.-t)*(1.-t))*(AA+BB+CC);

            inner_sum += GAUSS_W[j] * AP * AP;


        }
        inner_sum = 0.5 * (v2b - v2a) * inner_sum;
        outer_sum += GAUSS_W[i] * inner_sum * sin_theta;
    }

    double answer = 0.5 * (v1b - v1a) * outer_sum;

    // Normalize by Pi (Eqn. 16).
    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    answer /= M_PI_2; /* Form factor P(q) */

    // Multiply by contrast^2 and volume^2
    // volume of octahedron
    const double volume = (4./3.) * length_a * length_b * length_c * (1. - 3 * (1. - t) * (1. - t) * (1. - t));
    answer *= square((sld - solvent_sld) * volume);

    // Convert from [1e-12 A-1] to [cm-1]
    answer *= 1.0e-4;

    return answer;
}

static void
Fq(double q,
    double *F1,
    double *F2,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;

    /* Integration limits to use in Gaussian quadrature */
    const double v1a = 0.0;
    const double v1b = M_PI_2;  /* theta integration limits */
    const double v2a = 0.0;
    const double v2b = M_PI_2;  /* phi integration limits */

    double outer_sum_F1 = 0.0;
    double outer_sum_F2 = 0.0;

    for (int i = 0; i < GAUSS_N; i++) {
        const double theta = 0.5 * (GAUSS_Z[i] * (v1b - v1a) + v1a + v1b);
        double sin_theta, cos_theta;
        SINCOS(theta, sin_theta, cos_theta);

        double inner_sum_F1 = 0.0;
        double inner_sum_F2 = 0.0;
        for (int j = 0; j < GAUSS_N; j++) {
            double phi = 0.5 * (GAUSS_Z[j] * (v2b - v2a) + v2a + v2b);
            double sin_phi, cos_phi;
            SINCOS(phi, sin_phi, cos_phi);

            /* HERE: Octahedron formula */
            const double Qx = q * sin_theta * cos_phi;
            const double Qy = q * sin_theta * sin_phi;
            const double Qz = q * cos_theta;
            const double qx = Qx * length_a;
            const double qy = Qy * length_b;
            const double qz = Qz * length_c;
            const double AA = 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qx)*sin(qy*(1.-t)-qx*t)+(qy+qx)*sin(qy*(1.-t)+qx*t))+
                                1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qx)*sin(qz*(1.-t)-qx*t)+(qz+qx)*sin(qz*(1.-t)+qx*t));

            const double BB = 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qy)*sin(qz*(1.-t)-qy*t)+(qz+qy)*sin(qz*(1.-t)+qy*t))+
                                1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qy)*sin(qx*(1.-t)-qy*t)+(qx+qy)*sin(qx*(1.-t)+qy*t));

            const double CC = 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t))+
                                1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));

            /* normalisation to 1. of AP at q = 0. Division by a Factor 4/3. */
            const double AP = 6. / (1. - 3 * (1. - t) * (1. - t) * (1. - t)) * (AA + BB + CC);

            inner_sum_F1 += GAUSS_W[j] * AP;
            inner_sum_F2 += GAUSS_W[j] * AP * AP;

        }
        inner_sum_F1 = 0.5 * (v2b - v2a) * inner_sum_F1;
        inner_sum_F2 = 0.5 * (v2b - v2a) * inner_sum_F2;
        outer_sum_F1 += GAUSS_W[i] * inner_sum_F1 * sin_theta;
        outer_sum_F2 += GAUSS_W[i] * inner_sum_F2 * sin_theta;
    }

    outer_sum_F1 *= 0.5 * (v1b - v1a);
    outer_sum_F2 *= 0.5 * (v1b - v1a);

    // Normalize by Pi (Eqn. 16).
    // The factor 2 appears because the theta integral has been defined between
    // 0 and pi/2, instead of 0 to pi.
    outer_sum_F1 /= M_PI_2;
    outer_sum_F2 /= M_PI_2;

    // Multiply by contrast and volume
    // volume of octahedron
    const double s = (sld - solvent_sld) * (4./3.) * (length_a * length_b * length_c) * (1. - 3 * (1. - t) * (1. - t) * (1. - t));

    /* Convert from [1e-12 A-1] to [cm-1] and account for SLD units (1e-6/Ang^2) */
    *F1 = 1e-2 * s * outer_sum_F1;
    *F2 = 1e-4 * s * s * outer_sum_F2;
}


static double
Iqabc(double qa, double qb, double qc,
    double sld,
    double solvent_sld,
    double length_a,
    double b2a_ratio,
    double c2a_ratio,
    double t)
{
    const double length_b = length_a * b2a_ratio;
    const double length_c = length_a * c2a_ratio;

    /* HERE: Octahedron formula */
    const double qx = qa * length_a;
    const double qy = qb * length_b;
    const double qz = qc * length_c;
    const double AA = 1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qx)*sin(qy*(1.-t)-qx*t)+(qy+qx)*sin(qy*(1.-t)+qx*t))+
                                1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qx)*sin(qz*(1.-t)-qx*t)+(qz+qx)*sin(qz*(1.-t)+qx*t));

    const double BB = 1./(2*(qz*qz-qx*qx)*(qz*qz-qy*qy))*((qz-qy)*sin(qz*(1.-t)-qy*t)+(qz+qy)*sin(qz*(1.-t)+qy*t))+
                                1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qy)*sin(qx*(1.-t)-qy*t)+(qx+qy)*sin(qx*(1.-t)+qy*t));

    const double CC = 1./(2*(qx*qx-qy*qy)*(qx*qx-qz*qz))*((qx-qz)*sin(qx*(1.-t)-qz*t)+(qx+qz)*sin(qx*(1.-t)+qz*t))+
                                1./(2*(qy*qy-qz*qz)*(qy*qy-qx*qx))*((qy-qz)*sin(qy*(1.-t)-qz*t)+(qy+qz)*sin(qy*(1.-t)+qz*t));

    /* normalisation to 1. of AP at q = 0. Division by a Factor 4/3. */
    const double AP = 6. / (1. - 3 * (1. - t) * (1. - t) * (1. - t)) * (AA + BB + CC);

    /* Multiply by contrast and volume */
    const double s = (sld - solvent_sld) * (4./3.) * (length_a * length_b * length_c) * (1. - 3 * (1. - t) * (1. - t) * (1. - t));

    /* Convert from [1e-12 A-1] to [cm-1] */
    return 1.0e-4 * square(s * AP);
}
