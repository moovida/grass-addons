/*
 * random
 *
 *	RNG for landslide size simulator
 *
 * Copyright (c) 2004 by Colin P. Stark
 *
 * $Header: /local/home/meander/cstark/work/slides_dev/model/random.c,v 1.1
 *2004/11/26 18:10:55 cstark Exp $
 *
 * $Log: random.c,v $
 * Revision 1.1  2004/11/26 18:10:55  cstark
 * Initial revision
 *
 */

static char rcsId[] = "$Id: random.c,v 1.1 2004/11/26 18:10:55 cstark Exp $";
#define VERSION 1.1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "random2.h"

#define TRUE    -1
#define FALSE   0

#define GCC_RNG TRUE
/* #define CC_RNG TRUE */

static long MarsagliaRandom(UniSave* uniData);
static void MarsagliaSeed(UniSave* uniData, unsigned int);

static long random_irpi(UniSave* uniData);
static void srandom_irpi(UniSave* uniData, unsigned int);

static long (*randomProc)(UniSave*) = MarsagliaRandom;
static void (*seedProc)(UniSave*, unsigned int) = MarsagliaSeed;



static long MarsagliaRandom(UniSave* uniData)
{
    float luni;
    luni = uniData->u[uniData->ui] - uniData->u[uniData->uj];
    if (luni < 0.0) {
        luni += 1.0;
    }
    uniData->u[uniData->ui] = luni;
    if (--uniData->ui == 0) {
        uniData->ui = 97;
    }
    if (--uniData->uj == 0) {
        uniData->uj = 97;
    }
    if ((uniData->c -= uniData->cd) < 0.0) {
        uniData->c += uniData->cm;
    }
    if ((luni -= uniData->c) < 0.0) {
        luni += 1.0;
    }
    return (long)((luni * (1 << 24)));
}

static void Rstart(UniSave* uniData, int i, int j, int k, int l)
{
    int ii, jj, m;
    float s, t;
    for (ii = 1; ii <= 97; ii++) {
        s = 0.0;
        t = 0.5;
        for (jj = 1; jj <= 24; jj++) {
            m = ((i * j % 179) * k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53 * l + 1) % 169;
            if (l * m % 64 >= 32) {
                s += t;
            }
            t *= 0.5;
        }
        uniData->u[ii] = s;
    }
    uniData->c = 362436.0 / 16777216.0;
    uniData->cd = 7654321.0 / 16777216.0;
    uniData->cm = 16777213.0 / 16777216.0;
    uniData->ui = 97; /*  There is a bug in the original Fortran version */
    uniData->uj = 33; /*  of UNI -- i and j should be SAVEd in UNI()     */
}

static void MarsagliaSeed(UniSave* uniData, unsigned int ijkl)
{
    int i, j, k, l, ij, kl;
    if (ijkl > 900000000) {
        fprintf(stderr, "SeedUni: ijkl = %d -- out of range\n\n", ijkl);
    }

    ij = ijkl / 30082;
    kl = ijkl - (30082 * ij);

    i = ((ij / 177) % 177) + 2;
    j = (ij % 177) + 2;
    k = ((kl / 169) % 178) + 1;
    l = kl % 169;

    if ((i <= 0) || (i > 178)) {
        fprintf(stderr, "SeedUni: i = %d -- out of range\n\n", i);
    }
    if ((j <= 0) || (j > 178)) {
        fprintf(stderr, "SeedUni: j = %d -- out of range\n\n", j);
    }
    if ((k <= 0) || (k > 178)) {
        fprintf(stderr, "SeedUni: k = %d -- out of range\n\n", k);
    }
    if ((l < 0) || (l > 168)) {
        fprintf(stderr, "SeedUni: l = %d -- out of range\n\n", l);
    }
    if (i == 1 && j == 1 && k == 1) {
        fprintf(stderr, "SeedUni: 1 1 1 not allowed for 1st 3 seeds");
    }
    Rstart(uniData, i, j, k, l);
}

double Uniform(UniSave* uniData, double mean, double std_devn)
{
    /*
     * Return a uniform random variate of zero mean and unit variance
     *   - rember that variance of a uniform pdf is  (b-a)^2/12
     */
    std_devn *= sqrt(12.0);
    return (mean - std_devn / 2) +
           std_devn * (double)randomProc(uniData) / ((double)uniData->randMax + 1);
}

double SimpleUniform(UniSave* uniData, double min, double max)
{
    return min + (max - min) * (double)randomProc(uniData) / ((double)uniData->randMax + 1);
}

double Gaussian(UniSave* uniData, double mean, double std_devn)
{
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;
    if (iset == 0) {
        do {
            v1 = 2.0 * ((double)(randomProc(uniData) & uniData->randMax) / (double)uniData->randMax) -
                 1.0;
            v2 = 2.0 * ((double)(randomProc(uniData) & uniData->randMax) / (double)uniData->randMax) -
                 1.0;
            r = v1 * v1 + v2 * v2;
        } while (r >= 1.0);
        fac = sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac * std_devn + mean;
    }
    else {
        iset = 0;
        return gset * std_devn + mean;
    }
}

#ifdef M_PI
#else
#define M_PI (3.141592653589793115997963468544185161590576171875)
#endif

double Cauchy(UniSave* uniData, double mean, double half_width)
{
    return (mean + half_width * tan(SimpleUniform(uniData, -M_PI / 2, M_PI / 2)));
}

/*
 *----------------------------------------------------------------------
 *
 * SymmStable --
 *
 *
   Stable Random Number Generators (McCulloch, Aug96)

   Written and submitted by J. Huston McCulloch for public,
   non-commercial use.
   Contact author at mcculloch.2@osu.edu; 614-292-0382.

   This file contains two stable random number generator subroutines
     RNDSSTA -- generates vector of standard symmetric stable random
               numbers with alpha in [.1,2] and beta = 0.
     RNDSTA -- generates vector of standard stable random numbers with
               alpha in [.1,2] and arbitrary beta in [-1,1].
   At the head of the file is a demonstration driver program that
    calls RNDSTA, lets you type in alpha and beta, and makes nice
    graphs and histograms of the results.  Values of alpha in [1.5, 1.9],
    with beta in [-.3, .3], simulate many financial asset returns.
    Values of alpha in [1, 1.2] with beta = -1 often simulate icicles.
   See the comments at the beginning of RNDSTA and RNDSSTA for details
    of these procs.
   See McCulloch "Financial Applications of Stable Distributions",
    forthcoming 1996 in Volume 14 of the _Handbook of Statistics_ for
    a survey of pertinent literature.


   Returns rX1 vector of iid standard symmetric stable pseudo-random
   numbers with characteristic exponent alpha in [.1,2], and skewness
   parameter beta = 0, using method of Chambers, Mallows and Stuck
   (JASA 1976, 340-4).
   Encoded in GAUSS by J. Huston McCulloch, Ohio State University Econ Dept.
   (mcculloch.2@osu.edu), 12/95, directly from the CMS equation (4.1),
   as simplified when beta = 0.
   Each r.v. has log characteristic function
        log E exp(ixt) = -abs(t)^alpha
   When alpha = 2, the distribution is Gaussian, with variance 2.
   When alpha = 1, the distribution is standard Cauchy.
   If alpha > 1, Ex = 0.  In all cases, the median is 0.
   This proc uses 2r uniform pseudo-random numbers from RNDU.

 * Results:
 *
 *
 * Side effects:
 *
 *
 *----------------------------------------------------------------------
 */

double SymmStable(UniSave* uniData, double mean, double alpha, double half_width)
{
    double phi, w;
    w = -log(SimpleUniform(uniData, DBL_MIN, 1));
    phi = (SimpleUniform(uniData, DBL_MIN, 1) - 0.5) * M_PI;
    if (alpha == 2.0) {
        return mean + (2.0 * sqrt(w) * sin(phi));
    }
    else if (alpha == 1.0) {
        return mean + (tan(phi));
    }
    return mean +
           half_width *
               (pow((cos((1.0 - alpha) * phi) / w), (1.0 / alpha - 1.0)) *
                sin(alpha * phi) / pow(cos(phi), (1.0 / alpha)));
}

double Holtsmark(UniSave* uniData, double mean, double half_width)
{
    return SymmStable(uniData, mean, 1.5, half_width);
}

/*
 *----------------------------------------------------------------------
 *
 * Stable --
 *
 *
   Stable Random Number Generators (McCulloch, Aug96)

   Written and submitted by J. Huston McCulloch for public,
   non-commercial use.
   Contact author at mcculloch.2@osu.edu; 614-292-0382.

   This file contains two stable random number generator subroutines
     RNDSSTA -- generates vector of standard symmetric stable random
               numbers with alpha in [.1,2] and beta = 0.
     RNDSTA -- generates vector of standard stable random numbers with
               alpha in [.1,2] and arbitrary beta in [-1,1].
   At the head of the file is a demonstration driver program that
    calls RNDSTA, lets you type in alpha and beta, and makes nice
    graphs and histograms of the results.  Values of alpha in [1.5, 1.9],
    with beta in [-.3, .3], simulate many financial asset returns.
    Values of alpha in [1, 1.2] with beta = -1 often simulate icicles.
   See the comments at the beginning of RNDSTA and RNDSSTA for details
    of these procs.
   See McCulloch "Financial Applications of Stable Distributions",
    forthcoming 1996 in Volume 14 of the _Handbook of Statistics_ for
    a survey of pertinent literature.

   Returns rX1 vector x of iid standard stable pseudo-random
   numbers with characteristic exponent alpha in [.1,2], and
   skewness parameter beta in [-1,1].  Based on the method of Chambers,
   Mallows and Stuck (JASA 1976, 340-4).  Encoded in GAUSS by
   J. Huston McCulloch, Ohio State Univ. Econ. Dept.
   (mcculloch.2@osu.edu), 12/95.
     The CMS method was applied in such a way that x will have the
   standard (delta =0, c = 1) characteristic function
    log E exp(ixt) = -abs(t)^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2)),
                               for alpha /= 1,
                   = -abs(t)*(1+i*beta*(2/pi)*sign(t)*log(abs(t))),
                               for alpha = 1.
   With this parameterization, Ex = delta = 0 when alpha > 1, so that
   if the distribution is positively skewed (beta > 0),
   med(x) must become very negative as alpha approaches 1 from above.
   For alpha < 1, delta is the "focus of stability", so that if beta > 0,
   med(x) must become very positive as alpha approaches 1 from below.
   If beta /= 0, there is therefore a discontinuity in the distribution
   as a function of alpha as alpha passes 1.
     CMS remove this discontinuity by subtracting
                zeta = beta*tan(pi*alpha/2)
   (equivalent to their -tan(alpha*phi0)) from x for alpha /= 1
   in their program RSTAB (also known as GGSTA in IMSL).
   The result is a random number whose
   distribution is a continuous function of alpha, but whose location
   parameter has no known useful interpretation other than computational
   convenience.  The present program restores the standard parameterization
   by using the CMS (4.1), but with zeta added back in (ie with their
   initial tan(alpha*phi0) deleted).
     When alpha is precisely unity, x is computed from the CMS (2.4).
   Rather than using the CMS D2 and exp2 functions to compensate for ill-
   conditioning of (4.1) when alpha is very near 1, the present program
   merely fudges these cases by computing x from (2.4) and adjusting for
   zeta when alpha is within 1.e-8 of 1.  This should make no difference
   for simulation results with samples of size less than approximately
   10^8, and then only when the desired alpha is within 1.e-8 of 1 but
   not equal to 1.
     Companion proc RNDSSTA computes symmetric stable random variables
   with beta = 0. It is faster, and the above issues do not arise in
   these cases.
     When alpha = 2, the distribution is Gaussian, with variance 2, regardless
   of beta.  When alpha = 1 and beta = 0, the distribution is standard Cauchy.
   When alpha = .5 and beta = 1, the distribution is Levy (reciprocal
   Chi-squared(1)).



 *
 * Results:
 *
 *
 * Side effects:
 *
 *
 *----------------------------------------------------------------------
 */

double Stable(UniSave* uniData, double alpha, double beta, double r)
{
    double phi, w, zeta, x, cosphi, bphi;
    /*
     *   Bounds:
     *    0.1<=alpha<1
     *    -1<=beta<=1
     */
    w = -log(SimpleUniform(uniData, DBL_MIN, 1));
    phi = (SimpleUniform(uniData, DBL_MIN, 1) - 0.5) * M_PI;
    cosphi = cos(phi);
    bphi = (M_PI / 2.0) + beta * phi;
    x = (2.0 / M_PI) *
        (bphi * tan(phi) - beta * log((M_PI / 2.0) * w * cosphi / bphi));
    zeta = beta * tan((M_PI * alpha) / 2.0);
    return r * (x + zeta);
}

double DoubleLevy(UniSave* uniData, double m, double r)
{
    return (SimpleUniform(uniData, -1, 1) >= 0 ? 1 : -1) * Stable(uniData, 0.5, -1, r) + m;
}

double Pareto(UniSave* uniData, double c, double alpha, double m)
{
    double x, y;

    /*
     * p(x) ~ U(0,1)
     */
    x = SimpleUniform(uniData, 0.0, 1.0);
    /*
     * Simple Pareto prob density:
     *
     * p(y|alpha,c) ~ alpha*c^alpha* y^(-alpha-1)
     *
     *
     * Truncated Pareto prob density (pdf):
     *
     * p(y|alpha,c,m) ~  alpha*c^alpha
     *                  --------------- y^(-alpha-1)
     *                  [1-(c/m)^alpha]
     *
     * Truncated Pareto probability (cdf):
     *
     * P(Y=y|alpha,c,m) ~ [1-(c/y)^alpha]
     *                    ---------------
     *                    [1-(c/m)^alpha]
     *
     * Transformation:
     *
     * If p(x) ~ U(0,1) and x samples P(Y=y|.)
     * then p(y) ~ dPareto(c,m)
     * given:  y(x) = c[1-x*(1-(c/m))^alpha]^(-1/alpha)
     *
     */
    y = c * pow(1.0 - x * (1.0 - pow(c / m, alpha)), -1.0 / alpha);

    return y;
}

double StdExponential(UniSave* uniData)
{
    static float q[8] = {0.6931472f, 0.9333737f, 0.9888778f, 0.9984959f,
                         0.9998293f, 0.9999833f, 0.9999986f, .9999999f};
    static long i;
    static float sexpo, a, u, ustar, umin;
    static float *q1 = q;
    a = 0.0;
    u = (float)SimpleUniform(uniData, 1e-20, 1);

    //! TODO make sure below substitute is correct
    //     goto S30;
    // S20:
    //     a += *q1;
    // S30:
    //     u *= 2.0;
    //     if (u < 1.0)
    //         goto S20;
    u *= 2.0;
    while (u < 1.0) {
        a += *q1;
        u *= 2.0;
    }

    u -= 1.0;
    if (u <= *q1) {
        sexpo = a + u;
        return sexpo;
    }
    i = 1;
    ustar = (float)SimpleUniform(uniData, 0, 1);
    umin = ustar;
    do {
        ustar = (float)SimpleUniform(uniData, 0, 1);
        if (ustar < umin) {
            umin = ustar;
        }
        i++;
    } while (u > *(q + i - 1));
    sexpo = a + umin * (*q1);
    return sexpo;
}

double Exponential(UniSave* uniData, double r)
{
    return StdExponential(uniData) * r;
}

double DoubleExponential(UniSave* uniData, double m, double r)
{
    return (SimpleUniform(uniData, -1, 1) >= 0 ? 1 : -1) * StdExponential(uniData) * r /
               sqrt(2) +
           m;
}

double HalfExponential(UniSave* uniData, double m, double r)
{
    m = 0;
    return StdExponential(uniData) * r / sqrt(2);
}

float SimpleGamma(UniSave* uniData, float a)
{
    float d, c, x, v, u;
    d = (float)(a - 1.0 / 3.0);
    c = (float)(1 / sqrt(9.0 * d));
    for (;;) {
        do {
            x = (float)Gaussian(uniData, 0, 1);
            v = 1.f + c * x;
        } while (v <= 0.);
        v = v * v * v;
        u = (float)SimpleUniform(uniData, 0, 1);
        if (u < 1. - .0331 * (x * x) * (x * x))
            return (d * v);
        if (log(u) < 0.5 * x * x + d * (1. - v + log(v)))
            return (d * v);
    }
}

double Gamma(UniSave* uniData, double mean, double stddevn)
{
    double alpha = pow(mean / stddevn, 2.0);
    double theta = (stddevn * stddevn) / mean;
    if (alpha < 1.0) {
        return theta * SimpleGamma(uniData, (float)(alpha + 1)) *
               pow(SimpleUniform(uniData, 0, 1), 1.0 / alpha);
        // fprintf(stderr,"alpha=%g, theta=%g, var=%g\n",alpha,theta,gamvar);
    }
    else {
        return mean * SimpleGamma(uniData, (float)(alpha)) / alpha;
        // fprintf(stderr,"alpha=%g, theta=%g, var=%g\n",alpha,theta,gamvar);
    }
}

void Seed_RNG(UniSave* uniData, unsigned int seed)
{
    seedProc(uniData, seed);
    return;
}

void Init_RNG(UniSave* uniData, int RNG, unsigned int seed)
{
    if (RNG == 1) {
        /*
         * For Marsaglia RNG
         */
        uniData->randMax = 16777215;
        randomProc = MarsagliaRandom;
        seedProc = MarsagliaSeed;
    }
    else {
#ifdef CC_RNG
        /*
         *  For SUNSWpro cc
         */
        randMax = 0x7fffffff;
        randomProc = random_irpi;
        seedProc = srandom_irpi;
#else
        /*
         *  For Linux gcc
         */
        uniData->randMax = RAND_MAX;
        randomProc = random_irpi;
        seedProc = srandom_irpi;
#endif
    }
    seedProc(uniData, seed);
    return;
}

static long random_irpi(UniSave* uniData)
{
    return (long)rand();
}

static void srandom_irpi(UniSave* uniData, unsigned int seed)
{
    srand(seed);
}
