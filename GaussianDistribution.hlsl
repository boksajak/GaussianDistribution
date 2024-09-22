#define SQRT_TWO 1.41421356237
#define TWO_PI 6.28318530718

// Fast inverse error function approximation
// Source: "A handy approximation for the error function and its inverse" by Sergei Winitzki.
// Code: https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
float erfinv(float x)
{
    float tt1, tt2, lnx, sgn;
    sgn = (x < 0.0f) ? -1.0f : 1.0f;
 
    x = (1.0f - x) * (1.0f + x);
    lnx = log(x);
 
    tt1 = 2.0f / (3.14159265359f * 0.147f) + 0.5f * lnx;
    tt2 = 1.0f / (0.147f) * lnx;
 
    return (sgn * sqrt(-tt1 + sqrt(tt1 * tt1 - tt2)));
}

// Precise inverse error function approximation using polynomials
// Source: "Approximating the erfinv function"
// https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
float MBG_erfinv(float x)
{
    float w, p;
    w = -log((1.0f - x) * (1.0f + x));
    if (w < 5.000000f)
    {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    }
    else
    {
        w = sqrt(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

// Generates a normally distributed sample using a CDF inversion method with fast inverse error function approximation
// Parameter u is a random number in the range [0;1)
float sampleNormalDistributionInvCDFFast(float u, float mean, float standardDeviation)
{
    return mean + SQRT_TWO * standardDeviation * erfinv(2.0f * u - 1.0f);
}

// Generates a normally distributed sample using a CDF inversion method with precise inverse error function approximation
// Parameter u is a random number in the range [0;1)
float sampleNormalDistributionInvCDFPrecise(float u, float mean, float standardDeviation)
{
    return mean + SQRT_TWO * standardDeviation * MBG_erfinv(2.0f * u - 1.0f);
}

// Generates two normally distributed samples using a Box-Muller method
// Parameter u are two random numbers in the range [0;1)
float2 sampleNormalDistributionBoxMuller(float2 u, float mean, float standardDeviation)
{
    const float a = standardDeviation * sqrt(-2.0f * log(1.0f - u.x));
    const float b = TWO_PI * u.y;
    return float2(cos(b), sin(b)) * a + mean;
}