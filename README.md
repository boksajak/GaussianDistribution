# Gaussian (Normal) Distribution Sampler in HLSL

## Example Usage

```cpp
#include "GaussianDistribution.hlsl"

// Initialize your random numbers generator
// It should generate uniformly distributed random numbers in the interval <0,1)
RngState rng = InitializeRNG(randomSeed)

// Setup the gaussian distribution parameters
float gaussianMean = 0.0f;
float gaussianStandardDeviation = 3.0f;

// OPTION 1: Use CDF inversion method with FAST erfinv() approximation
// Generate a uniform random number
float u = rand(rng);

// Generate a normally distributed sample
float sample = sampleNormalDistributionInvCDFFast(u, gaussianMean, gaussianStandardDeviation);

// OPTION 2: Use CDF inversion method with PRECISE erfinv() approximation
// Generate a uniform random number
float u = rand(rng);

// Generate a normally distributed sample
float sample = sampleNormalDistributionInvCDFPrecise(u, gaussianMean, gaussianStandardDeviation);
     
// OPTION 3: Use Box Muller method
// Generate two uniform random numbers
float2 u = float2(rand(rng), rand(rng));

// Generate two normally distributed samples
float2 samples = sampleNormalDistributionBoxMuller(u, gaussianMean, gaussianStandardDeviation);

```