#include <stdio.h>
#include <math.h>


// Function prototypes
double ln(float x);
double sqr(float x);
double cdf(double x);
double black_scholes_call(float S, float K, float sigma, float r, float t);

int main() {
    double call_option_price = black_scholes_call(100, 100, .3, .05, 1);

    printf("Black-Scholes Call Option Price: %f\n", call_option_price);

    return 0;
}

// Simple math function
double ln(float x) {
    return log(x);
}

double sqr(float x) {
    return x * x;
}

double cdf(double x) {
    return 0.5 * (1 + erf(x * M_SQRT1_2));
}

// Black-Scholes Model
// C = Call option price
// S = Current stock (or other underlying) price
// K = Strike price
// r = Risk-free interest rate
// t = Time to maturity
// N = A normal distribution
double black_scholes_call(float S, float K, float sigma, float r, float t) {
    double d1 = (ln(S / K) + ((r + (sqr(sigma) / 2)) * t)) / (sigma * sqrt(t));
    double d2 = d1 - (sigma * sqrt(t));

    double C = S * cdf(d1) - exp(-r * t) * K * cdf(d2);

    return C;
}