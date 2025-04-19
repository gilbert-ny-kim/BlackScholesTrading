#include <stdio.h>
#include <stdlib.h>
#include <mysql.h>
#include <math.h>
#include <time.h>


#define DAYS 252  // Number of trading days in a year


// Function prototypes
double ln(double x);
double sqr(double x);
double cdf(double x);
double black_scholes_call(float S, float K, float sigma, float r, float t);
double brownian_motion(double S, double mu, double sigma, double dt, double path[]);
double box_muller();

int main() {
    // === Simulation Parameters ===
    const double S0 = 100.0;      // Initial stock price
    const double mu = 0.05;       // Expected return
    const double sigma = 0.3;     // Volatility
    const double dt = 1.0 / DAYS; // Time step (1 trading day)
    double path[DAYS];            // Array to store simulated path

    // === Random Seed Setup ===
    srand(time(NULL));      // Seed random number generator

    // === Generate GBM Path ===
    brownian_motion(S0, mu, sigma, dt, path);

    // === Black-Scholes Pricing ===
    const double K = 100.0;   // Strike price
    const double r = 0.03;    // Risk-free interest rate
    const double T = 1.0;     // Time to maturity (1 year)

    double call_price = black_scholes_call(S0, K, sigma, r, T);
    printf("Black-Scholes Call Option Price: %.5f\n", call_price);

    // === MySQL Connection ===
    MYSQL *conn = mysql_init(NULL);
    if (!conn || !mysql_real_connect(conn, "localhost", "gilbertnykim", "wnsdn08@Lotte",
                                     "BlackScholesTrading", 0, NULL, 0)) {
        fprintf(stderr, "MySQL connection failed: %s\n", mysql_error(conn));
        if (conn) mysql_close(conn);
        return EXIT_FAILURE;
    }

    // === Clear Table ===
    if (mysql_query(conn, "DELETE FROM stock_path")) {
        fprintf(stderr, "Failed to clear table: %s\n", mysql_error(conn));
        mysql_close(conn);
        return EXIT_FAILURE;
    }

    // === Insert Path Data ===
    char query[256];
    for (int i = 0; i < DAYS; i++) {
        snprintf(query, sizeof(query),
                 "INSERT INTO stock_path (day, price) VALUES (%d, %.5f);", i, path[i]);

        if (mysql_query(conn, query)) {
            fprintf(stderr, "INSERT failed (day %d): %s\n", i, mysql_error(conn));
        }
    }

    mysql_close(conn);

    return 0;
}

// Simple math function
double ln(double x) {
    return log(x);
}

double sqr(double x) {
    return x * x;
}

double cdf(double x) {
    return 0.5 * (1 + erf(x * M_SQRT1_2));
}

// Box-Muller Transform
double box_muller() {
    double U1 = (double)rand() / RAND_MAX;  // Uniform (0,1]
    double U2 = (double)rand() / RAND_MAX;

    return sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);  // Z ~ N(0,1)
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

// Brownian Motion Simulation
// S = Current stock price
// mu = Excepted return
// sigma = Volatility
// dt = Time step
double brownian_motion(double S, double mu, double sigma, double dt, double path[]) {
    path[0] = S;
    
    for (int i = 1; i < DAYS; i++) {
        double Z = box_muller(); // Generate a standard normal random variable
        double dW = sqrt(dt) * Z;

        S *= exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);

        path[i] = S;
    }
    
    return S;
}
