#include "ChooserOption.h"
#include <cmath>

double CumulativeIntegral(double x)
{
    return exp(-pow(x, 2) / 2);
}
double Pricing_Integral(double a, double b, double N)
{
    double h = (b - a) / N;
    double sum = 0;
    for (double i = 0; i < N + 1; ++i)
    {
        double multiplier = 2.0;
        if (i == 0 || i == N)
        {
            multiplier = 1.0;
        }
        sum += multiplier * CumulativeIntegral((a + (i * h)));
    }
    return (h / 2.0) * sum;
}
double N(double z) { //Cumulative Normal Distribution function

    const double pi = 4.0 * atan(1.0);
    double integral_value = Pricing_Integral(0, z, 1000.0);
    return integral_value / sqrt(2.0 * pi);
}
double D(double x, double r, double sigma, double T)
{
    return (log(x) + (r + ((.5 * pow(sigma, 2)) * T))) / (sigma * sqrt(T));
}






double ChooserOption::BSAnalyticalPrice(double S0, double r, double sigma)
{
    const double d = (log(S0 / m_K) + (r + .5 * pow(sigma, 2)) * m_T) / (sigma * sqrt(m_T));
    const double y = (log(S0 / m_K) + (r * m_T) + (.5 * pow(sigma, 2) * m_Tc)) / (sigma * sqrt(m_Tc));

    const double price = S0 * N(d) - m_K * exp(-r * m_T) * N(d - sigma * sqrt(m_T)) - S0 * N(-y) + m_K * exp(-r * m_T) * N(-y + sigma * sqrt(m_Tc));
    return price;
}

