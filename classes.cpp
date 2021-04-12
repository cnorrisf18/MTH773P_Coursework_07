#include "classes.h"

#include <cmath>

#include <vector>
using namespace std;
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


double UpAndOutEuropeanCall::BSAnalyticalPrice(double S0, double r, double sigma)
{
    const double d1 = D((S0 / m_K), r, sigma, m_T);  const double d2 = d1 - (sigma * sqrt(m_T));
    const double d3 = D((S0 / m_B), r, sigma, m_T);  const double d4 = d3 - (sigma * sqrt(m_T));
    const double d5 = D((S0 / m_B), -r, sigma, m_T);  const double d6 = d5 - (sigma * sqrt(m_T));
    const double d7 = D(((S0 * m_K) / (m_B * m_B)), -r, sigma, m_T);  const double d8 = d7 - (sigma * sqrt(m_T));

    return S0 * ((N(d1) - N(d3) - pow((m_B / S0), ((1.0 + ((2.0 * r)) / pow(sigma, 2.0)))) * (N(d6) - N(d8)))) -
        ((m_K * exp(-r * m_T)) * (N(d2) - N(d4) - pow((m_B / S0), ((-1.0 + ((2.0 * r)) / pow(sigma, 2.0)))) * (N(d5) - N(d7))));

}

double UpAndOutEuropeanCall::payoff(double S)
{
    if (S > m_B || m_K > m_B) //if either value is greater than the barrier, the option is worthless
    {
        return 0;
    }
    else if (S - m_K > 0) //basic call payoff function: max(S-K, 0)
    {
        return S - m_K;
    }
    else
    {
        return 0;
    }
}
double UpAndOutEuropeanCall::BinomialTreePrice(double S0, double r, double sigma)
{
    const int N = 90; //90 seems to return the closest to actual price for our example
    const double time_step = m_T / N;
    const double u = exp((r - (pow(sigma, 2.0) / 2.0) * time_step + sigma * sqrt(time_step)));
    const double d = exp((r - pow(sigma, 2.0) / 2.0) * time_step - sigma * sqrt(time_step));
    const double R = exp(r * time_step) - 1.0;
    const double qu = (1 + R - d) / (u - d);
    const double qd = 1 - qu;

    vector<double> v(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        const double ST = S0 * pow(u, i) * pow(d, N - i);
        v[i] = payoff(ST);
    }
    for (int n = N - 1; n >= 0; --n)
    {
        for (int i = 0; i <= n; ++i)
        {
            v[i] = (qu * v[i + 1] + qd * v[i]) / (1.0 + R);
        }
    }

    return v[0];
}

double UpAndOutEuropeanCall::FiniteDifferenceExplicitPrice(double S0, double r, double sigma)
{
	const int i_max = 36000;
	const int j_max = 800;
    //the program becomes unstable after j_max = 905, which would give C a value of ~ .039070846. The minimum value for i_max that satisfies S = 1000 steps
    //is 38950, which returns 4.94847 for the price

	//Specific to barrier call options
	const double S_max = 300;
	const double S_min = 0;

	const double dt = m_T / i_max;
	const double dS = (S_max - S_min) / j_max;

	vector<double> v(j_max + 1);

	double S = S_min;
	for (int j = 0; j <= j_max; ++j)
	{
		v[j] = payoff(S);
		S += dS;
	}

	
	double t = m_T;
	for (int i = i_max; i > 0; --i)
	{
		vector<double> u(j_max + 1);
		double S = S_min;
		for (int j = 1; j < j_max; ++j)
		{
			S += dS;

			const double a = -sigma * sigma * S * S / 2.0;
			const double b = -r * S;
			const double c = r;

			const double A = dt / dS * (b / 2.0 - a / dS);
			const double B = 1.0 - dt * c + 2.0 * dt * a / dS / dS;
			const double C = -dt / dS * (b / 2.0 + a / dS);

			u[j] = A * v[j - 1] + B * v[j] + C * v[j + 1];
		}

		t -= dt;

		//
		//	barrier call boundary conditions
		//
		u[0] = 0.0;
		u[j_max] = S <= m_B ? S_max - m_K * exp(-r * (m_T - t)) : 0.0;

		v = u;
	}
	const int left = floor(S0 - S_min) / dS;
	const double w = (S0 - (S_min + left * dS)) / dS;
	const double price = v[left] * (1.0 - w) + v[left + 1] * w;

    return price;
}



