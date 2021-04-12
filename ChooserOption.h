#pragma once
#include "Option.h"
class ChooserOption : public Option
{
public:
	ChooserOption(double T, double K, double Tc) : m_T(T), m_K(K), m_Tc(Tc) {};
	double BSAnalyticalPrice(double S0, double r, double sigma);
	double FiniteDifferenceExplicitPrice(double S0, double r, double sigma);
	double Payoff(double S, double sigma, double r) { return FiniteDifferenceExplicitPrice(S, sigma, r) };
private:
	const double m_T, m_K, m_Tc;
};