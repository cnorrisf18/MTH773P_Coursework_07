#pragma once

class UpAndOutEuropeanCall {
public:
	UpAndOutEuropeanCall(double T, double K, double B) : m_T(T), m_K(K), m_B(B) {};
	double payoff(double S);
	double BSAnalyticalPrice(double S0, double r, double sigma);
	double BinomialTreePrice(double S0, double r, double sigma);
	double FiniteDifferenceExplicitPrice(double S0, double r, double sigma);
private:
	const double m_T, m_K, m_B;
};
