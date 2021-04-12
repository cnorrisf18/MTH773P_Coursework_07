#include <iostream>
#include "classes.h"

using namespace std;


/*
int main()
{
	const double K = 100;
	const double T = 1.0;
	const double B = 130; 
	const double S0 = 98;
	const double r = 0.05;
	const double sigma = 0.2;
	
	UpAndOutEuropeanCall SampleCall(T, K, B);
	
	//cout << SampleCall.FiniteDifferenceExplicitPrice(S0, r, sigma) << "\n" << SampleCall.BSAnalyticalPrice(S0, r, sigma);

	ChooserOption SampleChooser(T, K, 0.5);
	//base option price is 13.8513

	for (int t = 0; t < 101; ++t)
	{
		cout << ChooserOption(T, K, double(t) / 100).BSAnalyticalPrice(100, r, sigma) << endl;
		//The price grows higher as t gets closer to T, which makes sense because you're more likely to know whether you want the option to be call or put as you get closer to T. 
	}
	//cout << SampleChooser.BSAnalyticalPrice(100.0, r, sigma);
}
*/


#include <iostream>
#include "BSModel01.h"
#include "Option.h"
#include "BSEq.h"
#include "CNMethod.h"
#include "ChooserOption.h"

int main()
{
	double S0 = 100.0, r = 0.05, sigma = 0.2;
	BSModel Model(S0, r, sigma);

	double T = 1. / 12., K = 100.0, zl = 0.0, zu = 2.0 * S0;
	Put EuropeanPut(K, T, zl, zu);

	int imax = 200, jmax = 2000;

	BSEq BSPDE(&Model, &EuropeanPut);

	CNMethod Method(&BSPDE, imax, jmax);
	Method.SolvePDE();
	cout << "Price = " << Method.v(0.0, S0) << endl;

	return 0;
}
