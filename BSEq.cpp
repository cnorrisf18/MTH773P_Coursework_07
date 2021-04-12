#include "ParabPDE.h"
#include "BSModel01.h"
#include "Option.h"

class BSEq : public ParabPDE
{
public:
	BSEq(BSModel* PtrModel_, Option* PtrOption_);
	BSModel* PtrModel;
	Option* PtrOption;

	double a(double t, double z);
	double b(double t, double z);
	double c(double t, double z);
	double d(double t, double z);

	double f(double z);
	double fl(double t);
	double fu(double t);
};

BSEq::BSEq(BSModel* PtrModel_, Option* PtrOption_)
{
	PtrModel = PtrModel_; PtrOption = PtrOption_;
	T = PtrOption->T;
	xl = PtrOption->zl;
	xu = PtrOption->zu;
}

double BSEq::a(double t, double z)
{
	return -0.5 * pow(PtrModel->sigma * z, 2.0);
}

double BSEq::b(double t, double z)
{
	return -PtrModel->r * z;
}

double BSEq::c(double t, double z)
{
	return PtrModel->r;
}

double BSEq::d(double t, double z)
{
	return 0.0;
}

double BSEq::f(double z)
{
	return PtrOption->Payoff(z);
}

double BSEq::fl(double t)
{
	return PtrOption->LowerBdCond(PtrModel, t);
}

double BSEq::fu(double t)
{
	return PtrOption->UpperBdCond(PtrModel, t);
}

