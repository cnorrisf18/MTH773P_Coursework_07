#pragma once
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