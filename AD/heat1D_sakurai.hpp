// heat1D_sakurai.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"
#include "../include/libSPARSE.hpp"

class SKit
{
public:
    double U;           // U(t,x)
    double dUdx;        // ∂u/∂x
    double dUdx2;      // ∂²u/∂x²
    double dUdx3;      // ∂³u/∂x³
    double dUdt;        // ∂u/∂t
    double dUdt2;       // ∂²u/∂t²
    SKit &operator=(const double &x)
    {
        U=x;
        dUdx=dUdx2=dUdx3=dUdt=dUdt2=0;
    }
    operator double() const { return U; }
};

class Heat1D_Sakurai : public bSolver1DBase<SKit>
{
    SPARSE::Coef coef;
public:
};
