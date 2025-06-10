// pSolverFVBase.hpp    Solve some periodic soltuion by Finite Volume method
#pragma once
#include "pSolver1DBase.hpp"
class FVKit
{
public:
    double U;       //セル内平均
    double LFlux;   //左境界フラックス
    FVKit &operator=(const double &x)
    {
        U = x;
        LFlux = 0;
        return *this;
    }
    operator double() const { return U; }
};
class pSolverFVBase : public pSolver1DBase<FVKit>
{
    double CFL2;
public:
    pSolverFVBase(const size_t &nx, const double &dt, const std::function<double(double)> &U0) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverFVBase..." << std::endl;
    }
    void Initialize(void *parm) override
    {
        CFL2=.5*CFL*(1.-CFL);
    }
    virtual void getSlope()=0;
    void Step(double U_W = 0.0) override
    {
        time = (++nt) * dt;
        for (auto &P : Data)
        {
            auto &Q = *P.pre;   // 左隣
            auto dU = P[nt - 1].U - Q[nt - 1].U;
            auto dV = P[nt-1].LFlux-Q[nt-1].LFlux;
            P[nt].U = P[nt - 1].U - CFL * dU - CFL2 * dV * dx;
        }
    }
};