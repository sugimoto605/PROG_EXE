// pSolverFVBase.hpp    Periodic soltuion by FV method using Cellwise Linear Reconsruction
#pragma once
#include "pSolver1DBase.hpp"
class FVKit
{
public:
    double U;       //セル内平均
    double dU;      //セル内勾配 ∂U/∂x
    double RFlux;   //右境界フラックス/c
    FVKit &operator=(const double &x)
    {
        U = x;
        RFlux = dU = 0;
        return *this;
    }
    operator double() const { return U; }
};
class pSolverFVBase : public pSolver1DBase<FVKit>
{
protected:
    void Grad2Flux()    // Calculate Numerical Flux by Linear Reconstruction
    {
        for (auto &P : Data) P[nt].RFlux=P[nt].U+0.5*CFL*(1.-CFL)*P[nt].dU*dx;
    }
public:
    pSolverFVBase(const size_t &nx, const double &dt, const std::function<double(double)> &U0) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverFVBase..." << std::endl;
        //  格子点ではなくセルなので1個少ない
        Data.resize(Data.size()-1);
        Data.front().pre=&Data.back();
        Data.back().next=&Data.front();
        for (auto &P : Data) P.X += 0.5 * dx;
    }
    void Initialize(void *parm=nullptr) override
    {
        std::cout << "pSolverFVBase::Initialize()" << std::endl;
    }
    virtual void getFlux()=0;
    void Step(double U_W = 0.0) override
    {
        getFlux();
        time = (++nt) * dt;
        for (auto &P : Data)
            P[nt].U = P[nt - 1].U - CFL * (P[nt - 1].RFlux - (*P.pre)[nt - 1].RFlux);
    }
};