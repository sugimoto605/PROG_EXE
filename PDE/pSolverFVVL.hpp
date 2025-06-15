// pSolverFVVL.hpp    TVD periodic Finite Volume with van-Leer MC-Limitter
#pragma once
#include "../include/pSolverFVBase.hpp"
class pSolverFVVL : public pSolverFVBase
{
    double MinMod(const double a, const double b)
    {
        if (a * b < 0)
            return 0.;
        if (std::abs(a) < std::abs(b))
            return a;
        return b;
    }
    double MinMod(const double a, const double b,const double c)
    {
        return MinMod(MinMod(a, b),c);
    }
public:
    pSolverFVVL(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                        { return 0.; }) : pSolverFVBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverFVVL..." << std::endl;
    }
    virtual void getFlux() override
    {
        for (auto &P : Data)
        {
            auto &L = *P.pre;
            auto &R = *P.next;
            P[nt].dU = MinMod(2.*(R[nt].U - P[nt].U),2.*(P[nt].U - L[nt].U),.5*(R[nt].U-L[nt].U)) / dx;
        }
        Grad2Flux();
    }
};
