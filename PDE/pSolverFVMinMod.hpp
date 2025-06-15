// pSolverFVMinMod.hpp    TVD periodic Finite Volume with MinMod Gradient
#pragma once
#include "../include/pSolverFVBase.hpp"
class pSolverFVMinMod : public pSolverFVBase
{
    double MinMod(const double a,const double b)
    {
        if (a*b<0) return 0.;
        if (std::abs(a)<std::abs(b)) return a;
        return b;
    }
public:
    pSolverFVMinMod(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                    { return 0.; }) : pSolverFVBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverFVMinMod..." << std::endl;
    }
    virtual void getFlux() override
    {
        for (auto &P : Data)
        {
            auto& L=*P.pre;
            auto& R=*P.next;
            P[nt].dU = MinMod(P[nt].U - L[nt].U,R[nt].U-P[nt].U)/dx;
        }
        Grad2Flux();
    }
};
