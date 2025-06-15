// pSolverFVLW.hpp    Lax-Wendroff periodic Finite Volume
#pragma once
#include "../include/pSolverFVBase.hpp"
class pSolverFVLW : public pSolverFVBase
{
public:
    pSolverFVLW(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                    { return 0.; }) : pSolverFVBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverFVLW..." << std::endl;
    }
    virtual void getFlux() override
    {
        for(auto& P:Data) P[nt].dU =((*P.next)[nt].U-P[nt].U)/dx;
        Grad2Flux();
    }
};
