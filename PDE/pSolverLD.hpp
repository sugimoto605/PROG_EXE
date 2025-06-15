// pSolverLD.hpp
#pragma once
#include "../include/pSolver1DBase.hpp"
class pSolverLD : public pSolver1DBase<double>
{
public:
    pSolverLD(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                  { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverLD..." << std::endl;
    }
    void Step(double U_W = 0.0) override
    {
        for (auto &P : Data) P[nt+1] = P[nt] - CFL * (P[nt] - (*P.pre)[nt]);
        time = (++nt) * dt;
    }
};