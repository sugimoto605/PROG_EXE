// pSolverHD.hpp
#pragma once
#include "../include/pSolver1DBase.hpp"
class pSolverHD : public pSolver1DBase<double>
{
public:
    pSolverHD(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                  { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverHD..." << std::endl;
    }
    void Step() override
    {
        time = (++nt) * dt;
        for (auto &P : Data)
        {
            auto &Pn1 = *P.pre;   // 左隣
            auto &Pn2 = *Pn1.pre; // 左左隣
            P[nt] = (1. - 1.5 * CFL) * P[nt - 1] + 2. * CFL * Pn1[nt - 1] - 0.5 * CFL * Pn2[nt - 1];
        }
    }
};