// heat1D_explicit.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"

class Heat1DExplicit : public bSolver1DBase<double>
{
public:
    Heat1DExplicit(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                   { return 0.; }) : bSolver1DBase<double>(nx, dt, U0)
    {
        std::cout << "Starting Heat1DExplicit solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    void Step() override
    {
        Data[0][nt + 1] = LBC_value; // 下端の境界条件
        Data[Data.size() - 1][nt + 1] = RBC_value; // 上端の境界条件
        for (size_t i = 1; i < Data.size() - 1; ++i)
            Data[i][nt + 1] = (1.-2.*DN)* Data[i][nt] + DN * (Data[i - 1][nt] + Data[i + 1][nt]);
        time += dt;
        nt++;
    }
};
// heat1D_explicit.hpp