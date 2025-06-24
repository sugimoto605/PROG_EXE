// heat1D_implicit.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"
#include <Eigen/Sparse>
class Heat1DImplicit : public bSolver1DBase<double>
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
public:
    Heat1DImplicit(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                       { return 0.; }) : bSolver1DBase<double>(nx, dt, U0)
    {
        std::cout << "Starting Heat1DImplicit solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        if (parm)
        {
            K = *static_cast<double *>(parm); // 拡散係数の初期化
            DN = K * dt / dx / dx;            // 拡散数の再計算計算
        }
        CFL=0;
        Eigen::SparseMatrix<double> coef; // スパース行列を使用して係数行列を表現
        coef.resize(Data.size(), Data.size());
        for (int i = 0; i < Data.size(); i++)
        {
            if ((i == 0) || (i == Data.size() - 1)) // 境界条件の処理
                coef.insert(i, i) = 1.;              // 境界条件の行列成分を設定
            else
            {
                coef.insert(i, i - 1) = -(DN + 0.5 * CFL); // 左隣の成分
                coef.insert(i, i) = 1. + 2. * DN;          // 対角成分
                coef.insert(i, i + 1) = -(DN - 0.5 * CFL); // 右隣の成分
            }
        }
        coef.makeCompressed(); // スパース行列を圧縮
        solver.compute(coef);
        std::cout << "Heat1DImplicit solver initialized with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    // 時間ステップを実行するメソッド
    void Step() override
    {
        Eigen::VectorXd b(Data.size());
        for (int i = 0; i < Data.size(); i++)
            if (i == 0)                     b[i] = LBC_value; // 下端の境界条件
            else if (i == Data.size() - 1)  b[i] = RBC_value; // 上端の境界条件
            else                            b[i] = Data[i][nt]; // 現在の値を設定
        auto result=solver.solve(b); // 連立方程式を解く
        time = (++nt) * dt;
        for (int i = 0; i < Data.size(); i++) Data[i][nt] = result[i];
    }
};
// heat1D_implicit.hpp