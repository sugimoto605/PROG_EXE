// heat1D_implicit.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"
#include <Eigen/Sparse>
class Heat1DImplicit : public bSolver1DBase<double>
{
    Eigen::SparseMatrix<double> coef; // スパース行列を使用して係数行列を表現
public:
    Heat1DImplicit(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                       { return 0.; }) : bSolver1DBase<double>(nx, dt, U0)
    {
        std::cout << "Starting Heat1DImplicit solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        K = *static_cast<double*>(parm); // 拡散係数の初期化
        DN = K * dt / dx / dx; // 拡散数の再計算計算
        for (int i = 0; i < Data.size(); i++) coef.New(i, 0.0);
        for (int i = 0; i < Data.size(); i++)
        {
            if (i==0)
            {
                coef[i][i] = 1;
            }
            else if (i == Data.size() - 1)
            {
                coef[i][i] = 1;
            }
            else
            {
                coef[i][i - 1] = -(DN+0.5*CFL);       // 左隣の成分
                coef[i][i] = 1. + 2. * DN;  // 対角成分
                coef[i][i + 1] = -(DN-0.5*CFL);       // 右隣の成分
            }
        }
        std::cout << "Heat1DImplicit solver initialized with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    // 時間ステップを実行するメソッド
    void Step() override
    {
        auto tmp=coef;
        for (int i = 0; i < Data.size(); i++)
        {
            if (i==0)                       tmp[i][Data.size()] = LBC_value; // 下端の境界条件
            else if (i == Data.size() - 1)  tmp[i][Data.size()] = RBC_value; // 上端の境界条件
            else                            tmp[i][Data.size()] = Data[i][nt]; // 現在の値を設定
        }
        tmp.MakeIndex(); // インデックスを作成
        tmp.EraseUp();   // 上三角行列に変形
        tmp.EraseDown(); // 下三角行列に変形
        tmp.ClearIndex(); // インデックスをクリア
        time = (++nt) * dt;
        for (int i = 0; i < Data.size(); i++)
        {
            Data[i][nt] = tmp[i][Data.size()];
        }
    }
};
// heat1D_implicit.hpp