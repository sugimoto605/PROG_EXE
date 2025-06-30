// AD_P1D.hpp
//  Advection-Diffusion Periodic1D solver   t-2nd x-2nd Finite Difference
#pragma once
#include "../include/pSolver1DBase.hpp"
#include "../include/libSPARSE.hpp"
class AD_P1D : public pSolver1DBase<double>
{
    SPARSE::Coef coef;
    double Usum; // Uの合計値
public:
    AD_P1D(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                               { return 0.; }) : pSolver1DBase<double>(nx, dt)
    {
        I_func= U0; // 初期条件の関数を設定
        std::cout << "Starting AdvectionDiffusion Periodic1D solver with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        if (parm)                             // パラメータが与えられているなら Kを1.0から変更するつもりであろう
            K = *static_cast<double *>(parm); // 拡散係数の初期化
        DN = K * dt / dx / dx;                // 拡散数の計算
        for (int i = 0; i < Data.size(); i++)
            coef.New(i, 0.0);
        Usum = 0.0; // Uの合計値を初期化
        for (int i = 0; i < Data.size(); i++)
        {
            Usum += Data[i][0]; // Uの合計値を計算
            if (i < Data.size() - 1)
            {
                coef[i][i] = -(DN +CFL*.5); // 左隣の成分
                coef[i][i+1] = 2. + DN*2.;      // 対角成分
                if (i < Data.size() - 2)
                    coef[i][i + 2] = -(DN - CFL*.5); // 右隣の成分
                else
                    coef[i][0] = -(DN - CFL*.5); // 周期境界条件
            }
            else // 最後の行
                for (int j = 0; j < Data.size(); j++)
                    coef[i][j] = 1.0; // 総量保存
        }
        std::cout << "Heat1DImplicit solver initialized with nx=" << Data.size() << ", dt=" << dt << std::endl;
    }
    // 時間ステップを実行するメソッド
    void Step() override
    {
        auto tmp = coef;
        for (int i = 0; i < Data.size(); i++)
        {
            if (i < Data.size() - 1)
                tmp[i][Data.size()] = (DN+CFL*.5)*Data[i][nt]+(2.-DN*2.)*Data[i+1][nt]+(DN-CFL*.5)*Data[i+2][nt]; // 現在の値を設定
            else
                tmp[i][Data.size()] = Usum; // 最後の行は総量保存のために合計値を設定
        }
        tmp.MakeIndex();  // インデックスを作成
        tmp.EraseUp();    // 上三角行列に変形
        tmp.EraseDown();  // 下三角行列に変形
        tmp.ClearIndex(); // インデックスをクリア
        time = (++nt) * dt;
        for (int i = 0; i < Data.size(); i++)
        {
            Data[i][nt] = tmp[i][Data.size()];
        }
    }
};
// heat1DP_implicit.hpps