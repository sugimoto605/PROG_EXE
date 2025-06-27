// heat1D_sakurai.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

class SKit
{
public:
    double U;     // U(t,x)
    double dUdx;  // ∂u/∂x
    double dUdx2; // ∂²u/∂x²
    double dUdx3; // ∂³u/∂x³
    double dUdt;  // ∂u/∂t
    double dUdt2; // ∂²u/∂t²
    SKit &operator=(const double &x)
    {
        U = x;
        dUdx = dUdx2 = dUdx3 = dUdt = dUdt2 = 0;
        return *this;
    }
    operator double() const { return U; }
};

class Heat1D_Sakurai : public bSolver1DBase<SKit>
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> coef; // スパース行列を使用して係数行列を表現
public:
    Heat1D_Sakurai(const size_t &nx, const double &dt) : bSolver1DBase<SKit>(nx, dt)
    {
        std::cout << "Booting Heat1D_Sakurai" << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        if (!dLBC_func)
            throw std::runtime_error("dLBC_func (du/dt by B.C.) is not set.");
        if (!dRBC_func)
            throw std::runtime_error("dRBC_func (du/dt by B.C.) is not set.");
        if (!dI_func)
            throw std::runtime_error("dI_func (du/dt calculated from I.C.) is not set.");
        if (parm)
        {
            auto P = *static_cast<std::pair<double, double> *>(parm);
            C = P.first;  // 速度の初期化
            K = P.second; // 拡散係数の初期化
        }
        bSolver1DBase<SKit>::Initialize(parm); // 基底クラスの初期化
        for (auto &P : Data)
            P[nt].dUdt = dI_func(P.X);    // 初期条件の追加
        coef.resize(6 * Data.size() - 5, 6 * Data.size() - 5);
        size_t idx = 0;
        for (int i = 0; i < Data.size(); i++) // Data.size()=nx+1
        {
            if (i == 0)
            {
                coef.insert(idx, idx) = 1.; // LBC
                idx++;
                coef.insert(idx, idx - 1) = 1.;                 // *.1
                coef.insert(idx, idx) = +dx;                    // *.1
                coef.insert(idx, idx + 1) = +dx * dx / 2.;      // *.1
                coef.insert(idx, idx + 2) = +dx * dx * dx / 6.; // *.1
                coef.insert(idx, idx + 3) = -1;                 // *.1
                idx++;
                coef.insert(idx, idx - 1) = 1.;            // *.2
                coef.insert(idx, idx) = +dx;               // *.2
                coef.insert(idx, idx + 1) = +dx * dx / 2.; // *.2
                coef.insert(idx, idx + 3) = -1;            // *.2
                idx++;
                coef.insert(idx, idx - 1) = 1.; // *.3
                coef.insert(idx, idx) = +dx;    // *.3
                coef.insert(idx, idx + 3) = -1; // *.3
                idx++;
                coef.insert(idx, idx - 4) = 0.;     // alpha/gamma of 0 = alpha u + beta ∂u/∂x + gamma ∂²u/∂x² - u_t
                coef.insert(idx, idx - 3) = -C / K; // beta/gamma
                coef.insert(idx, idx - 2) = +1.;    // gamma/gamma
                idx++;
            }
            else if (i < Data.size() - 1)
            {
                coef.insert(idx, idx - 1) = 1.;                 // *.1
                coef.insert(idx, idx) = +dx;                    // *.1
                coef.insert(idx, idx + 1) = +dx * dx / 2.;      // *.1
                coef.insert(idx, idx + 2) = +dx * dx * dx / 6.; // *.1
                coef.insert(idx, idx + 5) = -1;                 // *.1
                idx++;
                coef.insert(idx, idx - 1) = 1.;            // *.2
                coef.insert(idx, idx) = +dx;               // *.2
                coef.insert(idx, idx + 1) = +dx * dx / 2.; // *.2
                coef.insert(idx, idx + 5) = -1;            // *.2
                idx++;
                coef.insert(idx, idx - 1) = 1.; // *.3
                coef.insert(idx, idx) = +dx;    // *.3
                coef.insert(idx, idx + 5) = -1; // *.3
                idx++;
                coef.insert(idx, idx - 4) = 1.;                 // +.1
                coef.insert(idx, idx    ) = -dt;                // +.1
                coef.insert(idx, idx + 1) = +dt * dt / 2.;      // +.1
                idx++;
                coef.insert(idx, idx - 1) = 1.;                 // +.2
                coef.insert(idx, idx    ) = -dt;                // +.2
                idx++;
                coef.insert(idx, idx - 6) = 0.;     // alpha/gamma of 0 = alpha u + beta ∂u/∂x + gamma ∂²u/∂x² - u_t
                coef.insert(idx, idx - 5) = -C / K; // beta/gamma
                coef.insert(idx, idx - 4) = +1.;    // 1
                coef.insert(idx, idx - 2) = -1./ K; // -1/gamma
                idx++;
            }
            else
            {
                coef.insert(idx, idx - 1) = 1.; // RBC
                idx++;
                coef.insert(idx, idx - 2) = 0.;     // alpha/gamma of 0 = alpha u + beta ∂u/∂x + gamma ∂²u/∂x² - u_t
                coef.insert(idx, idx - 1) = -C / K; // beta/gamma
                coef.insert(idx, idx) = +1.;        // gamma/gamma
            }
        }
        // 条件数とかみてみる
        Eigen::MatrixXd dense = Eigen::MatrixXd(coef);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
        double cond = svd.singularValues()(0) / svd.singularValues().tail(1)(0);
        std::cout << "Condition number: " << cond << std::endl;
        std::cout << "SVD rank: " << svd.rank() << " of " << dense.rows() << " x " << dense.cols() << std::endl;
        std::cout << "Heat1D_Sakurai det= " << dense.determinant() << std::endl;
        // Eigen::IOFormat CleanFmt(2, 0, ", ", "\n", "[", "]");
        // std::cout << dense.format(CleanFmt) << std::endl;
        coef.makeCompressed(); // スパース行列を圧縮
        solver.compute(coef);
        std::cout << "Heat1DImplicit solver initialized with matrix " << solver.rows() << " x " << solver.cols() << ", dt=" << dt << std::endl;
    }
    void Step() override
    {
        // 時間ステップを実行
        Eigen::VectorXd b(solver.rows());
        {
            size_t idx = 0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                if (i == 0) // 下端の境界条件
                {
                    b[idx++] = LBC_func(time); // 現在の値を設定
                    b[idx++] = 0.;
                    b[idx++] = 0.;
                    b[idx++] = 0.;
                    b[idx++] = dLBC_func(time) / K;
                }
                else if (i < Data.size() - 1)
                {
                    b[idx++] = 0.;
                    b[idx++] = 0.;
                    b[idx++] = 0.;
                    b[idx++] = Data[i][nt].U;
                    b[idx++] = Data[i][nt].dUdt;
                    b[idx++] = 0.;
                }
                else
                {
                    b[idx++] = RBC_func(time); // 上端の境界条件
                    b[idx++] = dRBC_func(time) / K;
                }
            }
        }
        auto result = solver.solve(b); // 連立方程式を解く
        time = (++nt) * dt;
        {
            size_t idx=0;
            for (int i = 0; i < Data.size(); i++)
            {
                if (i == 0)
                {
                    Data[i][nt].U = result[idx++];
                    Data[i][nt].dUdx = result[idx++];
                    Data[i][nt].dUdx2 = result[idx++];
                    Data[i][nt].dUdx3 = result[idx++];
                }
                else if (i < Data.size() - 1)
                {
                    Data[i][nt].U = result[idx++];
                    Data[i][nt].dUdx = result[idx++];
                    Data[i][nt].dUdx2 = result[idx++];
                    Data[i][nt].dUdx3 = result[idx++];
                    Data[i][nt].dUdt = result[idx++];
                    Data[i][nt].dUdt2 = result[idx++];  
                }
                else
                {
                    Data[i][nt].U = result[idx++];
                    Data[i][nt].dUdx = result[idx++];
                    Data[i][nt].dUdx2 = result[idx++];
                }
            }
        }
        // デバッグ用の出力
        std::cout << "STEP " << nt << " at time " << time << std::endl;
        std::cout << std::setw(12) << "X" 
            << std::setw(12) << "U"
            << std::setw(12) << "U_x"
            << std::setw(12) << "U_xx" 
            << std::setw(12) << "U_xxx"
            << std::setw(12) << "U_t"
            << std::setw(12) << "U_tt"
            << std::endl;
        for(size_t i = 0; i < Data.size(); ++i)
        {
            std::cout << std::fixed << std::setprecision(6)
                      << std::setw(12) << Data[i].X
                      << std::setw(12) << Data[i][nt].U
                      << std::setw(12) << Data[i][nt].dUdx
                      << std::setw(12) << Data[i][nt].dUdx2
                      << std::setw(12) << Data[i][nt].dUdx3
                      << std::setw(12) << Data[i][nt].dUdt
                      << std::setw(12) << Data[i][nt].dUdt2
                      << std::endl;
        }   
    };
    std::string what() const override
    {
        return "Heat1D_Sakurai: " + bSolver1DBase<SKit>::what();
    }
};
