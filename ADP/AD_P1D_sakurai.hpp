// AD_P1D_sakurai.hpp
//  Advection-Diffusion Periodic1D solver Sakurai Scheme
#undef CONSERVE
#pragma once
#include "../include/pSolver1DBase.hpp"
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
class AD_P1D_Sakurai : public pSolver1DBase<SKit>
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
public:
    AD_P1D_Sakurai(const size_t &nx, const double &dt) : pSolver1DBase<SKit>(nx, dt)
    {
        std::cout << "Booting AD_P1D_sakurai" << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        if (!dI_func) throw std::runtime_error("AD_P1D_Sakurai::Initialize::Requires dI_Func.");
        if (parm)
        {
            auto P = *static_cast<std::pair<double, double> *>(parm);
            C = P.first;  // 速度の初期化
            K = P.second; // 拡散係数の初期化
        }
        for (auto &P : Data) P[nt].dUdt = dI_func(P.X);         // 初期条件の適用
        pSolver1DBase<SKit>::Initialize(parm); // 基底クラスの初期化
        Eigen::SparseMatrix<double> coef;      // スパース行列を使用して係数行列を表現
        coef.resize(6 * Data.size(), 6 * Data.size());
        size_t idx = 0;
        for (int i = 0; i < Data.size(); i++) // Data.size()=nx+1
        {
            // Line-1
            coef.insert(idx, idx) = 1.;
            coef.insert(idx, idx + 1) = dx;                // *.1
            coef.insert(idx, idx + 2) = dx * dx / 2.;      // *.1
            coef.insert(idx, idx + 3) = dx * dx * dx / 6.; // *.1
            if (i < Data.size() - 1)
                coef.insert(idx, idx + 6) = -1; // *.1
            else
                coef.insert(idx, 0) = -1;
            idx++;
            // Line-2
            coef.insert(idx, idx) = 1.;               // *.2
            coef.insert(idx, idx + 1) = dx;           // *.2
            coef.insert(idx, idx + 2) = dx * dx / 2.; // *.2
            if (i < Data.size() - 1)
                coef.insert(idx, idx + 6) = -1; // *.2
            else
                coef.insert(idx, 1) = -1;
            idx++;
            // Line-3
            coef.insert(idx, idx) = 1;      // *.3
            coef.insert(idx, idx + 1) = dx; // *.3
            if (i < Data.size() - 1)
                coef.insert(idx, idx + 6) = -1; // *.3
            else
                coef.insert(idx, 2) = -1;
            idx++;
            // Line-4
            coef.insert(idx, idx - 3) = 1.;            // +.1
            coef.insert(idx, idx + 1) = -dt;           // +.1
            coef.insert(idx, idx + 2) = +dt * dt / 2.; // +.1
            idx++;
            // Line-5
            coef.insert(idx, idx) = 1.;      // +.2
            coef.insert(idx, idx + 1) = -dt; // +.2
            idx++;
            // Line-6
#ifndef CONSERVE
            coef.insert(idx, idx - 5) = 0.;      // alpha
            coef.insert(idx, idx - 4) = -C / K;  // beta
            coef.insert(idx, idx - 3) = 1;       // gamma
            coef.insert(idx, idx - 1) = -1. / K; // -1/gamma
#else
            for(int j=0;j<Data.size()-1;j++) coef.insert(idx,6*j)=1.;
#endif
            idx++;
        }
        coef.makeCompressed(); // スパース行列を圧縮
        solver.compute(coef);
        if (solver.info()!=Eigen::Success)
            throw std::runtime_error("Error::Can not solve.");
        // // 条件数とかみてみる
        // Eigen::MatrixXd dense = Eigen::MatrixXd(coef);
        // std::cout << "Solver completed" << std::endl;
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // double cond = svd.singularValues()(0) / svd.singularValues().tail(1)(0);
        // std::cout << "Condition number: " << cond << std::endl;
        // std::cout << "SVD rank: " << svd.rank() << " of " << dense.rows() << " x " << dense.cols() << std::endl;
        // std::cout << "AD_P1D_Sakurai det= " << dense.determinant() << std::endl;
        // // std::cout << dense << std::endl;
    }
    void Step() override
    {
        Eigen::VectorXd b(solver.rows());
        {
            size_t idx = 0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                b[idx++] = 0.;
                b[idx++] = 0.;
                b[idx++] = 0.;
                b[idx++] = Data[i][nt].U;
                b[idx++] = Data[i][nt].dUdt;
                b[idx++] = 0.;
            }
#ifdef CONSERVE
            double sum=0;
            for(int j=0;j<Data.size();j++) sum+=Data[j][nt].U;
            b[b.size()-1]=sum;
#endif
        }
        auto result = solver.solve(b); // 連立方程式を解く
        time = (++nt) * dt;
        {
            size_t idx = 0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                Data[i][nt].U     = result[idx++];
                Data[i][nt].dUdx  = result[idx++];
                Data[i][nt].dUdx2 = result[idx++];
                Data[i][nt].dUdx3 = result[idx++];
                Data[i][nt].dUdt  = result[idx++];
                Data[i][nt].dUdt2 = result[idx++];
            }
        }
    }
};
