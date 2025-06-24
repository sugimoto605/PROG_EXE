// heat1D_sakurai.hpp
#pragma once
#include "../include/bSolver1DBase.hpp"
#include "../include/libSPARSE.hpp"

class SKit
{
public:
    double U;           // U(t,x)
    double dUdx;        // ∂u/∂x
    double dUdx2;      // ∂²u/∂x²
    double dUdx3;      // ∂³u/∂x³
    double dUdt;        // ∂u/∂t
    double dUdt2;       // ∂²u/∂t²
    SKit &operator=(const double &x)
    {
        U=x;
        dUdx=dUdx2=dUdx3=dUdt=dUdt2=0;
        return *this;
    }
    operator double() const { return U; }
};

class Heat1D_Sakurai : public bSolver1DBase<SKit>
{
    SPARSE::Coef coef;
public:
    Heat1D_Sakurai(const size_t &nx, const double &dt, const std::function<double(double)> &U0 
    = [](double x) { return 0.; }) : bSolver1DBase<SKit>(nx, dt, U0)
    {
        std::cout << "Booting Heat1D_Sakurai" << std::endl;
    }
    void Initialize(void *parm = nullptr) override
    {
        // 初期分布から係数を計算
        SPARSE::Coef iCoef;
        for (size_t i = 0; i < 4*Data.size()-3; ++i) iCoef.New(i, 0.0);
        size_t SZ=iCoef.Count();
        std::cout << "Initializing Coef with size " << SZ << " dx= " << dx << std::endl;
        // 行列設定
        {
            size_t idx = 0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                if (i == 0)
                {
                    iCoef[idx][idx] = dx;                           //OK
                    iCoef[idx][idx + 1] = dx * dx / 2.0;            //OK
                    iCoef[idx][idx + 2] = dx * dx * dx / 6.0;       //OK
                    idx++;
                    iCoef[idx][idx - 1] = 1.0;                      //OK     
                    iCoef[idx][idx] = dx;                           //OK   
                    iCoef[idx][idx + 1] = dx * dx / 2.0;            //OK
                    iCoef[idx][idx + 2] = -1.0;                     //OK  
                    idx++;
                    iCoef[idx][idx - 1] = 1.0;                      //OK     
                    iCoef[idx][idx] = dx;                           //OK  
                    iCoef[idx][idx + 2] = -1.0;                     //OK 
                    idx++;
                    iCoef[idx][idx - 3] = -C;
                    iCoef[idx][idx - 2] = +K;
                    idx++;
                }
                else if (i < Data.size() - 1)
                {
                    iCoef[idx][idx - 1] = dx;
                    iCoef[idx][idx] = dx * dx / 2.0;
                    iCoef[idx][idx + 1] = dx * dx * dx / 6.0;
                    idx++;
                    iCoef[idx][idx - 2] = 1.0;
                    iCoef[idx][idx - 1] = dx;
                    iCoef[idx][idx] = dx * dx / 2.0;
                    iCoef[idx][idx + 2] = -1.0;
                    idx++;
                    iCoef[idx][idx - 2] = 1.0;
                    iCoef[idx][idx - 1] = dx;
                    iCoef[idx][idx + 1] = -1.0;
                    idx++;
                    iCoef[idx][idx - 4] = -C;
                    iCoef[idx][idx - 3] = +K;
                    iCoef[idx][idx - 1] = -1.0;
                    idx++;
                }
                else
                {
                    iCoef[idx][idx - 1] = -C;
                    iCoef[idx][idx] = +K;
                    idx++;
                }
            }
            std::cout << "array idx= " << idx << std::endl;
        }
        std::cout << "iCoef Dimension=" << iCoef.Dimension() << " Count=" << iCoef.Count() << std::endl;
        iCoef.MakeIndex();
        iCoef.PrintIndex();
        std::cout << "CheckIndex ******" << std::endl;
        iCoef.CheckIndex();
        iCoef.Print();
        std::cout << "Swap 3-6" << std::endl;
        iCoef.Swap(3, 6);
        iCoef.Print();
        iCoef.CheckIndex();
        // ベクトル設定
        {
            size_t idx = 0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                if (i == 0)
                {
                    iCoef[idx][SZ] = Data[i + 1][nt].U - Data[i][nt].U;
                    idx += 3;
                    iCoef[idx][SZ] = 0.0; // 本当は ∂u/∂t(X=0)-alpha*u(X=0)
                    idx++;
                }
                else if (i < Data.size() - 1)
                {
                    iCoef[idx][SZ] = Data[i + 1][nt].U - Data[i][nt].U;
                    idx += 3;
                    iCoef[idx][SZ] = 0.0; // 本当は -alpha*u(X=i)
                    idx++;
                }
                else
                {
                    iCoef[idx][SZ] = 0.0; // 本当は ∂u/∂t(X=1)-alpha*u(X=1)
                    idx++;
                }
            }
            std::cout << "vector idx= " << idx << std::endl;
        }
        std::cout << "iCoef Dimension=" << iCoef.Dimension() << " Count=" << iCoef.Count() << std::endl;
        iCoef.Print();
        // 係数決定
        auto jCoef=iCoef;
        iCoef.MakeIndex();
        std::cout << "MakeIndex=" << std::endl;
        iCoef.PrintIndex();
        if (!iCoef.EraseDown()) throw std::runtime_error("Failed to erase down.");
        std::cout << "EraseDown=" << std::endl;
        iCoef.Print();
        iCoef.EraseUp();
        iCoef.RemoveIndex();
        std::cout << "Solved=" << std::endl;
        iCoef.Print();
        {
            size_t idx=0;
            for (size_t i = 0; i < Data.size(); ++i)
            {
                if (i==0)
                {
                    Data[i][nt].dUdx= iCoef[idx++][SZ];
                    Data[i][nt].dUdx2 = iCoef[idx++][SZ];
                    Data[i][nt].dUdx3 = iCoef[idx++][SZ];
                }
                else if (i < Data.size() - 1)
                {
                    Data[i][nt].dUdx = iCoef[idx++][SZ];
                    Data[i][nt].dUdx2 = iCoef[idx++][SZ];
                    Data[i][nt].dUdx3 = iCoef[idx++][SZ];
                    Data[i][nt].dUdt = iCoef[idx++][SZ];
                }
                else
                {
                    Data[i][nt].dUdx = iCoef[idx++][SZ];
                    Data[i][nt].dUdx2 = iCoef[idx++][SZ];
                }
            }
            std::cout << "copy idx= " << idx << std::endl;
        }
        // 検算
        for(int i = Data.size()-1; i >=0; --i)
        {
            auto &P = Data[i];
            std::cout << "Point[" << i << "]: X=" << std::fixed << P.X 
                      << ", U=" << P[nt].U 
                      << ", dUdx=" << P[nt].dUdx 
                      << ", dUdx2=" << P[nt].dUdx2 
                      << ", dUdx3=" << P[nt].dUdx3 
                      << ", dUdt=" << P[nt].dUdt 
                      << std::endl;
            if (i == Data.size()-1)
            {
                //最後の行をチェック
                size_t idx=4*i;
                double S=jCoef[idx][idx-1]*Data[i][nt].dUdx+jCoef[idx][idx]*Data[i][nt].dUdx2-jCoef[idx][SZ];;
                std::cout << "Check[" << idx << "]S= " << S <<std::endl;
            }
            else
            {
                size_t idx = 4 * i+3;
                double S = jCoef[idx][idx - 4] * Data[i][nt].dUdx + jCoef[idx][idx-3] * Data[i][nt].dUdx2 + jCoef[idx][idx-1]*Data[i][nt].dUdt - jCoef[idx][SZ];
                double Y = jCoef[idx][idx - 4] * iCoef[idx-4][SZ] + jCoef[idx][idx - 3] * iCoef[idx-3][SZ] + jCoef[idx][idx - 1] * iCoef[idx-1][SZ] - jCoef[idx][SZ];
                std::cout << "Check[" << idx << "]S=" << S << " Y= " << Y << std::endl;
                idx--;
                std::cout << "Check[" << idx << "]" << std::endl;
                idx--;
                std::cout << "Check[" << idx << "]" << std::endl;
                idx--;
                std::cout << "Check[" << idx << "]" << std::endl;
            }
        }
    }
    void Step() override
    {
        // 時間ステップを実行
        // coef.MakeIndex();
        // coef.EraseDown();
        // coef.EraseUp();
        // coef.RemoveIndex();
        // coef.MakeIndex();
        // for (size_t i = 0; i < Data.size(); ++i)
        // {
        //     auto &P = Data[i];
        //     P[nt + 1].U = coef.Solve(P[nt].U, P[nt].dUdx, P[nt].dUdx2, P[nt].dUdt);
        //     P[nt + 1].dUdt = (P[nt + 1].U - P[nt].U) / dt;
        //     P[nt + 1].dUdx = (P[nt + 1].U - P[nt - 1].U) / (2 * dx);
        //     P[nt + 1].dUdx2 = (P[nt + 1].U - 2 * P[nt].U + P[nt - 1].U) / (dx * dx);
        //     P[nt + 1].dUdx3 = (P[nt + 1].dUdx2 - P[nt - 1].dUdx2) / (2 * dx);
        // }
        nt++;
        time += dt;
    };
    std::string what() const override
    {
        return "Heat1D_Sakurai: " + bSolver1DBase<SKit>::what();
    }
};
