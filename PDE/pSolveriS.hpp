// pSolveriS.hpp
//  implicit-Sakurai Scheme: Akira Sakurai, private communication, 1992
//
#pragma once
#include "../include/pSolver1DBase.hpp"
// 東京電機大学の桜井先生のやつ
// データは(U,∂U/∂t,∂U/∂x) U, Ut, Ux
//   x=0     u=U0, u'=dU0   U=ax^2+dU0*x+U0
//   |       |
//   x=D     u=U1, u'=dU1   U1=a*D^2+dU0*D+U0 a=(U1-U0)/D^2-dU0/D dU1=2*(U1-U0)/D-dU0
//
//  Implicit-Sakurai
//  A       B         Bにおける ∂u/∂x = 2*(B.U-A.U)/dX-A.Ux
//  *-------*         Bにおける ∂u/∂t = 2*(B.U-C.U)/dt-C.Ut
//          |         2*(B.U-C.U)/dt-C.Ut + C*[2*(B.U-A.U)/dX-A.Ux] = 0
//          |         B.U-C.U-dt*C.Ut/2 + CFL*(B.U-A.U)-CFL*dX*A.Ux/2 = 0
//          * C       B.U+ CFL*B.U = C.U+dt*C.Ut/2+CFL(*A.U+dX*A.Ux/2)
//                    B.U = [C.U+dt*C.Ut/2+CFL(A.U+dX*A.Ux/2)]/(1+CFL)

class SKit
{ // 桜井明キット
public:
    double U;
    double Ut;
    double Ux;
    SKit &operator=(const double &x)
    {
        U = x;
        Ut = Ux = 0;
        return *this;
    }
    operator double() const { return U; }
};
class pSolveriS : public pSolver1DBase<SKit>
{
public:
    pSolveriS(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                  { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverSi..." << std::endl;
    }
    void Initialize(void *parm) override
    {
        std::cout << "Prepareing initial values..." << std::endl;
        // 桜井法は, Ux, Ut が初期値で必要である
        for (auto &R : Data)
        {
            auto &P = *R.pre;
            auto &L = *P.pre;
            //  L   P   R
            P[nt].Ux = (R[nt].U - L[nt].U) * .5 / dx;
            P[nt].Ut = -C * P[nt].Ux;
        }
    }
    //陰的櫻井法
    void Step() override
    {
        //周期条件の場合, Dataの末尾がまるで境界条件
        Data.back()[nt+1]=5*Data.back()[nt];
        //解く
        for (size_t ite = 0; ite < 2; ite++)
            for (auto &P : Data)
            {
                auto &Q = *P.pre;
                double Ct = P[nt].U + .5 * dt * P[nt].Ut;
                double Cx = Q[nt + 1].U + .5 * dx * Q[nt + 1].Ux;
                P[nt + 1].U = (Ct + CFL * Cx) / (1. + CFL);
                P[nt + 1].Ux = 2. * (P[nt + 1].U - Q[nt + 1].U) / dx - Q[nt + 1].Ux;
                P[nt + 1].Ut = 2. * (P[nt + 1].U - P[nt].U) / dt - P[nt].Ut;
            }
        time = (++nt) * dt;
    }
};