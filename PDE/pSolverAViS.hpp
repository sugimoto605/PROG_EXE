// pSolverAViS.hpp
//  implicit-Sakurai Scheme: Akira Sakurai, private communication, 1992
//      with Artificial Viscosity
#pragma once
#include "pSolver1DBase.hpp"
// 東京電機大学の桜井先生のやつ
// データは(U,∂U/∂t,∂U/∂x) U, Ut, Ux
//   t=0     u=U0, u'=dU0   U=a*t^2+dU0*t+U0
//   |       |
//   t=D     u=U1=a*D^2+dU0*D+U0, a=(U1-U0)/D^2-dU0/D, u'=dU1=2*(U1-U0)/D-dU0
//
//  Implicit-Sakurai
//  A       B         Bにおける ∂u/∂x = 2*(B.U-A.U)/dX-A.Ux
//  *-------*         Bにおける ∂u/∂t = 2*(B.U-C.U)/dt-C.Ut
//          |         2*(B.U-C.U)/dt-C.Ut + C*[2*(B.U-A.U)/dX-A.Ux] = 0
//          |         B.U-C.U-dt*C.Ut/2 + CFL*(B.U-A.U)-CFL*dX*A.Ux/2 = 0
//          * C       B.U+ CFL*B.U = C.U+dt*C.Ut/2+CFL(*A.U+dX*A.Ux/2)
//                    B.U = [C.U+dt*C.Ut/2+CFL(A.U+dX*A.Ux/2)]/(1+CFL)
//  [人工粘性]基礎方程式を ∂U/∂t+c∂U/∂x=k ∂^2 U/∂x^2  にする.
//       k=eps*dx^n   (n=1,2,3,...)とでもしとけば，k→0 だからええやん，という考え方である
//  データは((U,∂U/∂t,∂U/∂x,∂^2U/∂x^2 ) (U, Ut, Ux, Uxx)
//   x=0     U=U0, U'=d1U0, U"=d2U0   U=a*x^3+d2U0*x^2/2+d1U0*x+U0
//   |       |
//   x=D     U=U1=a*D^3+d2U0*D^2/2+d1U0*D+U0, a=(U1-U0)/D^3-0.5*d2U0/D-d1U0/D^2
//           U'=d1U1=3*a*D^2+d2U0*D+d1U0=3*[(U1-U0)/D^3-0.5*d2U0/D-d1U0/D^2]*D^2+d2U0*D+d1U0
//                  =3*(U1-U0)/D-2*d1U0-0.5*d2U0*D
//           U"=d2U1=6*(U1-U0)/D^2-6*d1U0/D-2*d2U0
//  Implicit-Sakurai
//  A       B         Bにおける ∂u/∂t      = 2*(B.U-C.U)/dt-C.Ut
//  *-------*         Bにおける ∂u/∂x      = 3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx
//          |         Bにおける ∂^2 u/∂x^2 = 6*(B.U-A.U)/dx^2-6*A.Ux/dx-2*A.Uxx
//          |         Bにおいて ∂u/∂t+C*∂u/∂x=K*∂^2 u/∂x^2 を要求すれば
//          * C       2*(B.U-C.U)/dt-C.Ut+C*[3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx]
//                    = K*[6*(B.U-A.U)/dx^2-6*A.Ux/dx-2*A.Uxx]
// 2*(B.U-C.U)-C.Ut*dt+CFL*[3*(B.U-A.U)-2*A.Ux*dx-0.5*A.Uxx*dx^2]=K*dt/dx^2*[6*(B.U-A.U)-6*A.Ux*dx-2*A.Uxx*dx^2]
//      VTS=K*dt/dx^2 として
// B.U-C.U-0.5*C.Ut*dt+0.5*CFL*3*(B.U-A.U)-0.5*CFL*2*A.Ux*dx-0.5*CFL*0.5*A.Uxx*dx^2
//    =0.5*VTS*6*(B.U-A.U)-0.5*VTS*6*A.Ux*dx-0.5*VTS*2*A.Uxx*dx^2
// B.U =[(C.U+0.5*C.Ut*dt
//       +CFL*(1.5*A.U+  A.Ux*dx+0.25*A.Uxx*dx^2)
//       -VTS*(3  *A.U+3*A.Ux*dx+     A.Uxx*dx^2)]/(1+1.5*CFL-3*VTS)
// B.Ut =2*(B.U-C.U)/dt-C.Ut
// B.Ux =3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx
// B.Uxx=6*(B.U-A.U)/dx^2-6*A.Ux/dx-2*A.Uxx
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
class pSolverAViS : public pSolver1DBase<SKit>
{
    double AV;  // 人工粘性係数.
public:
    double kAV; //  AV=kAV*dx^2
    pSolverAViS(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                  { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverSi..." << std::endl;
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
    // AV淫的櫻井法

    // 陰的櫻井法
    void Step_org(double U_W = 0.0)
    {
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