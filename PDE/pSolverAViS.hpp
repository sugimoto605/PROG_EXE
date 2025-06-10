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
//  ううん，これ，拡散方程式を学んでからでは？
//          CFLも，そっちのを満たすものじゃねえとやばい感じ
//  Implicit-Sakurai
//  A       B         Bにおける ∂u/∂x = 2*(B.U-A.U)/dX-A.Ux
//  *-------*         Bにおける ∂u/∂t = 2*(B.U-C.U)/dt-C.Ut
//          |         2*(B.U-C.U)/dt-C.Ut + C*[2*(B.U-A.U)/dX-A.Ux] = 0
//          |         B.U-C.U-dt*C.Ut/2 + CFL*(B.U-A.U)-CFL*dX*A.Ux/2 = 0
//          * C       B.U+ CFL*B.U = C.U+dt*C.Ut/2+CFL(*A.U+dX*A.Ux/2)
//                    B.U = [C.U+dt*C.Ut/2+CFL(A.U+dX*A.Ux/2)]/(1+CFL)
//  [人工粘性]基礎方程式を ∂U/∂t+c∂U/∂x=k ∂^2 U/∂x^2  にする.
//       k=eps*dx^n   (n=1,2,3,...)とでもしとけば，k→0 だからええやん，という考え方である
//       やってみたが，ますます発散してわけがわからないのであった.
//  データは((U,∂U/∂t,∂U/∂x,∂^2U/∂x^2 ) (U, Ut, Ux, Uxx)
//   x=0     U=U0, U'=d1U0, U"=d2U0   U=a*x^3+d2U0*x^2/2+d1U0*x+U0
//   |       |                        U'=3*a*x^2+d2U0*x+d1U0  U"=6*a*x+d2U0
//   |       |
//   x=D     U=U1=a*D^3+d2U0*D^2/2+d1U0*D+U0, a=(U1-U0)/D^3-d1U0/D^2-0.5*d2U0/D                       OK
//           U'=d1U1=3*a*D^2+d2U0*D+d1U0=3*[(U1-U0)/D^3-d1U0/D^2-0.5*d2U0/D]*D^2+d2U0*D+d1U0
//                  =3*(U1-U0)/D-2*d1U0 -0.5*d2U0*D                                                   OK
//           U"=d2U1=6*[(U1-U0)/D^3-d1U0/D^2-0.5*d2U0/D]*D+d2U0
//                  =6*(U1-U0)/D^2 -6*d1U0/D -2*d2U0                                                  OK
//  Implicit-Sakurai
//  A       B         Bにおける ∂u/∂t      = 2*(B.U-C.U)/dt-C.Ut                             OK
//  *-------*         Bにおける ∂u/∂x      = 3*(B.U-A.U)/dx  -2*A.Ux   -0.5*A.Uxx*dx         OK
//          |         Bにおける ∂^2 u/∂x^2 = 6*(B.U-A.U)/dx^2-6*A.Ux/dx-  2*A.Uxx            OK
//          |         Bにおいて ∂u/∂t+C*∂u/∂x=K*∂^2 u/∂x^2 を要求すれば
//          * C       2*(B.U-C.U)/dt-C.Ut   +C*[3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx]
//                        = K*[6*(B.U-A.U)/dx^2-6*A.Ux/dx              -2*A.Uxx]
// 2*(B.U-C.U)-C.Ut*dt+dt*C/dx*[3*(B.U-A.U)-2*A.Ux*dx-0.5*A.Uxx*dx^2]=dt*K/dx^2*[6*(B.U-A.U)-6*A.Ux*dx-2*A.Uxx*dx^2]
// 2*(B.U-C.U)-C.Ut*dt    +CFL*[3*(B.U-A.U)-2*A.Ux*dx-0.5*A.Uxx*dx^2]=VTS*[6*(B.U-A.U)-6*A.Ux*dx-2*A.Uxx*dx^2]
//      VTS=K*dt/dx^2 として
// B.U-C.U-0.5*C.Ut*dt+1.5*CFL*(B.U-A.U)-CFL*A.Ux*dx-0.25*CFL*A.Uxx*dx^2
//       =3*VTS*(B.U-A.U)-3*VTS*A.Ux*dx-VTS*A.Uxx*dx^2
// B.U =[ (C.U+0.5*C.Ut*dt)
//       +CFL*(1.5*A.U+  A.Ux*dx+0.25*A.Uxx*dx^2)
//       -VTS*(3  *A.U+3*A.Ux*dx+     A.Uxx*dx^2)]/(1+1.5*CFL-3*VTS)
// B.Ut =2*(B.U-C.U)/dt-C.Ut
// B.Ux =3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx
// B.Uxx=6*(B.U-A.U)/dx^2-6*A.Ux/dx-2*A.Uxx
class AVSKit
{ // 桜井明キット
public:
    double U;
    double Ut;
    double Ux;
    double Uxx;
    AVSKit &operator=(const double &x)
    {
        U = x;
        Ut = Ux = Uxx = 0;
        return *this;
    }
    operator double() const { return U; }
};
class pSolverAViS : public pSolver1DBase<AVSKit>
{
    double VTS; // 人工粘性係数(Viscous Time Scale)
public:
    pSolverAViS(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                    { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverAViS..." << std::endl;
        VTS=0.;
    }
    // Initialize(K) で, 人工粘性係数が K*dx*dx になる
    void Initialize(void *parm) override
    {
        VTS = *static_cast<double *>(parm);
        double KV=VTS*dx*dx;
        VTS = KV*dt/dx/dx;
        std::cout << "VTS= " << VTS << " KV= " << KV << std::endl;
        // 桜井法は, Ux, Ut が初期値で必要である
        for (auto &R : Data)
        {
            auto &P = *R.pre;
            auto &L = *P.pre;
            //  L   P   R
            P[nt].Uxx = (R[nt].U - 2. * P[nt].U + L[nt].U) / dx / dx;
            P[nt].Ux = (R[nt].U - L[nt].U) * .5 / dx;
            P[nt].Ut = -C * P[nt].Ux+KV*P[nt].Uxx;
        }
    }
    // AV淫的櫻井法
    void Step(double U_W = 0.0) override
    {
        for (size_t ite = 0; ite < 2; ite++)
        {
            for (auto &P : Data)
            {
                std::cout << std::fixed << "time " << nt << " x= " << P.X << std::endl;
                std::cout << "[nt  ]  Ut [ " << std::setw(12) << P[nt].Ut << " ] ";
                std::cout << "  Ux [ " << std::setw(12) << P[nt].Ux << " ] ";
                std::cout << "  Uxx[ " << std::setw(12) << P[nt].Uxx << " ] ";
                std::cout << "  U  [ " << std::setw(12) << P[nt].U << " ] ";
                std::cout << std::endl;
                auto &Q = *P.pre;
                double Ct = P[nt].U + .5 * dt * P[nt].Ut;
                double Cx = 1.5 * Q[nt + 1].U +      Q[nt + 1].Ux * dx + 0.25 * Q[nt + 1].Uxx * dx * dx;
                double Cv = 3.0 * Q[nt + 1].U + 3. * Q[nt + 1].Ux * dx + Q[nt + 1].Uxx * dx * dx;
                // B.U =[ (C.U+0.5*C.Ut*dt)
                //       +CFL*(1.5*A.U+  A.Ux*dx+0.25*A.Uxx*dx^2)
                //       -VTS*(3  *A.U+3*A.Ux*dx+     A.Uxx*dx^2)]/(1+1.5*CFL-3*VTS)
                P[nt + 1].U = (Ct + CFL * Cx - VTS * Cv) / (1. + 1.5 * CFL - 3. * VTS);
                // B.Ut =2*(B.U-C.U)/dt-C.Ut
                // B.Ux =3*(B.U-A.U)/dx-2*A.Ux-0.5*A.Uxx*dx
                // B.Uxx=6*(B.U-A.U)/dx^2-6*A.Ux/dx-2*A.Uxx
                P[nt + 1].Ut  = 2. * (P[nt + 1].U - P[nt].U) / dt - P[nt].Ut;
                P[nt + 1].Ux  = 3. * (P[nt + 1].U - Q[nt + 1].U) / dx      - 2. * Q[nt + 1].Ux     - 0.5 * Q[nt + 1].Uxx *dx;
                P[nt + 1].Uxx = 6. * (P[nt + 1].U - Q[nt + 1].U) / dx / dx - 6. * Q[nt + 1].Ux / dx - 2. * Q[nt + 1].Uxx;
                std::cout << "[nt+1]  Ut [ " << std::setw(12) << P[nt+1].Ut << " ] ";
                std::cout << "  Ux [ " << std::setw(12) << P[nt+1].Ux << " ] ";
                std::cout << "  Uxx[ " << std::setw(12) << P[nt+1].Uxx << " ] ";
                std::cout << "  U  [ " << std::setw(12) << P[nt+1].U << " ] ";
                std::cout << std::endl;
            }
            throw(std::runtime_error("AHO"));
        }
        time = (++nt) * dt;
        throw(std::runtime_error("AHO"));
    }
};