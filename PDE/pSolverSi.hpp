// pSolverSi.hpp
#pragma once
#include "pSolver1DBase.hpp"
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
{//櫻井明キット
public:
    double U;
    double Ut;
    double Ux;
    SKit& operator=(const double& x){
        U=x;Ut=Ux=0;
        return *this;
    }
    operator double() const { return U;}
};
class pSolverSi : public pSolver1DBase<SKit>
{
public:
    pSolverSi(const size_t &nx, const double &dt, const std::function<double(double)> &U0 = [](double x)
                                                  { return 0.; }) : pSolver1DBase(nx, dt, U0)
    {
        std::cout << "Booting pSolverSi..." << std::endl;
    }
    //隠的櫻井法
    void Step(double U_W = 0.0) override
    {
        for(size_t ite=0;ite<2;ite++)
            for (auto &P : Data)
            {
                auto& Q=*P.pre;
                double Ct=P[nt].U+.5*dt*P[nt].Ut;
                double Cx=Q[nt+1].U+.5*dx*Q[nt+1].Ux;
                P[nt+1].U=(Ct+CFL*Cx)/(1.+CFL);
                P[nt+1].Ux=2.*(P[nt+1].U-Q[nt+1].U)/dx-Q[nt+1].Ux;
                P[nt+1].Ut=2.*(P[nt+1].U-P[nt].U)/dt-P[nt].Ut;
            }
        time = (++nt) * dt;
    }
};