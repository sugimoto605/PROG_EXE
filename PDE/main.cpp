// main.cpp
#include <iostream>
#include "pCONV.hpp"
#include "pSolverLD.hpp"
#include "pSolverHD.hpp"
#include "pSolveriS.hpp"
#include "pSolverAViS.hpp"
#include "pSolverFVLW.hpp"
#include "pSolverFVMinMod.hpp"  //MinMod-Limitter
#include "pSolverFVVL.hpp"      //van-Leer Limitter
int main()
{
    size_t nx = 20;          // 空間分割数 100(text), 40, 20, 8
    double dt = 0.8 / nx;    // 時間刻み幅 CFL=0.8 に固定しとこ. つまり 0.2/40
    size_t ntmax = 5. / dt;  // 最大時間ステップ数
    size_t ntsave = 5. / dt; // データ保存間隔
    auto I_SIN = [](double x)
    { return std::pow(std::sin(2.0 * M_PI * x), 5); }; // 初期条件として正弦波^5を設定
    auto I_LEVEQUE = [](double x)
    {
        double b=200;
        double xc=x-std::floor(x+0.2)-0.3;
        double y=exp(-b*xc*xc);
        if ((xc>0.3)&&(xc<0.5)) y+=1;
        return y;
    }; // Lebequeのテキストの初期関数
    auto I_DC = [](double x)
    { return (x > 0.8) ? 0 : std::max(2. * x - 0.6, 0.); };
    //--------------------------------------------------
    std::string filename = "Data/test03/TEXT_FVMM20.dat";
    pSolverFVMinMod mySolver(nx, dt, I_LEVEQUE);
    //--------------------------------------------------
    mySolver.Initialize();
    for (size_t nt = 0; nt < ntmax; nt++)
    {
        mySolver.Step(); // 時間ステップを実行
        if ((nt + 1) % ntsave == 0)
            mySolver.Write(filename);
    }
    return 0;
}

// mySolver.StepLQ();       // エクスプリシットでの時間ステップ
// mySolver.StepHQ();       // エクスプリシットでの時間ステップ
//  mySolver.Step_Implicit(); // インプリシットでの時間ステップ
// #else
//     std::cout << "This is a test build. No PDE solver executed." << std::endl;
//     pvector<double> Data(7);
//     std::cout << "SIZE: " << Data.size() << std::endl;
//     for (size_t n = 0; auto &v : Data)
//         v = n++;
//     for (int i = -10; i < 10; i++)
//         std::cout << "Data[" << i << "] = " << Data[i] << std::endl;
// #endif
