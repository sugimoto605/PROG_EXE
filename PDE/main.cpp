// main.cpp
#include <iostream>
#include "pCONV.hpp"
#include "pSolverLD.hpp"
#include "pSolverHD.hpp"
#include "pSolverSi.hpp"
int main()
{
#ifndef _TEST
    size_t nx=40; // 空間分割数
    double dt=0.005; // 時間刻み幅
    size_t ntmax=200; // 最大時間ステップ数
    size_t ntsave=200; // データ保存間隔
#ifdef _OLD
    pCONV mySolver(nx, dt, [](double x)
                   { return std::pow(std::sin(2.0 * M_PI * x), 5); }); // 初期条件として正弦波を設定
#else
    pSolverLD mySolver(nx, dt, [](double x)
                       { return std::pow(std::sin(2.0 * M_PI * x), 5); }); // 初期条件として正弦波を設定
#endif
    mySolver.Write("Data/test03/PSOLVER_LD.dat"); // 初期状態の保存
    for(size_t nt=0;nt<ntmax;nt++)
    {
        mySolver.Step(); // 時間ステップを実行
        if ((nt+1) % ntsave == 0) mySolver.Write();
    }
#else
    std::cout << "This is a test build. No PDE solver executed." << std::endl;
    pvector<double> Data(7);
    std::cout << "SIZE: " << Data.size() << std::endl;
    for (size_t n=0;auto &v : Data) v=n++;
    for(int i=-10;i<10;i++) std::cout << "Data[" << i << "] = " << Data[i] << std::endl;
#endif
    return 0;
}

// mySolver.StepLQ();       // エクスプリシットでの時間ステップ
// mySolver.StepHQ();       // エクスプリシットでの時間ステップ
//  mySolver.Step_Implicit(); // インプリシットでの時間ステップ
