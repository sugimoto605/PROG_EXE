// main.cpp
#include <iostream>
#include "CONV.hpp"
#include "pvector.hpp"
#define _TEST
int main()
{
#ifndef _TEST
    size_t nx=10; // 空間分割数
    double dt=0.1; // 時間刻み幅
    size_t ntmax=50; // 最大時間ステップ数
    size_t ntsave=10; // データ保存間隔
    CONV mySolver(nx,dt);
    mySolver.Write("Data/test03/PDE_02.dat"); // 初期状態の保存
    for(size_t nt=0;nt<ntmax;nt++)
    {
        // mySolver.Step();       // エクスプリシットでの時間ステップ
        mySolver.Step_Implicit(); // インプリシットでの時間ステップ
        if ((nt+1) % ntsave == 0) mySolver.Write();
    }
#else
    std::cout << "This is a test build. No PDE solver executed." << std::endl;
    pvector Data(7);
    std::cout << "SIZE: " << Data.size() << std::endl;
    for (size_t n=0;auto &v : Data) v=n++;
    for(int i=-10;i<10;i++) std::cout << "Data[" << i << "] = " << Data[i] << std::endl;
#endif
    return 0;
}